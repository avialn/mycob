version 1.0

task Minimap2 {

    input {
        File? fastq_1
        File? fastq_2
        File ref
        Int threads
        String docker
    }

    command <<<
        set -ex -o pipefail

        minimap2 -t ~{threads} \
        -ax sr -c \
        --secondary=yes \
        ~{ref} \
        ~{fastq_1} ~{fastq_2} > file.sam

        matched_count=$(grep -v '^@' file.sam | cut -f 3 | sort | uniq -c | awk '$2 != "*" {print $1}')

        if [ -z "$matched_count" ]; then
            echo "No aligned reads" > proceed.txt
        else
            echo "yes" > proceed.txt
        fi

    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File file_sam = "file.sam"
        String proceed = read_string("proceed.txt")
    }
}

task Samtools {

    input {
        File ref
        File sam
        Int min_cov=10
        String docker

    }

    command <<<
        set -ex -o pipefail

        samtools view -bo file.bam ~{sam}
        samtools sort -o sorted.bam file.bam
        samtools index sorted.bam

        # count reads aligned to ref
        samtools view -c sorted.bam > total_count.txt
        samtools view -F 4 sorted.bam | cut -f 3 | sort | uniq -c | sort -rnk1 > matched_count.txt

        # drop "/n" in in the reference fasta
        awk '/^>/ {if (NR!=1) {printf("\n%s\n", $0)} \
        else {printf("%s\n",$0)}; next} { printf("%s",$0)} END {printf("\n")}' \
        ~{ref} > ref.fasta

        # filter only reads aligned to reference
        header=$(grep -o -P '^>\K[^ ]+' ~{ref})
        samtools view -h -F 4 sorted.bam $header > filtered.sam

        samtools view -bo filtered.bam filtered.sam
        samtools sort filtered.bam -o sorted.bam
        samtools index sorted.bam

        samtools mpileup -uf ~{ref} sorted.bam | bcftools call -mv -Oz -o file.bcf
        bcftools index file.bcf
        bcftools consensus -f ~{ref} file.bcf > consensus.fasta
        bcftools view -Ov -o file.vcf file.bcf
        bcftools view file.bcf | bcftools filter -i 'DP>~{min_cov}' > filtered.vcf

        rm file.bam file.bcf file.bcf.csi filtered.bam filtered.sam ref.fasta sorted.bam sorted.bam.bai

   >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File? file_vsf = "file.vcf"
        File? filtered_vsf = "filtered.vcf"
        File? consensus_fasta = "HA_consensus.fasta"
        File? matched_count_txt = "matched_count.txt"
        File? total_count_txt = "total_count.txt"
    }

}

task Snpeff {

    input {
        File snpeff_config
        File snpeff_db
        File? vcf_file
        File ref
        String docker
    }

    command <<<
        set -ex -o pipefail

        unzip ~{snpeff_db}

        ref_name=$(basename "~{ref}" .fasta)

        echo $ref_name > ref_name.txt

        header=$(grep -o -P '^>\K[^ ]+' ~{ref})

        java -Xmx8g -jar /opt/snpEff/snpEff.jar \
        -c ~{snpeff_config} \
        -dataDir "${PWD}/data" \
        "$header" -v \
        "~{vcf_file}" \
        -onlyProtein -noMotif -noNextProt \
        -noStats -classic -no-downstream \
        -no-intergenic -no-intron -no-upstream -no-utr \
        > ann.vcf

        java -jar /opt/snpEff/SnpSift.jar extractFields -v -e '#' ann.vcf EFF \
        | tail -n +2 \
        | sed '{/^#/d;s/|/\t/g;s/(/\t/;s/)//}' \
        | awk -F "\t" '{print $7,":",$5}' \
        | sed 's/ : /:/' \
        | tr '\n' ',' \
        | sed 's/,$/\n/' > ann.csv

        rm -r data

    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
       File snpeff_vcf = "ann.vcf"
       File snpeff_csv = "ann.csv"
       String ref_name = read_string("ref_name.txt")
    }

}

task ParseSnpeff {

    input {
        File snpeff_csv
        String ref_name
        File? antigenic_frame
        String docker
    }

    command <<<
        set -ex -o pipefail

        python3 <<CODE

        import pandas as pd
        import json
        from typing import Dict, Optional, Tuple, List
        import os.path
        from pathlib import Path
        import re
        import csv

        default_frame=[128, 129, 156, 157, 158, 159, 160, 162, 163, 164, 165, \
               166, 167, 187, 188, 189, 190, 191, 192, 193, 194, 195, \
               196, 197, 198, 169, 170, 171, 172, 173, 206, 207, 208, \
               238, 239, 240, 141, 142, 143, 144, 145, 224, 225, 74, 75, 76, 77, 78, 79]

        HA_pattern = re.compile(r"HA")

        if not os.path.isfile("~{antigenic_frame}"):
            antigenic_frame=default_frame
        else:
            read_frame = [i for i in Path("~{antigenic_frame}").read_text(encoding="utf-8").replace("\n", " ").split()]
            antigenic_frame = [int(item) for item in read_frame]

        def parse_mutations (
            snpeff_csv: Optional[str],
            antigenic_frame: Optional[List[int]]
        ):
            antigenic_frame.sort()
            with open(snpeff_csv, 'r') as file:
                reader = csv.reader(file)
                mut = []
                for row in reader:
                    mut.extend(row)
            if len(mut) == 0:
                mut = "-"
                ant_mut = "-"
            else:
                ant_mut = []
                for m in mut:
                    dig = int(re.search( r'(?<=:(\w|\*))[0-9]*', m).group())
                    if dig in antigenic_frame:
                        ant_mut.append(m)
                mut = " ".join(mut)
                if len(ant_mut) == 0:
                    ant_mut = "-"
                else:
                    ant_mut = " ".join(ant_mut)
            return mut, ant_mut

        def write_report (
            ref_name: Optional[str],
            snpeff_csv: Optional[str],
            antigenic_frame: Optional[List[str]]

        ):
            mut, ant_mut = parse_mutations (snpeff_csv, antigenic_frame)
            if HA_pattern.search(ref_name):
                with open("snpeff.json", 'w') as json_file:
                    data = {
                        f"mutations_to_{ref_name}_snpeff" : mut,
                        f"antigenic_site_mutations_to_{ref_name}_snpeff" : ant_mut
                    }
                    json.dump(data, json_file, indent=4)
            else:
                with open("snpeff.json", 'w') as json_file:
                    data = {
                        f"mutations_to_{ref_name}_snpeff" : mut
                    }
                    json.dump(data, json_file, indent=4)


        write_report ("~{ref_name}", "~{snpeff_csv}", antigenic_frame)

        CODE
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
       File snpeff_json = "snpeff.json"
    }
}







