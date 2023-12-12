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
            echo "0" > matched_count.txt
        else
            echo "yes" > proceed.txt
            echo $matched_count > matched_count.txt
        fi

    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File file_sam = "file.sam"
        String proceed = read_string("proceed.txt")
        File matched_count_txt = "matched_count.txt"
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
        samtools view -F 4 sorted.bam | cut -f 3 | sort | uniq -c | sort -rnk1 > count.txt

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
        File? count_txt = "count.txt"
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
    }

}


