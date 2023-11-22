version 1.0

task CheckInput {

    input {
        File fastq
        Int min_reads
        String docker
    }

    command <<<
        set -euxo pipefail

        count_in="$(zcat ~{fastq} | wc -l)"
        count_in=$((count_in / 4))

        if (($count_in>~{min_reads})); then
          proceed="yes"
        else
          proceed="no"
        fi

        echo $proceed > proceed.txt
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        String proceed = read_string("proceed.txt")
    }
}

task FastQC {

    input {
        File? fastq
        String docker
    }

    command <<<
        set -ex -o pipefail

        chmod 755 /opt/fastqc/fastqc
        JAVA=$(which java)
        /opt/fastqc/fastqc ~{fastq} -java=$JAVA --noextract --outdir "."
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File summary_html = select_first(glob("*.html"))
        File zip = select_first(glob("*.zip"))
    }
}

task Trimmomatic {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        Int threads
        String docker
    }

    command <<<
        set -ex -o pipefail

        java -jar /opt/trimmomatic/trimmomatic-0.39.jar PE \
        -threads ~{threads} \
        ~{fastq_1} ~{fastq_2} \
        R1_trimmed.fastq.gz R1_unpaired.fastq.gz \
        R2_trimmed.fastq.gz R2_unpaired.fastq.gz \
        ILLUMINACLIP:/opt/trimmomatic/adapters/TruSeq3-PE-2.fa:4:30:5:2:keepBothReads \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:6:15 MINLEN:36 >> ~{sample_name}_trimmomatic_stats.txt 2>&1

        grep "Input Read Pairs:" ~{sample_name}_trimmomatic_stats.txt | \
        sed 's/\(Input Read Pairs: [0-9]*\).*/\1/' > input_read_pairs.txt

        grep "Input Read Pairs:" ~{sample_name}_trimmomatic_stats.txt | \
        sed -n 's/.*\(Both Surviving: [0-9]* (.*%)\).*/\1/p' | awk -F ' Forward' '{print $1}' > both_surviving.txt
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File trim_fastq_1 = "R1_trimmed.fastq.gz"
        File trim_fastq_2 = "R2_trimmed.fastq.gz"
        File summary_txt = "~{sample_name}_trimmomatic_stats.txt"
        String input_reads_trim = read_string("input_read_pairs.txt")
        String both_surviving_trim = read_string("both_surviving.txt")
    }
}

task Cutadapt {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        File? primer_left
        File? primer_right
        String docker
    }

    command <<<
        set -ex -o pipefail

        . /opt/cutadapt/bin/activate

        /opt/cutadapt/bin/cutadapt \
        --cores=0 \
        -a file:~{primer_left} \
        -A file:~{primer_right} \
        -o cut_R1.fastq.gz -p cut_R2.fastq.gz \
        --discard-untrimmed \
        ~{fastq_1} \
        ~{fastq_2} >> ~{sample_name}_cutadapt_stats.txt 2>&1

        count_in=$(grep "Total read pairs processed" ~{sample_name}_cutadapt_stats.txt \
        | grep -o '[0-9,]*' | tr -d ',')
        echo "Input Read Pairs: $count_in" > input_read_pairs.txt

        count_out=$(grep "passing filters" ~{sample_name}_cutadapt_stats.txt \
        | cut -d ":" -f2 | sed 's/^[[:space:]]*//' | tr -d ',')
        echo "Both Surviving: $count_out" > both_surviving.txt
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File cut_fastq_1 = "cut_R1.fastq.gz"
        File cut_fastq_2 = "cut_R2.fastq.gz"
        File summary_txt = "~{sample_name}_cutadapt_stats.txt"
        String input_reads_cut = read_string("input_read_pairs.txt")
        String both_surviving_cut = read_string("both_surviving.txt")
    }
}

task HostFilter {

    input {
        File? fastq_1
        File? fastq_2
        String sample_name
        File index_tar
        Int threads
        String docker
    }

    command <<<
        set -euxo pipefail

        count_in="$(zcat ~{fastq_1} ~{fastq_2} | wc -l)"
        count_in=$((count_in / 8))
        jq --null-input --arg count_in "$count_in" '{"reads_with_host":$count_in}' > ~{sample_name}_hostfilter_stats.txt
        echo "Input Read Pairs: $count_in" > input_read_pairs.txt

        db_name=$(basename "~{index_tar}" .bowtie2.tar)

        tar xf '~{index_tar}' -C /tmp

        /opt/bowtie2/bowtie2 \
        -p ~{threads} \
        -x /tmp/$db_name/$db_name \
        -1 ~{fastq_1} \
        -2 ~{fastq_2} \
        --un-conc-gz \
        sample_host_removed \
        > sample_mapped_and_unmapped.sam

        mv sample_host_removed.1 sample_host_removed_R1.fastq.gz
        mv sample_host_removed.2 sample_host_removed_R2.fastq.gz

        count_out="$(zcat sample_host_removed_R{1,2}.fastq.gz | wc -l)"
        count_out=$((count_out / 8))
        jq --null-input --arg count_out "$count_out" '{"host_filtered_reads":$count_out}' >> ~{sample_name}_hostfilter_stats.txt

        percentage=$((100 * count_out / count_in))
        echo "Both Surviving: $count_out ($percentage%)" > both_surviving.txt
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File summary_txt = "~{sample_name}_hostfilter_stats.txt"
        File host_filtered_fastq_1 = "sample_host_removed_R1.fastq.gz"
        File host_filtered_fastq_2 = "sample_host_removed_R2.fastq.gz"
        String input_reads_host = read_string("input_read_pairs.txt")
        String both_surviving_host = read_string("both_surviving.txt")
    }
}

task PreprocessingQC {

    input {
        String sample_name
        String trim_input_reads
        String trim_both_surviving
        String? cut_input_reads
        String? cut_both_surviving
        String host_input_reads
        String host_both_surviving
        String docker
    }

    command <<<
        set -ex -o pipefail

        python3 <<CODE
        import json

        data = {
            "trimmomatic": "~{trim_input_reads}, ~{trim_both_surviving}",
            "cutadapt": "~{cut_input_reads}, ~{cut_both_surviving}",
            "host_filtering": "~{host_input_reads}, ~{host_both_surviving}"
        }

        with open("~{sample_name}_preprocessing_qc.json", "w") as json_file:
            json.dump(data, json_file, indent=4)
        CODE
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File report_json = "~{sample_name}_preprocessing_qc.json"
    }
}
