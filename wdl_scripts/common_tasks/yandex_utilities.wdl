version 1.0

task trimm {
    input {
        String sample_id
        File fastq_1
        File fastq_2
        Int lines_number
        Int compression_level
        Int max_retries
        Int? minlen = 36
        String docker
    }

    command <<<
        set -e
        /opt/fastp/fastp \
            -i ~{fastq_1} \
            -I ~{fastq_2} \
            --out1 outR1.fastq.gz \
            --out2 outR2.fastq.gz \
            --qualified_quality_phred 15 \
            --unqualified_percent_limit 40 \
            --length_required ~{minlen} \
            --detect_adapter_for_pe \
            --cut_tail \
            --cut_tail_window_size 6 \
            --cut_mean_quality 15 \
            --dont_eval_duplication \
            --compression 2 \
            --trim_poly_g \
            --thread 16
    >>>
	runtime {
		docker: "~{docker}"
		maxRetries: max_retries
	}
	output {
		File R1_file = "outR1.fastq.gz"
		File R2_file = "outR2.fastq.gz"
		File fastp_json = "fastp.json"
		File fastp_html = "fastp.html"
    }
}

task Minimap2 {
    input {
        File? fastq_1
        File? fastq_2
        File ref
        Int threads
        String docker
    }

    command <<<
        set -e
        minimap2 -t ~{threads} \
        -ax sr -c \
        --secondary=yes \
        ~{ref} \
        ~{fastq_1} ~{fastq_2} > file.sam

        # STATS counting
        grep -v '^@' file.sam | wc -l > total_count.txt
        grep -v '^@' file.sam | cut -f 3 | sort | uniq -c | awk '$2 != "*" {print $1, $2}' > matched_count.txt
        #sed "s/:/#/g" matched_count.txt > modified_matched_count.txt # for transfer only, ":" in fasta's headers interpreted as a path 
    >>>
    runtime {
        docker: "~{docker}"
    }
    output {
        File file_sam = "file.sam"
        File total_count_txt = "total_count.txt"
        File matched_count_txt = "matched_count.txt"
    }
}

task Minimap2Parse {
    input {
        String docker
        File matched_count
    }

    command <<<
        set -e
        python3 <<CODE

        import pandas as pd
        import json


        CODE
    >>>
    runtime {
        docker: "~{docker}"
    }
    output {
        File? matched_count_json = "output.json"
    }
}

