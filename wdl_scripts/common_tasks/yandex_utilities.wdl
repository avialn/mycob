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
        String total_count = read_string("total_count.txt")
        File matched_count_txt = "matched_count.txt"
    }
}

task Minimap2Parse {
    input {
        String docker
        File matched_reads_count
        String total_reads_count
    }

    command <<<
        set -e
        python3 <<CODE

        import pandas as pd
        import json
        
        def create_matched_count_json(df):
            """creates dict.json from dataframe of matched reads"""
            dict_json = {}
            if df.empty:
                dict_json["matched_reads_count"] = 0
            else:
                dict_json = dict(zip (df["fasta_name"], df["reads_count"]))
            return dict_json

        def create_merged_json(count, dict_json):
            """merges total_count.json and matched_count.json to merged.json"""
            return {
                "total_reads": int(count),
                "matched_reads" : dict_json
            }

        def write_json(dict_json, file_json):
            """writes file.json from data.json"""
            with open(file_json, "w") as file:
                json.dump(dict_json, file, indent =4)
          
        matched_count_df = pd.read_csv(
            "~{matched_reads_count}", header=None, sep = " ", names=["reads_count", "fasta_name"]
        )
        sorted_matched_count_df = matched_count_df.sort_values(by="reads_count", ascending=False)
        matched_count_json = create_matched_count_json(sorted_matched_count_df)
        merged_json = create_merged_json(~{total_reads_count}, matched_count_json)
        write_json(merged_json, "matched_reads_count.json")
        CODE
    >>>
    runtime {
        docker: "~{docker}"
    }
    output {
        File matched_count_json = "matched_reads_count.json"
    }
}

