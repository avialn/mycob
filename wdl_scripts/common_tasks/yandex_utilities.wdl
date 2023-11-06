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
            -i ${fastq_1} \
            -I ${fastq_2} \
            --out1 outR1.fastq.gz \
            --out2 outR2.fastq.gz \
            --qualified_quality_phred 15 \
            --unqualified_percent_limit 40 \
            --length_required ${minlen} \
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
		File R2_file = "*outR2.fastq.gz"
		File fastp_json = "fastp.json"
		File fastp_html = "fastp.html"
    }
}