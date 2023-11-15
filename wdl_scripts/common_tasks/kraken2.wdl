version 1.0

task Kraken2 {
    input {
        File fastq_1
        File fastq_2
        String sample_name
        File kraken2_classifier
        Int threads
        String docker
    }

    command <<<
        set -ex -o pipefail

        DB_NAME=$(basename "~{kraken2_classifier}" .tar.gz)

        tar xf ~{kraken2_classifier}

        kraken2 --db $DB_NAME \
        --threads ~{threads} \
        --report ~{sample_name}_report.txt \
        --output ~{sample_name}_output.tsv \
        --paired ~{fastq_1} ~{fastq_2}

        rm -r $DB_NAME
    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File report_txt = "~{sample_name}_report.txt"
        File output_tsv = "~{sample_name}_output.tsv"
    }
}

task Bracken {
    input {
        String sample_name
        File kraken_report
        File kraken2_classifier
        String level
        Int threshold = 0 # 10 by default - in cases where there's a largely diverse sample, 10 the number of reads at the strain/genome level required, not at the genus/species levels. 
        String docker
    }

    command <<<
        set -ex -o pipefail

        DB_NAME=$(basename "~{kraken2_classifier}" .tar.gz)
        tar xf ~{kraken2_classifier}

        cd $DB_NAME
        BRACKEN_LENGTH=$(ls -v1 *.kmer_distrib | tail -1 | sed -e "s@^database@@" -e "s@mers.kmer_distrib@@")
        cd ..

        bracken -d $DB_NAME \
        -i ~{kraken_report} \
        -o report.txt \
        -r $BRACKEN_LENGTH \
        -l ~{level} \
        -t ~{threshold}

        head -n 1 report.txt > ~{sample_name}_sorted_bracken_report.txt
        tail -n +2 report.txt | sort -k7 -n -r -t$'\t' >> ~{sample_name}_sorted_bracken_report.txt

        rm -r $DB_NAME
    >>>

    runtime {
        docker:"~{docker}"
    }

     output {
        File report_txt = "~{sample_name}_sorted_bracken_report.txt"
    }
}

task Krona {
    input {
        String sample_name
        File report
        String docker
    }

    command <<<
        set -ex -o pipefail

        ktImportTaxonomy -o ~{sample_name}_krona.html ~{report}

    >>>

    runtime {
        docker: "~{docker}"
    }

    output {
        File report_html = "~{sample_name}_krona.html"
    }
}
