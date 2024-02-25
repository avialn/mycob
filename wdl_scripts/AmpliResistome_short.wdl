version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run wdl_scripts/AmpliResistome_short.wdl -i inputs.json -o options.json
#https://github.com/chanzuckerberg/czid-workflows/blob/main/workflows/amr/run.wdl
#https://chanzuckerberg.zendesk.com/hc/en-us/articles/15091031482644-AMR-Pipeline-Workflow

import "./common_tasks/preprocessing.wdl" as preprocessing
import "./common_tasks/kraken2.wdl" as kraken2

import "./common_tasks/amr_tasks.wdl" as tasks
import "./common_tasks/yandex_utilities.wdl" as Utils

workflow AmpliResistome {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        Int threads = 8
        Boolean cut_primers = false
        File? primer_left = "primers/left_primers.fasta"
        File? primer_right = "primers/right_primers.fasta"
        Boolean transcriptome_filtering = true # false -> index_genome, true -> index_transcriptome
        File index_genome = "db/bowtie2/GRCh38_ERCC.bowtie2.tar"
        File index_transcriptome = "db/bowtie2/GRCh38_transcriptome.bowtie2.tar"
        File card_json = "/home/admin/mycob/db/amr/CARD/card.json"
        File kraken2_standard_8gb = "db/kraken2/k2_standard_08gb_20230605.tar.gz"
        Int min_contig_length = 500
        String kraken_level = "S" #(U)nclassified, (R)oot, (D)omain, (K)ingdom (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
        File minimap_ref = "/home/admin/antares2.fasta"
    }

    call preprocessing.FastQC as fastqc_row_R1 {
        input:
            fastq = fastq_1,
            docker = "fastqc:4.4"
    }

    call preprocessing.FastQC as fastqc_row_R2 {
        input:
            fastq = fastq_2,
            docker = "fastqc:4.4"
    }

    call preprocessing.Trimmomatic as trimmomatic {
        input:
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            sample_name = sample_name,
            threads = threads,
            docker = "trimmomatic:4.5"
    }

    if (cut_primers) {
        call preprocessing.Cutadapt as cutadapt {
            input:
                fastq_1 = trimmomatic.trim_fastq_1,
                fastq_2 = trimmomatic.trim_fastq_2,
                sample_name = sample_name,
                primer_left = primer_left,
                primer_right = primer_right,
                docker = "cutadapt:4.4"
        }
    }

    call preprocessing.HostFilter as host_filter {
        input:
            fastq_1 = if cut_primers then cutadapt.cut_fastq_1 else trimmomatic.trim_fastq_1,
            fastq_2 = if cut_primers then cutadapt.cut_fastq_2 else trimmomatic.trim_fastq_2,
            sample_name = sample_name,
            index_tar = if transcriptome_filtering then index_transcriptome else index_genome,
            threads = threads,
            docker = "bowtie2:4.4"
    }

    call preprocessing.PreprocessingQC as preprocessing_qc {
        input:
            sample_name = sample_name,
            trim_input_reads = trimmomatic.input_reads_trim,
            trim_both_surviving = trimmomatic.both_surviving_trim,
            cut_input_reads = if defined(cutadapt.input_reads_cut) then cutadapt.input_reads_cut else "cut_primers is false",
            cut_both_surviving = if defined(cutadapt.both_surviving_cut) then cutadapt.both_surviving_cut else "no primers were cut",
            host_input_reads = host_filter.input_reads_host,
            host_both_surviving = host_filter.both_surviving_host,
            docker = "python:4.4"

    }

    call preprocessing.FastQC as fastqc_trimed_R1 {
        input:
            fastq = host_filter.host_filtered_fastq_1,
            docker = "fastqc:4.4"
    }

    call preprocessing.FastQC as fastqc_trimed_R2 {
        input:
            fastq = host_filter.host_filtered_fastq_2,
            docker = "fastqc:4.4"
    }

    #call kraken2.Kraken2 as kraken2 {
    #    input:
    #        fastq_1 = host_filter.host_filtered_fastq_1,
    #        fastq_2 = host_filter.host_filtered_fastq_2,
    #        sample_name = sample_name,
    #        kraken2_classifier = kraken2_standard_8gb,
    #        threads = threads,
    #        docker = "staphb/kraken2:latest"
    #}

    #call kraken2.Bracken as bracken {
    #    input:
    #        sample_name = sample_name,
    #        kraken_report = kraken2.report_txt,
    #        kraken2_classifier = kraken2_standard_8gb,
    #        level = kraken_level,
    #        docker = "nanozoo/bracken:2.8--dcb3e47"
    #}

    #call kraken2.Krona as krona_kraken {
    #    input:
    #    sample_name = sample_name,
    #    report = kraken2.report_txt,
    #    docker = "krona:2.8.1"
    #}

    call tasks.RgiBwt as rig_bwt {
        input:
            fastq_1 = host_filter.host_filtered_fastq_1,
            fastq_2 = host_filter.host_filtered_fastq_2,
            card_json = card_json,
            docker = "rgi:1.1"
    }

    call Utils.Minimap2 as minimap {
        input:
            fastq_1 = host_filter.host_filtered_fastq_1,
            fastq_2 = host_filter.host_filtered_fastq_2,
            ref = minimap_ref,
            threads = threads,
            docker = "staphb/minimap2:latest"
    }

    call Utils.Minimap2Parse as minimap_parse {
        input:
            matched_reads_count = minimap.matched_count_txt,
            total_reads_count = minimap.total_count,
            docker = "python:4.4"
    }

    output {
        # preprocessing: fastqc, trimmomatic, cutadapt, host-filtering
        File fastqc_row_R1_html = fastqc_row_R1.summary_html
        File fastqc_row_R2_html = fastqc_row_R2.summary_html
        File fastqc_trimed_R1_html = fastqc_trimed_R1.summary_html
        File fastqc_trimed_R2_html = fastqc_trimed_R2.summary_html
        File preprocessing_qc_json = preprocessing_qc.report_json

        #rgi_bwt
        File allele_mapping_data_txt = rig_bwt.kma_amr_results_txt
        File artifacts_mapping_stats_txt = rig_bwt.artifacts_mapping_stats
        File gene_mapping_data_txt = rig_bwt.gene_mapping_data
        File overall_mapping_stats_txt = rig_bwt.overall_mapping_stats
        File reference_mapping_stats_txt = rig_bwt.reference_mapping_stats
        File allele_mapping_data_json = rig_bwt.kma_amr_results_json

        # kraken
        #File kraken_txt = kraken2.report_txt
        #File bracken_txt = bracken.report_txt

        #minimap
        String minimap_total_reads = minimap.total_count # host_filtered x 2
        File minimap_matched_reads_json = minimap_parse.matched_count_json
    }
}
