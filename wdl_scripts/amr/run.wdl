version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run run.wdl -i inputs.json -o options.json
#https://github.com/chanzuckerberg/czid-workflows/blob/main/workflows/amr/run.wdl
#https://chanzuckerberg.zendesk.com/hc/en-us/articles/15091031482644-AMR-Pipeline-Workflow

import "amr_tasks.wdl" as tasks
import "yandex_utilities.wdl" as Utils

workflow processing {

    input {
        Array[Array[String]] Files
	    Int max_retries = 1
	    Int compression_level = 5
        String TaskID
        String Attempt
        String TableID = TaskID + "-" + Attempt
        String SampleID = TableID
        String SampleName = TableID
        String sample_name = SampleName
        Int threads = 8
        Boolean cut_primers = false
        File? adapter_f = "/home/cromwell/mycob-ref/rsv_full/left_primers.fasta"
        File? adapter_r = "/home/cromwell/mycob-ref/rsv_full/right_primers.fasta"
        File index_tar = "/home/cromwell/mycob-ref/rsv_full/human_bowtie2_index/GRCh38_ERCC.bowtie2.tar"
        File card_json = "/home/cromwell/mycob-ref/amr/CARD/card.json"
        Int min_contig_length = 500
        Int lines_number = length(Files)
    }


    call Utils.trimm as trimm {
        input:
            fastq_1 = Files[0][0],
            fastq_2 = Files[0][1],
            lines_number = lines_number,
            sample_id = SampleID,
            minlen = 36,
            compression_level = compression_level,
            max_retries = max_retries,
            docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastp:0.23.4"
    }

    String fastq_1 = trimm.R1_file
    String fastq_2 = trimm.R2_file


    call tasks.FastQC as row_R1_FastQC {
        input:
        fastq = Files[0][0],
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call tasks.FastQC as row_R2_FastQC {
        input:
        fastq = Files[0][1],
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call tasks.Trimmomatic as trimmomatic {
        input:
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        threads = threads,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/trimmomatic:0.39"
    }

    if (cut_primers) {
        call tasks.Cutadapt as cutadapt {
            input:
            fastq_1 = trimmomatic.trim_fastq_1,
            fastq_2 = trimmomatic.trim_fastq_2,
            adapter_f = adapter_f,
            adapter_r = adapter_r,
            docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/cutadapt:4.4"
        }
    }

    call tasks.FastQC as trimmed_R1_FastQC {
        input:
        fastq = if cut_primers then cutadapt.cut_fastq_1 else trimmomatic.trim_fastq_1,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call tasks.FastQC as trimmed_R2_FastQC {
        input:
        fastq = if cut_primers then cutadapt.cut_fastq_2 else trimmomatic.trim_fastq_2,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call tasks.Bowtie2Filter as bowtie2_filter {
        input:
        fastq_1 = if cut_primers then cutadapt.cut_fastq_1 else trimmomatic.trim_fastq_1,
        fastq_2 = if cut_primers then cutadapt.cut_fastq_2 else trimmomatic.trim_fastq_2,
        index_tar = index_tar,
        threads = threads,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/bowtie2:2.5.1"
    }

    call tasks.SamtoolsFilter as samtools_filter {
        input:
        bowttie2_output_sam = bowtie2_filter.output_sam,
        threads = threads,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/samtools:1.9"
    }

    call tasks.RgiBwt as rig_bwt {
        input:
        fastq_1 = samtools_filter.host_filtered_fastq_1,
        fastq_2 = samtools_filter.host_filtered_fastq_2,
        card_json = card_json,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/rgi:3.2.6"
    }

    call tasks.Flash as flash {
        input:
        fastq_1 = samtools_filter.host_filtered_fastq_1,
        fastq_2 = samtools_filter.host_filtered_fastq_2,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/flash:1.2.11"
    }

    call tasks.SeqtkToFa as seqtk_to_fa {
        input:
        flash_output_fq = flash.output_fq,
        docker ="cr.yandex/crpl2lv1lkr7g21e6q8g/seqtk:1.4-r122"
    }

    call tasks.Spades as spades {
        input:
        seqtk_output_fa  = seqtk_to_fa.output_fa,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/spades:3.15.5"
    }

    call tasks.SeqtkFilter as seqtk_filter {
        input:
        contigs = spades.contigs_fa,
        min_contig_length = min_contig_length,
        docker ="cr.yandex/crpl2lv1lkr7g21e6q8g/seqtk:1.4-r122"
    }

    call tasks.RgiMain as rgi_main {
        input:
        contigs_filt = seqtk_filter.contigs_filtered,
        card_json = card_json,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/rgi:3.2.6"
    }

    call tasks.GeneCoverage as gene_coverage {
        input:
        main_amr_results = rgi_main.output_txt,
        main_output_json = rgi_main.output_json,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/python:3"
    }

    call tasks.Report as report {
        input:
        main_output = rgi_main.output_txt,
        kma_output = rig_bwt.kma_amr_results_txt,
        gene_coverage_tsv = gene_coverage.output_tsv,
        sample_name = sample_name,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/python:3"
    }

    output {
        File r1_row_fastqc_html = row_R1_FastQC.qc_report_html
        File r2_row_fastqc_html = row_R2_FastQC.qc_report_html
        File? cutadapt_summary_tsv = cutadapt.cutadapt_summary_tsv
        File r1_trimmed_fastqc_html = trimmed_R1_FastQC.qc_report_html
        File r2_trimmed_fastqc_html = trimmed_R2_FastQC.qc_report_html
        File number_row_reads = bowtie2_filter.row_reads_count
        File number_host_filtered_reads = samtools_filter.host_filtered_reads_count
        File allele_mapping_data_txt = rig_bwt.kma_amr_results_txt
        File artifacts_mapping_stats_txt = rig_bwt.artifacts_mapping_stats
        File gene_mapping_data_txt = rig_bwt.gene_mapping_data
        File overall_mapping_stats_txt = rig_bwt.overall_mapping_stats
        File reference_mapping_stats_txt = rig_bwt.reference_mapping_stats
        File allele_mapping_data_json = rig_bwt.kma_amr_results_json
        File spades_contigs_fa = spades.contigs_fa
        File? spades_scaffolds_fa = spades.scaffolds_fa
        Array[File] contig_amr_report = rgi_main.output_main
        File main_report_tsv = report.contigs_tsv
        File kma_report_tsv= report.reads_tsv
        File comprehensive_AMR_metrics_tsv = report.final_summary_tsv
        File bigtable_report_tsv = report.bigtable_tsv
        File primary_AMR_report_tsv = report.synthesized_report_tsv
    } 


}
