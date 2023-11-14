version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run wdl_scripts/FluGenotyping.wdl -i inputs.json -o options.json

import "./common_tasks/preprocessing.wdl" as preprocessing
import "./common_tasks/kraken2.wdl" as kraken2
import "./common_tasks/irma.wdl" as irma
import "./common_tasks/nextclade.wdl" as nextclade

workflow FluGenotyping {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        Boolean cut_primers = false
        File? primer_left = "primers/left_primers.fasta"
        File? primer_right = "primers/right_primers.fasta"
        Boolean transcriptome_filtering = true # false -> index_genome, true -> index_transcriptome
        File index_genome = "db/bowtie2/GRCh38_ERCC.bowtie2.tar"
        File index_transcriptome = "db/bowtie2/GRCh38_transcriptome.bowtie2.tar"
        File kraken2_standard_8gb = "db/kraken2/k2_standard_08gb_20230605.tar.gz"
        File kraken2_virus = "db/kraken2/k2_viral_20231009.tar.gz"
        String kraken_level = "S" #(U)nclassified, (R)oot, (D)omain, (K)ingdom (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
        Int threads = 8
        File reference_fasta_flu = "reference_fasta/FLU.fasta"
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

    call kraken2.Kraken2 as kraken2 {
        input:
        fastq_1 = host_filter.host_filtered_fastq_1,
        fastq_2 = host_filter.host_filtered_fastq_2,
        sample_name = sample_name,
        kraken2_classifier = kraken2_standard_8gb,
        threads = threads,
        docker = "staphb/kraken2:latest"
    }

    call kraken2.Bracken as bracken {
        input:
        sample_name = sample_name,
        kraken_report = kraken2.report_txt,
        kraken2_classifier = kraken2_standard_8gb,
        level = kraken_level,
        docker = "nanozoo/bracken:2.8--dcb3e47"
    }

    call kraken2.Krona as krona_kraken {
        input:
        sample_name = sample_name,
        report = kraken2.report_txt,
        docker = "krona:2.8.1"
    }

    call kraken2.Kraken2 as kraken2_vir {
        input:
        fastq_1 = host_filter.host_filtered_fastq_1,
        fastq_2 = host_filter.host_filtered_fastq_2,
        sample_name = sample_name,
        kraken2_classifier = kraken2_virus,
        threads = threads,
        docker = "staphb/kraken2:latest"
    }

    call kraken2.Bracken as bracken_vir {
        input:
        sample_name = sample_name,
        kraken_report = kraken2_vir.report_txt,
        kraken2_classifier = kraken2_virus,
        level = kraken_level,
        threshold = 0,
        docker = "nanozoo/bracken:2.8--dcb3e47"
    }

    call irma.Irma as irma_flu {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "FLU",
        docker = "irma:230723"
    }

    if (irma_flu.proceed == "yes") {

        call irma.IrmaQC as irma_qc_flu {
            input:
            sample_name = sample_name,
            module_irma = "FLU",
            irma_qc = irma_flu.irma_qc,
            reference_fasta = reference_fasta_flu,
            irma_fasta = irma_flu.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_flu {
            input:
            sample_name = sample_name,
            module_irma = "FLU",
            fasta_files = irma_flu.irma_fasta,
            db = "influenzaDB",
            threads = threads,
            docker = "blast:4.5"
        }

        call irma.Report as report_flu {
            input:
            sample_name = sample_name,
            module_irma = "FLU",
            blast_res = blast_flu.blast_res,
            irma_type = irma_flu.flu_type,
            docker = "python:4.4"
        }

        call nextclade.Nextclade as nextclade_flu {
            input:
            ha_fasta = irma_flu.ha_fasta,
            na_fasta = irma_flu.na_fasta,
            docker = "nextclade:3.3"
        }

        call nextclade.NextcladeParse as nextclade_parse_flu {
            input:
            HA_tsv = nextclade_flu.HA_nextclade,
            NA_tsv = nextclade_flu.NA_nextclade,
            report = report_flu.report,
            antigenic_frame = "/home/admin/cromwell_flu/inputs/VD_6/antigenic_frame.txt",
            docker = "python:4.4"
        }
    }

    output {
        File fastqc_row_R1_html = fastqc_row_R1.summary_html
        File fastqc_row_R2_html = fastqc_row_R2.summary_html
        File fastqc_trimed_R1_html = fastqc_trimed_R1.summary_html
        File fastqc_trimed_R2_html = fastqc_trimed_R2.summary_html
        File preprocessing_qc_json = preprocessing_qc.report_json
        File? flu_qc_json = irma_qc_flu.irma_qc_json
        File? flu_json = report_flu.report
        File kraken_txt = kraken2.report_txt
        File bracken_txt = bracken.report_txt
        File krona_kraken_html = krona_kraken.report_html
        File kraken_virus_txt = kraken2_vir.report_txt
        File bracken_virus_txt = bracken_vir.report_txt
    }
}
