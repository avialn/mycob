version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run wdl_scripts/VirGenotyping.wdl -i inputs.json -o options.json

import "./common_tasks/preprocessing.wdl" as preprocessing
import "./common_tasks/kraken2.wdl" as kraken2
import "./common_tasks/irma.wdl" as irma
import "./common_tasks/nextclade.wdl" as nextclade

workflow VirGenotyping {

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
        File reference_fasta_measles = "reference_fasta/MEASLES.fasta"
        File reference_fasta_adeno = "reference_fasta/ADENO.fasta"
        File reference_fasta_pneumo = "reference_fasta/PNEUMO.fasta"
        File reference_fasta_boca = "reference_fasta/BOCA.fasta"
        File reference_fasta_corona = "reference_fasta/CORONA.fasta"
        File reference_fasta_respiro = "reference_fasta/RESPIRO.fasta"
        File reference_fasta_rhino = "reference_fasta/RHINO.fasta"
        File reference_fasta_rubula = "reference_fasta/RUBULA.fasta"
        File reference_fasta_cov = "reference_fasta/CoV.fasta"
        File reference_fasta_metapneumo = "reference_fasta/METAPMEUMO.fasta"
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


    call irma.Irma as irma_measles {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "MEASLES",
        docker = "irma:230723"
    }

    if (irma_measles.proceed == "yes") {

        call irma.IrmaQC as irma_qc_measles {
            input:
            sample_name = sample_name,
            module_irma = "MEASLES",
            irma_qc = irma_measles.irma_qc,
            reference_fasta = reference_fasta_measles,
            irma_fasta = irma_measles.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_measles {
            input:
            sample_name = sample_name,
            module_irma = "MEASLES",
            fasta_files = irma_measles.irma_fasta,
            db = "MeaslesDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_measles {
            input:
            sample_name = sample_name,
            module_irma = "MEASLES",
            blast_res = blast_measles.blast_res,
            irma_type = irma_qc_measles.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_adeno {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ADENO",
        docker = "irma:230723"
    }

    if (irma_adeno.proceed == "yes") {

        call irma.IrmaQC as irma_qc_adeno {
            input:
            sample_name = sample_name,
            module_irma = "ADENO",
            irma_qc = irma_adeno.irma_qc,
            reference_fasta = reference_fasta_adeno,
            irma_fasta = irma_adeno.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_adeno {
            input:
            sample_name = sample_name,
            module_irma = "ADENO",
            fasta_files = irma_adeno.irma_fasta,
            db = "MastadenovirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_adeno {
            input:
            sample_name = sample_name,
            module_irma = "ADENO",
            blast_res = blast_adeno.blast_res,
            irma_type = irma_qc_adeno.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_pneumo {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "PNEUMO",
        docker = "irma:230723"
    }

    if (irma_pneumo.proceed == "yes") {

        call irma.IrmaQC as irma_qc_pneumo {
            input:
            sample_name = sample_name,
            module_irma = "PNEUMO",
            irma_qc = irma_pneumo.irma_qc,
            reference_fasta = reference_fasta_pneumo,
            irma_fasta = irma_pneumo.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_pneumo {
            input:
            sample_name = sample_name,
            module_irma = "PNEUMO",
            fasta_files = irma_pneumo.irma_fasta,
            db = "OrthopneumovirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_pneumo {
            input:
            sample_name = sample_name,
            module_irma = "PNEUMO",
            blast_res = blast_pneumo.blast_res,
            irma_type = irma_qc_pneumo.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_boca {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "BOCA",
        docker = "irma:230723"
    }

    if (irma_boca.proceed == "yes") {

        call irma.IrmaQC as irma_qc_boca {
            input:
            sample_name = sample_name,
            module_irma = "BOCA",
            irma_qc = irma_boca.irma_qc,
            reference_fasta = reference_fasta_boca,
            irma_fasta = irma_boca.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_boca {
            input:
            sample_name = sample_name,
            module_irma = "BOCA",
            fasta_files = irma_boca.irma_fasta,
            db = "BocaparvovirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_boca {
            input:
            sample_name = sample_name,
            module_irma = "BOCA",
            blast_res = blast_boca.blast_res,
            irma_type = irma_qc_boca.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_corona {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "CORONA",
        docker = "irma:230723"
    }

    if (irma_corona.proceed == "yes") {

        call irma.IrmaQC as irma_qc_corona {
            input:
            sample_name = sample_name,
            module_irma = "CORONA",
            irma_qc = irma_corona.irma_qc,
            reference_fasta = reference_fasta_corona,
            irma_fasta = irma_corona.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_corona {
            input:
            sample_name = sample_name,
            module_irma = "CORONA",
            fasta_files = irma_corona.irma_fasta,
            db = "CoronavirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_corona {
            input:
            sample_name = sample_name,
            module_irma = "CORONA",
            blast_res = blast_corona.blast_res,
            irma_type = irma_qc_corona.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_respiro {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "RESPIRO",
        docker = "irma:230723"
    }

    if (irma_respiro.proceed == "yes") {

        call irma.IrmaQC as irma_qc_respiro {
            input:
            sample_name = sample_name,
            module_irma = "RESPIRO",
            irma_qc = irma_respiro.irma_qc,
            reference_fasta = reference_fasta_respiro,
            irma_fasta = irma_respiro.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_respiro {
            input:
            sample_name = sample_name,
            module_irma = "RESPIRO",
            fasta_files = irma_respiro.irma_fasta,
            db = "RespirovirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_respiro {
            input:
            sample_name = sample_name,
            module_irma = "RESPIRO",
            blast_res = blast_respiro.blast_res,
            irma_type = irma_qc_respiro.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_rhino {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "RHINO",
        docker = "irma:230723"
    }

    if (irma_rhino.proceed == "yes") {

        call irma.IrmaQC as irma_qc_rhino {
            input:
            sample_name = sample_name,
            module_irma = "RHINO",
            irma_qc = irma_rhino.irma_qc,
            reference_fasta = reference_fasta_rhino,
            irma_fasta = irma_rhino.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_rhino {
            input:
            sample_name = sample_name,
            module_irma = "RHINO",
            fasta_files = irma_rhino.irma_fasta,
            db = "RhinovirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_rhino {
            input:
            sample_name = sample_name,
            module_irma = "RHINO",
            blast_res = blast_rhino.blast_res,
            irma_type = irma_qc_rhino.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_rubula {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "RUBULA",
        docker = "irma:230723"
    }

    if (irma_rubula.proceed == "yes") {

        call irma.IrmaQC as irma_qc_rubula {
            input:
            sample_name = sample_name,
            module_irma = "RUBULA",
            irma_qc = irma_rubula.irma_qc,
            reference_fasta = reference_fasta_rubula,
            irma_fasta = irma_rubula.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_rubula {
            input:
            sample_name = sample_name,
            module_irma = "RUBULA",
            fasta_files = irma_rubula.irma_fasta,
            db = "OrthorubulavirusDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_rubula {
            input:
            sample_name = sample_name,
            module_irma = "RUBULA",
            blast_res = blast_rubula.blast_res,
            irma_type = irma_qc_rubula.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_metapneumo {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "METAPMEUMO",
        docker = "irma:230723"
    }

    if (irma_metapneumo.proceed == "yes") {

        call irma.IrmaQC as irma_qc_metapneumo {
            input:
            sample_name = sample_name,
            module_irma = "METAPMEUMO",
            irma_qc = irma_metapneumo.irma_qc,
            reference_fasta = reference_fasta_metapneumo,
            irma_fasta = irma_metapneumo.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_metapneumo {
            input:
            sample_name = sample_name,
            module_irma = "METAPMEUMO",
            fasta_files = irma_metapneumo.irma_fasta,
            db = "MetapneumoDB",
            threads = threads,
            docker = "blast:230723"
        }

        call irma.Report as report_metapneumo {
            input:
            sample_name = sample_name,
            module_irma = "METAPMEUMO",
            blast_res = blast_metapneumo.blast_res,
            irma_type = irma_qc_metapneumo.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.Irma as irma_cov {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "CoV",
        docker = "irma:230723"
    }

    if (irma_cov.proceed == "yes") {

        call irma.IrmaQC as irma_qc_cov {
            input:
            sample_name = sample_name,
            module_irma = "CoV",
            irma_qc = irma_cov.irma_qc,
            reference_fasta = reference_fasta_cov,
            irma_fasta = irma_cov.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Report as report_cov {
            input:
            sample_name = sample_name,
            module_irma = "CoV",
            irma_type = irma_qc_cov.irma_type,
            docker = "python:4.4"
        }
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
        File? measles_qc_json = irma_qc_measles.irma_qc_json
        File? measles_json = report_measles.report
        File? adeno_qc_json = irma_qc_adeno.irma_qc_json
        File? adeno_json = report_adeno.report
        File? pneumo_qc_json = irma_qc_pneumo.irma_qc_json
        File? pneumo_json = report_pneumo.report
        File? boca_qc_json = irma_qc_boca.irma_qc_json
        File? boca_json = report_boca.report
        File? corona_qc_json = irma_qc_corona.irma_qc_json
        File? corona_json = report_corona.report
        File? respiro_qc_json = irma_qc_respiro.irma_qc_json
        File? respiro_json = report_respiro.report
        File? rhino_qc_json = irma_qc_rhino.irma_qc_json
        File? rhino_json = report_rhino.report
        File? rubula_qc_json = irma_qc_rubula.irma_qc_json
        File? rubula_json = report_rubula.report
        File? metapneumo_qc_json = irma_qc_metapneumo.irma_qc_json
        File? metapneumo_json = report_metapneumo.report
        File? cov_qc_json = irma_qc_cov.irma_qc_json
        File? cov_json = report_cov.report
        File? flu_qc_json = irma_qc_flu.irma_qc_json
        File? flu_json = report_flu.report
        File kraken_txt = kraken2.report_txt
        File bracken_txt = bracken.report_txt
        File krona_kraken_html = krona_kraken.report_html
        File kraken_virus_txt = kraken2_vir.report_txt
        File bracken_virus_txt = bracken_vir.report_txt
    }
}
