version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run wdl_scripts/GutGenotyping.wdl -i inputs.json -o options.json

import "./common_tasks/preprocessing.wdl" as preprocessing
import "./common_tasks/kraken2.wdl" as kraken2
import "./common_tasks/irma.wdl" as irma


workflow GutGenotyping {

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
        File reference_fasta_entero = "reference_fasta/entero.fasta"
        File reference_fasta_adeno = "reference_fasta/adeno.fasta"
        File reference_fasta_astrovirus = "reference_fasta/astrovirus.fasta"
        File reference_fasta_caliciviridae = "reference_fasta/caliciviridae.fasta"
        File reference_fasta_cosavirus = "reference_fasta/cosavirus.fasta"
        File reference_fasta_rotavirus = "reference_fasta/rotavirus.fasta"
        File reference_fasta_salivirus = "reference_fasta/salivirus.fasta"
        File reference_fasta_torovirus = "reference_fasta/torovirus.fasta"
        File rotavirus_metadata = "rotavirus_metadata.csv"
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
        docker = "nanozoo/bracken:2.8--dcb3e47"
    }


    call irma.IrmaGut as irma_entero {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ENTERO",
        docker = "irma:v1.1.3"
    }

    if (irma_entero.proceed == "yes") {

        call irma.IrmaQC as irma_qc_entero {
            input:
            sample_name = sample_name,
            module_irma = "ENTERO",
            irma_qc = irma_entero.irma_qc,
            reference_fasta = reference_fasta_entero,
            irma_fasta = irma_entero.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_entero {
            input:
            sample_name = sample_name,
            module_irma = "ENTERO",
            fasta_files = irma_entero.irma_fasta,
            db = "EnteroDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_entero {
            input:
            sample_name = sample_name,
            module_irma = "ENTERO",
            blast_res = blast_entero.blast_res,
            irma_type = irma_qc_entero.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_salivirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "SALIVIRUS",
        docker = "irma:v1.1.3"
    }

    if (irma_salivirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_salivirus {
            input:
            sample_name = sample_name,
            module_irma = "SALIVIRUS",
            irma_qc = irma_salivirus.irma_qc,
            reference_fasta = reference_fasta_salivirus,
            irma_fasta = irma_salivirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_salivirus {
            input:
            sample_name = sample_name,
            module_irma = "SALIVIRUS",
            fasta_files = irma_salivirus.irma_fasta,
            db = "SalivirusDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_salivirus {
            input:
            sample_name = sample_name,
            module_irma = "SALIVIRUS",
            blast_res = blast_salivirus.blast_res,
            irma_type = irma_qc_salivirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_astrovirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ASTROVIRUS",
        docker = "irma:v1.1.3"
    }

    if (irma_astrovirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_astrovirus {
            input:
            sample_name = sample_name,
            module_irma = "ASTROVIRUS",
            irma_qc = irma_astrovirus.irma_qc,
            reference_fasta = reference_fasta_astrovirus,
            irma_fasta = irma_astrovirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_astrovirus {
            input:
            sample_name = sample_name,
            module_irma = "ASTROVIRUS",
            fasta_files = irma_astrovirus.irma_fasta,
            db = "AstrovirusDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_astrovirus {
            input:
            sample_name = sample_name,
            module_irma = "ASTROVIRUS",
            blast_res = blast_astrovirus.blast_res,
            irma_type = irma_qc_astrovirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_caliciviridae {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "CALICIVIRIDAE",
        docker = "irma:v1.1.3"
    }

    if (irma_caliciviridae.proceed == "yes") {

        call irma.IrmaQC as irma_qc_caliciviridae {
            input:
            sample_name = sample_name,
            module_irma = "CALICIVIRIDAE",
            irma_qc = irma_caliciviridae.irma_qc,
            reference_fasta = reference_fasta_caliciviridae,
            irma_fasta = irma_caliciviridae.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_caliciviridae {
            input:
            sample_name = sample_name,
            module_irma = "CALICIVIRIDAE",
            fasta_files = irma_caliciviridae.irma_fasta,
            db = "CaliciviridaeDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_caliciviridae {
            input:
            sample_name = sample_name,
            module_irma = "CALICIVIRIDAE",
            blast_res = blast_caliciviridae.blast_res,
            irma_type = irma_qc_caliciviridae.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_cosavirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "COSAVIRUS",
        docker = "irma:v1.1.3"
    }

    if (irma_cosavirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_cosavirus {
            input:
            sample_name = sample_name,
            module_irma = "COSAVIRUS",
            irma_qc = irma_cosavirus.irma_qc,
            reference_fasta = reference_fasta_cosavirus,
            irma_fasta = irma_cosavirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_cosavirus {
            input:
            sample_name = sample_name,
            module_irma = "COSAVIRUS",
            fasta_files = irma_cosavirus.irma_fasta,
            db = "CosavirusDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_cosavirus {
            input:
            sample_name = sample_name,
            module_irma = "COSAVIRUS",
            blast_res = blast_cosavirus.blast_res,
            irma_type = irma_qc_cosavirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_torovirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "TOROVIRUS",
        docker = "irma:v1.1.3"
    }

    if (irma_torovirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_torovirus {
            input:
            sample_name = sample_name,
            module_irma = "TOROVIRUS",
            irma_qc = irma_torovirus.irma_qc,
            reference_fasta = reference_fasta_torovirus,
            irma_fasta = irma_torovirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_torovirus {
            input:
            sample_name = sample_name,
            module_irma = "TOROVIRUS",
            fasta_files = irma_torovirus.irma_fasta,
            db = "TorovirusDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_torovirus {
            input:
            sample_name = sample_name,
            module_irma = "TOROVIRUS",
            blast_res = blast_torovirus.blast_res,
            irma_type = irma_qc_torovirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_rotavirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ROTAVIRUS",
        docker = "irma:v1.1.3"
    }

    if (irma_rotavirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_rotavirus {
            input:
            sample_name = sample_name,
            module_irma = "ROTAVIRUS",
            irma_qc = irma_rotavirus.irma_qc,
            reference_fasta = reference_fasta_rotavirus,
            irma_fasta = irma_rotavirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_rotavirus {
            input:
            sample_name = sample_name,
            module_irma = "ROTAVIRUS",
            fasta_files = irma_rotavirus.irma_fasta,
            db = "RotavirusDB",
            threads = threads,
            docker = "blast:2.15.0"
        }

        call irma.Report as report_rotavirus {
            input:
            sample_name = sample_name,
            module_irma = "ROTAVIRUS",
            blast_res = blast_rotavirus.blast_res,
            irma_type = irma_rotavirus.rotavirus_type,
            rotavirus_metadata = rotavirus_metadata,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_adeno {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ADENO",
        docker = "irma:v1.1.3"
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
            docker = "blast:2.15.0"
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


    output {
        File fastqc_row_R1_html = fastqc_row_R1.summary_html
        File fastqc_row_R2_html = fastqc_row_R2.summary_html
        File fastqc_trimed_R1_html = fastqc_trimed_R1.summary_html
        File fastqc_trimed_R2_html = fastqc_trimed_R2.summary_html
        File preprocessing_qc_json = preprocessing_qc.report_json
        File? entero_qc_json = irma_qc_entero.irma_qc_json
        File? entero_json = report_entero.report
        File? adeno_qc_json = irma_qc_adeno.irma_qc_json
        File? adeno_json = report_adeno.report
        File? salivirus_qc_json = irma_qc_salivirus.irma_qc_json
        File? salivirus_json = report_salivirus.report
        File? astrovirus_qc_json = irma_qc_astrovirus.irma_qc_json
        File? astrovirus_json = report_astrovirus.report
        File? caliciviridae_qc_json = irma_qc_caliciviridae.irma_qc_json
        File? caliciviridae_json = report_caliciviridae.report
        File? cosavirus_qc_json = irma_qc_cosavirus.irma_qc_json
        File? cosavirus_json = report_cosavirus.report
        File? torovirus_qc_json = irma_qc_torovirus.irma_qc_json
        File? torovirus_json = report_torovirus.report
        File? rotavirus_qc_json = irma_qc_rotavirus.irma_qc_json
        File? rotavirus_json = report_rotavirus.report
        File kraken_txt = kraken2.report_txt
        File bracken_txt = bracken.report_txt
        File krona_kraken_html = krona_kraken.report_html
        File kraken_virus_txt = kraken2_vir.report_txt
        File bracken_virus_txt = bracken_vir.report_txt
    }
}
