version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run run.wdl -i inputs.json

import "../tasks/preprocessing.wdl" as preprocessing
import "../tasks/kraken2.wdl" as kraken2
import "../tasks/irma.wdl" as irma

workflow GutGenotyping {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        Boolean cut_primers = false
        File? primer_left = "/home/admin/cromwell_all/cut_primers/left_primers.fasta"
        File? primer_right = "/home/admin/cromwell_all/cut_primers/right_primers.fasta"
        Boolean transcriptome_filtering = true # false will take index_genome, true - index_transcriptome
        File index_genome = "/home/admin/db/human_bowtie2_index/GRCh38_ERCC.bowtie2.tar"
        File index_transcriptome = "/home/admin/db/human_bowtie2_index/GRCh38_transcriptome.bowtie2.tar"
        File kraken2_standard_8gb = "/home/admin/db/kraken2/k2_standard_08gb_20231009.tar.gz"
        File kraken2_virus = "/home/admin/db/kraken2/k2_viral_20231009.tar.gz"
        String kraken_level = "S" #(U)nclassified, (R)oot, (D)omain, (K)ingdom (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
        Int threads = 8

        File reference_fasta_entero = "/home/admin/cromwell_gut/refseq/entero.fasta"
        File reference_fasta_astrovirus = "/home/admin/cromwell_gut/refseq/astrovirus.fasta"
        File reference_fasta_kobuvirus = "/home/admin/cromwell_gut/refseq/kobuvirus.fasta"
        File reference_fasta_mastadeno = "/home/admin/cromwell_gut/refseq/mastadeno.fasta"
        File reference_fasta_norovirus = "/home/admin/cromwell_gut/refseq/norovirus.fasta"
        File reference_fasta_parechovirus = "/home/admin/cromwell_gut/refseq/parechovirus.fasta"
        File reference_fasta_picobina = "/home/admin/cromwell_gut/refseq/picobina.fasta"
        File reference_fasta_rosavirus = "/home/admin/cromwell_gut/refseq/rosavirus.fasta"
        File reference_fasta_rotavirus = "/home/admin/cromwell_gut/refseq/rotavirus.fasta"
        File reference_fasta_salivirus = "/home/admin/cromwell_gut/refseq/salivirus.fasta"
        File reference_fasta_sapovirus = "/home/admin/cromwell_gut/refseq/sapovirus.fasta"
        File reference_fasta_torovirus = "/home/admin/cromwell_gut/refseq/torovirus.fasta"
        File reference_fasta_cosavirus = "/home/admin/cromwell_gut/refseq/cosavirus.fasta"

        File rotavirus_metadata = "/home/admin/cromwell_gut/make_metadata_rotavirus/rotavirus_metadata.csv"
    }


    #call preprocessing.FastQC as fastqc_row_R1 {
    #    input:
    #    fastq = fastq_1,
    #    docker = "fastqc:4.4"
    #}

    #call preprocessing.FastQC as fastqc_row_R2 {
    #    input:
    #    fastq = fastq_2,
    #    docker = "fastqc:4.4"
    #}

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

    #call preprocessing.FastQC as fastqc_trimed_R1 {
    #    input:
    #    fastq = host_filter.host_filtered_fastq_1,
    #    docker = "fastqc:4.4"
    #}

    #call preprocessing.FastQC as fastqc_trimed_R2 {
    #    input:
    #    fastq = host_filter.host_filtered_fastq_2,
    #    docker = "fastqc:4.4"
    #}

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

    #call kraken2.Krona as krona_kraken {
    #    input:
    #    sample_name = sample_name,
    #    report = kraken2.report_txt,
    #    docker = "krona:2.8.1"
    #}

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
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
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
            docker = "blast:2.15.0gut"
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
        module_irma = "SALI",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_salivirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_salivirus {
            input:
            sample_name = sample_name,
            module_irma = "SALI",
            irma_qc = irma_salivirus.irma_qc,
            reference_fasta = reference_fasta_salivirus,
            irma_fasta = irma_salivirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_salivirus {
            input:
            sample_name = sample_name,
            module_irma = "SALI",
            fasta_files = irma_salivirus.irma_fasta,
            db = "SaliDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_salivirus {
            input:
            sample_name = sample_name,
            module_irma = "SALI",
            blast_res = blast_salivirus.blast_res,
            irma_type = irma_qc_salivirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_astrovirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ASTRO",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_astrovirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_astrovirus {
            input:
            sample_name = sample_name,
            module_irma = "ASTRO",
            irma_qc = irma_astrovirus.irma_qc,
            reference_fasta = reference_fasta_astrovirus,
            irma_fasta = irma_astrovirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_astrovirus {
            input:
            sample_name = sample_name,
            module_irma = "ASTRO",
            fasta_files = irma_astrovirus.irma_fasta,
            db = "AstroDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_astrovirus {
            input:
            sample_name = sample_name,
            module_irma = "ASTRO",
            blast_res = blast_astrovirus.blast_res,
            irma_type = irma_qc_astrovirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_kobuvirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "KOBU",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_kobuvirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_kobuvirus {
            input:
            sample_name = sample_name,
            module_irma = "KOBU",
            irma_qc = irma_kobuvirus.irma_qc,
            reference_fasta = reference_fasta_kobuvirus,
            irma_fasta = irma_kobuvirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_kobuvirus {
            input:
            sample_name = sample_name,
            module_irma = "KOBU",
            fasta_files = irma_kobuvirus.irma_fasta,
            db = "KobuDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_kobuvirus {
            input:
            sample_name = sample_name,
            module_irma = "KOBU",
            blast_res = blast_kobuvirus.blast_res,
            irma_type = irma_qc_kobuvirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_cosavirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "COSA",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_cosavirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_cosavirus {
            input:
            sample_name = sample_name,
            module_irma = "COSA",
            irma_qc = irma_cosavirus.irma_qc,
            reference_fasta = reference_fasta_cosavirus,
            irma_fasta = irma_cosavirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_cosavirus {
            input:
            sample_name = sample_name,
            module_irma = "COSA",
            fasta_files = irma_cosavirus.irma_fasta,
            db = "CosaDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_cosavirus {
            input:
            sample_name = sample_name,
            module_irma = "COSA",
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
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
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
            docker = "blast:2.15.0gut"
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
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
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
            docker = "blast:2.15.0gut"
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

    call irma.IrmaGut as irma_mastadeno {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "MASTADENO",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_mastadeno.proceed == "yes") {

        call irma.IrmaQC as irma_qc_mastadeno {
            input:
            sample_name = sample_name,
            module_irma = "MASTADENO",
            irma_qc = irma_mastadeno.irma_qc,
            reference_fasta = reference_fasta_mastadeno,
            irma_fasta = irma_mastadeno.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_mastadeno {
            input:
            sample_name = sample_name,
            module_irma = "MASTADENO",
            fasta_files = irma_mastadeno.irma_fasta,
            db = "MastadenoDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_mastadeno {
            input:
            sample_name = sample_name,
            module_irma = "MASTADENO",
            blast_res = blast_mastadeno.blast_res,
            irma_type = irma_qc_mastadeno.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_norovirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "NOROVIRUS",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_norovirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_norovirus {
            input:
            sample_name = sample_name,
            module_irma = "NOROVIRUS",
            irma_qc = irma_norovirus.irma_qc,
            reference_fasta = reference_fasta_norovirus,
            irma_fasta = irma_norovirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_norovirus {
            input:
            sample_name = sample_name,
            module_irma = "NOROVIRUS",
            fasta_files = irma_norovirus.irma_fasta,
            db = "NorovirusDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_norovirus {
            input:
            sample_name = sample_name,
            module_irma = "NOROVIRUS",
            blast_res = blast_norovirus.blast_res,
            irma_type = irma_qc_norovirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_parechovirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "PARECHO",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_parechovirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_parechovirus {
            input:
            sample_name = sample_name,
            module_irma = "PARECHO",
            irma_qc = irma_parechovirus.irma_qc,
            reference_fasta = reference_fasta_parechovirus,
            irma_fasta = irma_parechovirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_parechovirus {
            input:
            sample_name = sample_name,
            module_irma = "PARECHO",
            fasta_files = irma_parechovirus.irma_fasta,
            db = "ParechoDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_parechovirus {
            input:
            sample_name = sample_name,
            module_irma = "PARECHO",
            blast_res = blast_parechovirus.blast_res,
            irma_type = irma_qc_parechovirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_picobina {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "PICOBINA",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_picobina.proceed == "yes") {

        call irma.IrmaQC as irma_qc_picobina {
            input:
            sample_name = sample_name,
            module_irma = "PICOBINA",
            irma_qc = irma_picobina.irma_qc,
            reference_fasta = reference_fasta_picobina,
            irma_fasta = irma_picobina.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_picobina {
            input:
            sample_name = sample_name,
            module_irma = "PICOBINA",
            fasta_files = irma_picobina.irma_fasta,
            db = "PicobinaDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_picobina {
            input:
            sample_name = sample_name,
            module_irma = "PICOBINA",
            blast_res = blast_picobina.blast_res,
            irma_type = irma_qc_picobina.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_rosavirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "ROSA",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_rosavirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_rosavirus {
            input:
            sample_name = sample_name,
            module_irma = "ROSA",
            irma_qc = irma_rosavirus.irma_qc,
            reference_fasta = reference_fasta_rosavirus,
            irma_fasta = irma_rosavirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_rosavirus {
            input:
            sample_name = sample_name,
            module_irma = "ROSA",
            fasta_files = irma_rosavirus.irma_fasta,
            db = "RosaDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_rosavirus {
            input:
            sample_name = sample_name,
            module_irma = "ROSA",
            blast_res = blast_rosavirus.blast_res,
            irma_type = irma_qc_rosavirus.irma_type,
            docker = "python:4.4"
        }
    }

    call irma.IrmaGut as irma_sapovirus {
        input:
        trim_R1 = host_filter.host_filtered_fastq_1,
        trim_R2 = host_filter.host_filtered_fastq_2,
        module_irma = "SAPOVIRUS",
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/irma:0.6.4gut"
    }

    if (irma_sapovirus.proceed == "yes") {

        call irma.IrmaQC as irma_qc_sapovirus {
            input:
            sample_name = sample_name,
            module_irma = "SAPOVIRUS",
            irma_qc = irma_sapovirus.irma_qc,
            reference_fasta = reference_fasta_sapovirus,
            irma_fasta = irma_sapovirus.irma_fasta,
            docker = "python:4.4"
        }

        call irma.Blast as blast_sapovirus {
            input:
            sample_name = sample_name,
            module_irma = "SAPOVIRUS",
            fasta_files = irma_sapovirus.irma_fasta,
            db = "SapovirusDB",
            threads = threads,
            docker = "blast:2.15.0gut"
        }

        call irma.Report as report_sapovirus {
            input:
            sample_name = sample_name,
            module_irma = "SAPOVIRUS",
            blast_res = blast_sapovirus.blast_res,
            irma_type = irma_qc_sapovirus.irma_type,
            docker = "python:4.4"
        }
    }


    output {
        #File fastqc_row_R1_html = fastqc_row_R1.summary_html
        #File fastqc_row_R2_html = fastqc_row_R2.summary_html
        #File fastqc_trimed_R1_html = fastqc_trimed_R1.summary_html
        #File fastqc_trimed_R2_html = fastqc_trimed_R2.summary_html
        File preprocessing_qc_json = preprocessing_qc.report_json

        File? entero_qc_json = irma_qc_entero.irma_qc_json
        File? entero_json = report_entero.report

        File? adeno_qc_json = irma_qc_mastadeno.irma_qc_json
        File? adeno_json = report_mastadeno.report

        File? salivirus_qc_json = irma_qc_salivirus.irma_qc_json
        File? salivirus_json = report_salivirus.report

        File? astrovirus_qc_json = irma_qc_astrovirus.irma_qc_json
        File? astrovirus_json = report_astrovirus.report

        File? rosavirus_qc_json = irma_qc_rosavirus.irma_qc_json
        File? rosavirus_json = report_rosavirus.report

        File? cosavirus_qc_json = irma_qc_cosavirus.irma_qc_json
        File? cosavirus_json = report_cosavirus.report

        File? torovirus_qc_json = irma_qc_torovirus.irma_qc_json
        File? torovirus_json = report_torovirus.report

        File? rotavirus_qc_json = irma_qc_rotavirus.irma_qc_json
        File? rotavirus_json = report_rotavirus.report

        File? kobuvirus_qc_json = irma_qc_kobuvirus.irma_qc_json
        File? kobuvirus_json = report_kobuvirus.report

        File? norovirus_qc_json = irma_qc_norovirus.irma_qc_json
        File? norovirus_json = report_norovirus.report

        File? parechovirus_qc_json = irma_qc_parechovirus.irma_qc_json
        File? parechovirus_json = report_parechovirus.report

        File? picobina_qc_json = irma_qc_picobina.irma_qc_json
        File? picobina_json = report_picobina.report

        File? sapovirus_qc_json = irma_qc_sapovirus.irma_qc_json
        File? sapovirus_json = report_sapovirus.report

        File kraken_txt = kraken2.report_txt
        File bracken_txt = bracken.report_txt
        #File krona_kraken_html = krona_kraken.report_html
        File kraken_virus_txt = kraken2_vir.report_txt
        File bracken_virus_txt = bracken_vir.report_txt
    }
}

