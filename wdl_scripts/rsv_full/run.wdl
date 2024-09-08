version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run VirGenotyping.wdl -i inputs.json -o options.json

import "preprocessing.wdl" as preprocessing
import "kraken2.wdl" as kraken2
import "irma.wdl" as irma
import "nextclade.wdl" as nextclade
import "yandex_utilities.wdl" as Utils
import "summary_report.wdl" as summary_report
import "sarscov2.wdl" as sarscov2
#import "../common_tasks/preprocessing.wdl" as preprocessing
#import "../common_tasks/kraken2.wdl" as kraken2
#import "../common_tasks/irma.wdl" as irma
#import "../common_tasks/nextclade.wdl" as nextclade
#import "../common_tasks/yandex_utilities.wdl" as Utils

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
        Boolean cut_primers = false
        File? primer_left = "/home/cromwell/mycob-ref/rsv_full/left_primers.fasta"
        File? primer_right = "/home/cromwell/mycob-ref/rsv_full/right_primers.fasta"
        Boolean transcriptome_filtering = true # false will take index_genome, true - index_transcriptome
        File index_genome = "/home/cromwell/mycob-ref/rsv_full/human_bowtie2_index/GRCh38_ERCC.bowtie2.tar"
        File index_transcriptome = "/home/cromwell/mycob-ref/rsv_full/human_bowtie2_index/GRCh38_transcriptome.bowtie2.tar"
        File kraken2_standard_8gb = "/home/cromwell/mycob-ref/rsv_full/kraken2/k2_standard_08gb_20231009.tar.gz"
        File kraken2_virus = "/home/cromwell/mycob-ref/rsv_full/kraken2/k2_viral_20231009.tar.gz"
        String kraken_level = "S" #(U)nclassified, (R)oot, (D)omain, (K)ingdom (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
        Int threads = 8
        File reference_fasta_measles = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/MEASLES/reference/consensus.fasta"
        File reference_fasta_adeno = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/ADENO/reference/consensus.fasta"
        File reference_fasta_pneumo = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/PNEUMO/reference/consensus.fasta"
        File reference_fasta_boca = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/BOCA/reference/consensus.fasta"
        File reference_fasta_corona = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/CORONA/reference/consensus.fasta"
        File reference_fasta_respiro = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/RESPIRO/reference/consensus.fasta"
        File reference_fasta_rhino = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/RHINO/reference/consensus.fasta"
        File reference_fasta_rubula = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/RUBULA/reference/consensus.fasta"
        File reference_fasta_cov = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/CoV/reference/consensus.fasta"
        File reference_fasta_metapneumo = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/METAPMEUMO/reference/consensus.fasta"
        File reference_fasta_flu = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/consensus.fasta"
        File HA_ref = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/HA_wis_67_2022.fasta"
        File NA_ref = "/home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/NA_wis_67_2022.fasta"
        File snpeff_config = "snpEff_5.2/snpEff.config"
        File snpeff_db = "snpEff_5.2/data.zip"        
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

    call preprocessing.FastQC as fastqc_row_R1 {
        input:
        fastq = Files[0][0],
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call preprocessing.FastQC as fastqc_row_R2 {
        input:
        fastq = Files[0][1],
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }


    call preprocessing.Trimmomatic as trimmomatic {
        input:
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        sample_name = sample_name,
        threads = threads,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/trimmomatic:0.39"
    }


    if (cut_primers) {
        call preprocessing.Cutadapt as cutadapt {
            input:
            fastq_1 = trimmomatic.trim_fastq_1,
            fastq_2 = trimmomatic.trim_fastq_2,
            sample_name = sample_name,
            primer_left = primer_left,
            primer_right = primer_right,
            docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/cutadapt:4.4"
        }
    }

    call preprocessing.HostFilter as host_filter {
        input:
        fastq_1 = if cut_primers then cutadapt.cut_fastq_1 else trimmomatic.trim_fastq_1,
        fastq_2 = if cut_primers then cutadapt.cut_fastq_2 else trimmomatic.trim_fastq_2,
        sample_name = sample_name,
        index_tar = if transcriptome_filtering then index_transcriptome else index_genome,
        threads = threads,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/bowtie2:2.5.1"
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
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/python:3"
    }

    call preprocessing.FastQC as fastqc_trimed_R1 {
        input:
        fastq = host_filter.host_filtered_fastq_1,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call preprocessing.FastQC as fastqc_trimed_R2 {
        input:
        fastq = host_filter.host_filtered_fastq_2,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/fastqc:0.12.0"
    }

    call kraken2.Kraken2 as kraken2 {
        input:
        fastq_1 = host_filter.host_filtered_fastq_1,
        fastq_2 = host_filter.host_filtered_fastq_2,
        sample_name = sample_name,
        kraken2_classifier = kraken2_standard_8gb,
        threads = 1,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/kraken2:2.1.3"
    }

    call kraken2.Bracken as bracken {
        input:
        sample_name = sample_name,
        kraken_report = kraken2.report_txt,
        kraken2_classifier = kraken2_standard_8gb,
        level = kraken_level,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/bracken:2.8--dcb3e47"
    }

    #call kraken2.Krona as krona_kraken {
    #    input:
    #    sample_name = sample_name,
    #    report = kraken2.report_txt,
    #    docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/krona:2.8.1"
    #}

    call kraken2.Kraken2 as kraken2_vir {
        input:
        fastq_1 = host_filter.host_filtered_fastq_1,
        fastq_2 = host_filter.host_filtered_fastq_2,
        sample_name = sample_name,
        kraken2_classifier = kraken2_virus,
        threads = 1,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/kraken2:2.1.3"
    }

    call kraken2.Bracken as bracken_vir {
        input:
        sample_name = sample_name,
        kraken_report = kraken2_vir.report_txt,
        kraken2_classifier = kraken2_virus,
        level = kraken_level,
        docker = "cr.yandex/crpl2lv1lkr7g21e6q8g/bracken:2.8--dcb3e47"
    }


    output {
    String fastq_R1 = fastq_1
    String fastq_R2 = fastq_2
    }    
}
