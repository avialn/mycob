version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run wdl_scripts/FluGenotyping.wdl -i inputs.json -o options.json

import "./common_tasks/preprocessing.wdl" as preprocessing
import "./common_tasks/kraken2.wdl" as kraken2
import "./common_tasks/irma.wdl" as irma
import "./common_tasks/nextclade.wdl" as nextclade
import "./common_tasks/snpeff.wdl" as snpeff
import "./common_tasks/summary_report.wdl" as summary_report

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
        File HA_ref = "reference_fasta/HA_wis_67_2022.fasta"
        File NA_ref = "reference_fasta/NA_wis_67_2022.fasta"
        File snpeff_config = "snpEff_5.2/snpEff.config"
        File snpeff_db = "snpEff_5.2/data.zip"
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

    call snpeff.Minimap2 as minimap_flu_ha {
        input:
        fastq_1 = host_filter.host_filtered_fastq_1,
        fastq_2 = host_filter.host_filtered_fastq_2,
        ref = HA_ref,
        threads = threads,
        docker = "staphb/minimap2:latest"
    }

    if (minimap_flu_ha.proceed == "yes") {
        call snpeff.Samtools as samtools_flu_ha {
            input:
            ref = HA_ref,
            sam = minimap_flu_ha.file_sam,
            docker = "samtools_bcftoolss:1.19"
        }

        call snpeff.Snpeff as snpeff_ha {
            input:
            snpeff_config = snpeff_config,
            snpeff_db = snpeff_db,
            ref = HA_ref,
            vcf_file = samtools_flu_ha.filtered_vsf,
            docker = "snpeff:5.2"
        }

        call snpeff.ParseSnpeff as parse_snpeff_ha {
            input:
            snpeff_csv = snpeff_ha.snpeff_csv,
            ref_name = snpeff_ha.ref_name,
            docker = "python:4.4"
        }

    }

    call snpeff.MinimapQC as minimap_qc_ha {
        input:
        sample_name = sample_name,
        ref_name = minimap_flu_ha.ref_name,
        total_count = minimap_flu_ha.total_count,
        matched_count = minimap_flu_ha.matched_count,
        alignment_done = minimap_flu_ha.proceed,
        ref_length_bp = samtools_flu_ha.ref_length_bp,
        coverage_bp = samtools_flu_ha.coverage_bp,
        coverage_percentage = samtools_flu_ha.coverage_percentage,
        docker = "python:4.4"
    }

    call snpeff.Minimap2 as minimap_flu_na {
        input:
        fastq_1 = host_filter.host_filtered_fastq_1,
        fastq_2 = host_filter.host_filtered_fastq_2,
        ref = NA_ref,
        threads = threads,
        docker = "staphb/minimap2:latest"
    }

    if (minimap_flu_na.proceed == "yes") {
        call snpeff.Samtools as samtools_flu_na {
            input:
            ref = NA_ref,
            sam = minimap_flu_na.file_sam,
            docker = "samtools_bcftoolss:1.19"
        }

        call snpeff.Snpeff as snpeff_na {
            input:
            snpeff_config = snpeff_config,
            snpeff_db = snpeff_db,
            ref = NA_ref,
            vcf_file = samtools_flu_na.filtered_vsf,
            docker = "snpeff:5.2"
        }

        call snpeff.ParseSnpeff as parse_snpeff_na {
            input:
            snpeff_csv = snpeff_na.snpeff_csv,
            ref_name = snpeff_na.ref_name,
            docker = "python:4.4"
        }
    }

    call snpeff.MinimapQC as minimap_qc_na {
        input:
        sample_name = sample_name,
        ref_name = minimap_flu_na.ref_name,
        total_count = minimap_flu_na.total_count,
        matched_count = minimap_flu_na.matched_count,
        alignment_done = minimap_flu_na.proceed,
        ref_length_bp = samtools_flu_na.ref_length_bp,
        coverage_bp = samtools_flu_na.coverage_bp,
        coverage_percentage = samtools_flu_na.coverage_percentage,
        docker = "python:4.4"
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

        call nextclade.Nextclade as nextclade_flu_ha {
            input:
            fasta = irma_flu.ha_fasta,
            ref_name = "HA_wis_67_2022",
            docker = "nextclade:3.3"
        }

        if (nextclade_flu_ha.proceed == "yes") {

            call nextclade.NextcladeParse as nextclade_parse_flu_ha {
                input:
                nextclade_tsv = nextclade_flu_ha.nextclade_tsv,
                search_antigenic_mut = "yes",
                ref_name = "HA_wis_67_2022",
                docker = "python:4.4"
            }
        }

        call nextclade.Nextclade as nextclade_flu_na {
            input:
            fasta = irma_flu.na_fasta,
            ref_name = "NA_wis_67_2022",
            docker = "nextclade:3.3"
        }

        if (nextclade_flu_na.proceed == "yes") {
            call nextclade.NextcladeParse as nextclade_parse_flu_na {
                input:
                nextclade_tsv = nextclade_flu_na.nextclade_tsv,
                search_antigenic_mut = "no",
                ref_name = "NA_wis_67_2022",
                docker = "python:4.4"
            }
        }
    }

    if (irma_flu.proceed == "yes" || minimap_flu_ha.proceed == "yes" || minimap_flu_na.proceed == "yes") {
        call summary_report.SummaryReport as summary_report_flu {
            input:
            sample_name = sample_name,
            irma_report = report_flu.report,
            snpeff_HA_report = parse_snpeff_ha.snpeff_json,
            snpeff_NA_report = parse_snpeff_na.snpeff_json,
            nextclade_HA_report = nextclade_parse_flu_ha.nextclade_json,
            nextclade_NA_report = nextclade_parse_flu_na.nextclade_json,
            docker = "python:4.4"
        }
    }

    output {

        # preprocessing: fastqc, trimmomatic, cutadapt, host-filtering
        File fastqc_row_R1_html = fastqc_row_R1.summary_html
        File fastqc_row_R2_html = fastqc_row_R2.summary_html
        File fastqc_trimed_R1_html = fastqc_trimed_R1.summary_html
        File fastqc_trimed_R2_html = fastqc_trimed_R2.summary_html
        File preprocessing_qc_json = preprocessing_qc.report_json

        # virus genotyping: irma, blasn
        File? irma_qc_json = irma_qc_flu.irma_qc_json
        File? irma_report_json = report_flu.report

        # virus genotyping: kraken, bracken, krona
        File kraken_txt = kraken2.report_txt
        File bracken_txt = bracken.report_txt
        File krona_kraken_html = krona_kraken.report_html
        File kraken_virus_txt = kraken2_vir.report_txt
        File bracken_virus_txt = bracken_vir.report_txt

        # finding aaSubstitutions: nextclade of irma's fasta
        File? HA_nextclade_tsv = nextclade_flu_ha.nextclade_tsv
        File? NA_nextclade_tsv = nextclade_flu_na.nextclade_tsv
        File? HA_nextclade_report_json = nextclade_parse_flu_ha.nextclade_json
        File? NA_nextclade_report_json = nextclade_parse_flu_ha.nextclade_json
        String? HA_nextclade_coverage_percentage = nextclade_flu_ha.coverage_percentage
        String? NA_nextclade_coverage_percentage = nextclade_flu_na.coverage_percentage

        # finding mutations (aaSubstitutions included): snpeff of minimap's fasta
        String minimap_count = minimap_flu_ha.total_count # host_filtered x 2
        String HA_minimap_matched_count = minimap_flu_ha.matched_count
        String NA_minimap_matched_count = minimap_flu_na.matched_count
        String? HA_ref_length_bp = samtools_flu_ha.ref_length_bp
        String? NA_ref_length_bp = samtools_flu_na.ref_length_bp
        String? HA_minimap_coverage_bp = samtools_flu_ha.coverage_bp
        String? NA_minimap_coverage_bp = samtools_flu_na.coverage_bp
        String? HA_minimap_coverage_percentage = samtools_flu_ha.coverage_percentage
        String? NA_minimap_coverage_percentage = samtools_flu_na.coverage_percentage
        File? HA_minimap_qc_json = minimap_qc_ha.qc_json
        File? NA_minimap_qc_json = minimap_qc_na.qc_json
        File? HA_vsf = samtools_flu_ha.file_vsf
        File? NA_vsf = samtools_flu_na.file_vsf
        File? HA_filtered_vsf = samtools_flu_ha.filtered_vsf # min_cov=5
        File? NA_filtered_vsf = samtools_flu_na.filtered_vsf # min_cov=5
        File? HA_consensus_fasta = samtools_flu_ha.consensus_fasta
        File? NA_consensus_fasta = samtools_flu_na.consensus_fasta
        File? HA_snpeff_vcf = snpeff_ha.snpeff_vcf
        File? NA_snpeff_vcf = snpeff_na.snpeff_vcf
        File? HA_snpeff_csv = snpeff_ha.snpeff_csv
        File? NA_snpeff_csv = snpeff_na.snpeff_csv
        File? HA_snpeff_report_json = parse_snpeff_ha.snpeff_json
        File? NA_snpeff_report_json = parse_snpeff_na.snpeff_json

        # final report: irma, nextclade, snpeff
        File? summary_report_json = summary_report_flu.report_json

    }
}
