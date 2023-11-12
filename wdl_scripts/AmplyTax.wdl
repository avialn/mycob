version 1.0

#java -Dconfig.file=backend.conf -jar tools/cromwell-65.jar run wdl_scripts/AmplyTax.wdl -i inputs.json -o options.json

#kraken2_db:
#kraken2_ncbi (archaeal,bacterial and fungal 16S/18S and ITS): https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-metagenomics/ncbi_16s_18s_28s_ITS/ncbi_16s_18s_28s_ITS_kraken2.tar.gz
#kraken2_standard_8gb (archaea,bacteria,viral,plasmid,human,UniVec_Core): https://benlangmead.github.io/aws-indexes/k2
#kraken2_pluspf_8gb (standard plus protozoa & fungi): https://benlangmead.github.io/aws-indexes/k2
#kraken2_16s: https://benlangmead.github.io/aws-indexes/k2

#dada2_db:
#microb_16s:
#silva_138 (bacteria and archaea, prokaryotic 16/18s): https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
#sbdi-gtdb (16s by conserved single-copy proteins): https://scilifelab.figshare.com/ndownloader/files/36980767
#rdp (bacteria and archaeal 16S rRNA, fungal 28s, genus_level): https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
#fungi_its:
#unite_fungi (all eukaryotes): https://unite.ut.ee/repository.php
#unite_all (fungi): https://unite.ut.ee/repository.php
#fungi_18s:
#silva_132 (eukaryotic 18S, v132 & v128): https://zenodo.org/record/1447330
#pr2 (18S protists, metazoa, fungi and plants): https://github.com/pr2database/pr2database/releases/tag/v5.0.0

import "./common_tasks/preprocessing.wdl" as preprocessing
import "./common_tasks/kraken2.wdl" as kraken2
import "./common_tasks/dada2.wdl" as dada2

workflow AmplyTax {

    input {
        File fastq_1
        File fastq_2
        String sample_name
        Boolean cut_primers = false
        File? primer_left = "primers/left_primers.fasta"
        File? primer_right = "primers/right_primers.fasta"
        Boolean transcriptome_filtering = false # false -> index_genome, true -> index_transcriptome
        File index_genome = "db/bowtie2/GRCh38_ERCC.bowtie2.tar"
        File index_transcriptome = "db/bowtie2/GRCh38_transcriptome.bowtie2.tar"

        #kraken/bracken classifiers
        File kraken2_ncbi = "db/kraken2/ncbi_16s_18s_28s_ITS_kraken2.tar.gz" # default
        File kraken2_standard_8gb = "db/kraken2/k2_standard_08gb_20230605.tar.gz"
        File kraken2_pluspf_8gb = "db/kraken2/k2_pluspf_08gb_20231009.tar.gz"
        File kraken2_16s = "db/kraken2/16S_Greengenes13.5_20200326.tar.gz"

        Boolean microb_16s = true
        #microb_16s dada2 classifiers
        File silva_138 = "db/dada2/silva_nr99_v138.1_wSpecies_train_set.fa.gz" #default
        File sbdi_gtdb = "db/dada2/36980767"
        File rdp = "db/dada2/rdp_train_set_18.fa.gz"

        Boolean fungi_its = false
        #fungi_ITS dada2 classifiers
        File unite_fungi = "db/dada2/sh_general_release_dynamic_18.07.2023.fasta"   #default
        File unite_all = "db/dada2/sh_general_release_dynamic_all_18.07.2023.fasta"

        Boolean fungi_18s = false
        #fungi_18s dada2 classifiers
        File silva_132 = "db/dada2/silva_132.18s.99_rep_set.dada2.fa.gz" #default
        File pr2= "db/dada2/pr2_version_5.0.0_SSU_dada2.fasta.gz"

        Int threads = 8
        Int trim_f = 21
        Int trim_r = 26
        Int trunc_q = 2
        Int trunc_f = 250 #V3V4 is 460 nts with primers, so trunc_f+trunc_r > amplicon's length + 12
        Int trunc_r = 210 #250 and 210 should be shorter than R1 and R1 lengths
        Int minOverlap = 12 #если сумма ридов не покрывают длину ампликона, можно увеличить trunc_f+trunc_r или уменьшить minOverlap
        String kraken_level = "S" #(U)nclassified, (R)oot, (D)omain, (K)ingdom (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
    }

    call preprocessing.CheckInput as check_input {
        input:
        fastq = fastq_1,
        min_reads = 500,
        docker = "resouer/ubuntu-bc:latest"
    }

    if (check_input.proceed == "yes") {

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
        call preprocessing.Cutadapt as cutadapt_kraken {
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
        fastq_1 = if cut_primers then cutadapt_kraken.cut_fastq_1 else trimmomatic.trim_fastq_1,
        fastq_2 = if cut_primers then cutadapt_kraken.cut_fastq_2 else trimmomatic.trim_fastq_2,
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
        cut_input_reads = if defined(cutadapt_kraken.input_reads_cut) then cutadapt_kraken.input_reads_cut else "cut_primers is false",
        cut_both_surviving = if defined(cutadapt_kraken.both_surviving_cut) then cutadapt_kraken.both_surviving_cut else "no primers were cut",
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
        kraken2_classifier = kraken2_ncbi,
        threads = threads,
        docker = "staphb/kraken2:latest"
      }

      call kraken2.Bracken as bracken {
        input:
        sample_name = sample_name,
        kraken_report = kraken2.report_txt,
        kraken2_classifier = kraken2_ncbi,
        level = kraken_level,
        docker = "nanozoo/bracken:2.8--dcb3e47"
      }

      call kraken2.Krona as krona_kraken {
        input:
        sample_name = sample_name,
        report = kraken2.report_txt,
        docker = "krona:2.8.1"
      }

      if (cut_primers) {
        call preprocessing.Cutadapt as cutadapt_dada2 {
            input:
            fastq_1 = fastq_1,
            fastq_2 = fastq_2,
            sample_name = sample_name,
            primer_left = primer_left,
            primer_right = primer_right,
            docker = "cutadapt:4.4"
        }
      }

      if (microb_16s) {
        call dada2.Dada2 as dada2_microb_16s {
            input:
            fastq_1 = if cut_primers then cutadapt_dada2.cut_fastq_1 else fastq_1,
            fastq_2 = if cut_primers then cutadapt_dada2.cut_fastq_2 else fastq_2,
            sample_name = sample_name,
            trim_f = trim_f,
            trim_r = trim_r,
            trunc_q = trunc_q,
            trunc_f = trunc_f,
            trunc_r = trunc_r,
            minOverlap = minOverlap,
            dada2_classifier = silva_138,
            docker = "dada2:1.26"
        }
      }

      if (fungi_18s) {
        call dada2.Dada2 as dada2_fungi_18s {
            input:
            fastq_1 = if cut_primers then cutadapt_dada2.cut_fastq_1 else fastq_1,
            fastq_2 = if cut_primers then cutadapt_dada2.cut_fastq_2 else fastq_2,
            sample_name = sample_name,
            trim_f = trim_f,
            trim_r = trim_r,
            trunc_q = trunc_q,
            trunc_f = trunc_f,
            trunc_r = trunc_r,
            minOverlap = minOverlap,
            dada2_classifier = pr2,
            docker = "dada2:1.26"
        }
      }

      if (fungi_its) {
        call dada2.Dada2 as dada2_fungi_its  {
            input:
            fastq_1 = if cut_primers then cutadapt_dada2.cut_fastq_1 else fastq_1,
            fastq_2 = if cut_primers then cutadapt_dada2.cut_fastq_2 else fastq_2,
            sample_name = sample_name,
            trim_f = trim_f,
            trim_r = trim_r,
            trunc_q = trunc_q,
            trunc_f = trunc_f,
            trunc_r = trunc_r,
            minOverlap = minOverlap,
            dada2_classifier = unite_fungi,
            docker = "dada2:1.26"
        }
      }
    }

    output {
        File? fastqc_row_R1_html = fastqc_row_R1.summary_html
        File? fastqc_row_R2_html = fastqc_row_R2.summary_html
        File? fastqc_trimed_R1_html = fastqc_trimed_R1.summary_html
        File? fastqc_trimed_R2_html = fastqc_trimed_R2.summary_html
        File? preprocessing_qc_json = preprocessing_qc.report_json
        File? kraken_txt = kraken2.report_txt
        File? bracken_txt = bracken.report_txt
        File? krona_kraken_html = krona_kraken.report_html
        File? dada2_microb_16s_tsv = dada2_microb_16s.seqtab_nochim_tsv
        File? stat_dada2_microb_16s_tsv = dada2_microb_16s.stat_tsv
        File? dada2_fungi_18s_tsv = dada2_fungi_18s.seqtab_nochim_tsv
        File? stat_dada2_fungi_18s_tsv = dada2_fungi_18s.stat_tsv
        File? dada2_fungi_its_tsv = dada2_fungi_its.seqtab_nochim_tsv
        File? stat_dada2_fungi_its_tsv = dada2_fungi_its.stat_tsv
    }

}
