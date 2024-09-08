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
}
