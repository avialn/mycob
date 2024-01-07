#!/bin/bash

sudo shutdown -h +360
export USER_ID=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/user-id)
export TASK_ID=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/task-id)
export VM_NAME=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/name)
export ATTEMPT=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/attempt)
export FILES=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/files)
export PANEL=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/panel)
export BUILD=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/build)
export TYPE_HS=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/type_hs)
export TYPE_VC=$(curl -H Metadata-Flavor:Google 169.254.169.254/computeMetadata/v1/instance/attributes/type_vc)
VER=${BUILD:2:2}
export AWS_ACCESS_KEY_ID=YCAJETbuUo6peFsb89GqhV-CH
export AWS_SECRET_ACCESS_KEY=YCMzfRDYoDO9hAzTVZgbSblGjHdftk-WxmZSYsqc
export AWS_DEFAULT_REGION=ru-central1
touch /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
watch -n 10 touch /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log &>/dev/null &
echo "UserID: ${USER_ID}" > /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "TaskID: ${TASK_ID}" >> /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "Attempt: ${ATTEMPT}" >> /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "Panel: ${PANEL}" >> /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "Build: ${BUILD}" >> /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "Type HS: ${TYPE_HS}" >> /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "Type VC: ${TYPE_VC}" >> /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
mkdir -p /home/cromwell/fastq
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/${PANEL}/run.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/irma.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/kraken2.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/nextclade.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/preprocessing.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/amr_tasks.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/dada2.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/yandex_utilities.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/summary_report.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/common_tasks/sarscov2.wdl /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/yandex_inputs.json /home/cromwell/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/yandex_options.json /home/cromwell/options.json |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-wdl/yandex.conf /home/cromwell/ |& tee -a /home/cromwell/aws.log
zip --junk-paths /home/cromwell/imports.zip /home/cromwell/irma.wdl /home/cromwell/kraken2.wdl /home/cromwell/nextclade.wdl /home/cromwell/preprocessing.wdl /home/cromwell/yandex_utilities.wdl /home/cromwell/amr_tasks.wdl /home/cromwell/dada2.wdl /home/cromwell/summary_report.wdl /home/cromwell/sarscov2.wdl


mkdir -p /home/cromwell/mycob-ref
mkdir -p /home/cromwell/mycob-ref/blast
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/left_primers.fasta /home/cromwell/mycob-ref/rsv_full/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/right_primers.fasta /home/cromwell/mycob-ref/rsv_full/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/antigenic_frame.txt /home/cromwell/mycob-ref/rsv_full/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/human_bowtie2_index/GRCh38_ERCC.bowtie2.tar /home/cromwell/mycob-ref/rsv_full/human_bowtie2_index/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/human_bowtie2_index/GRCh38_transcriptome.bowtie2.tar /home/cromwell/mycob-ref/rsv_full/human_bowtie2_index/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/kraken2/k2_standard_08gb_20231009.tar.gz /home/cromwell/mycob-ref/rsv_full/kraken2/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/kraken2/k2_viral_20231009.tar.gz /home/cromwell/mycob-ref/rsv_full/kraken2/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/kraken2/16S_Greengenes13.5_20200326.tar.gz /home/cromwell/mycob-ref/rsv_full/kraken2/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/kraken2/k2_pluspf_08gb_20231009.tar.gz /home/cromwell/mycob-ref/rsv_full/kraken2/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/kraken2/ncbi_16s_18s_28s_ITS_kraken2.tar.gz /home/cromwell/mycob-ref/rsv_full/kraken2/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/16S/dada2/silva_nr99_v138.1_wSpecies_train_set.fa.gz /home/cromwell/mycob-ref/16S/dada2/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/MEASLES/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/MEASLES/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/ADENO/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/ADENO/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/PNEUMO/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/PNEUMO/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/BOCA/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/BOCA/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/CORONA/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/CORONA/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/RESPIRO/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/RESPIRO/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/RHINO/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/RHINO/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/RUBULA/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/RUBULA/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/CoV/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/CoV/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/METAPMEUMO/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/METAPMEUMO/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/consensus.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/HA_wis_67_2022.fast /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/NA_wis_67_2022.fasta /home/cromwell/mycob-ref/rsv_full/IRMA_RES/modules/FLU/reference/ |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp s3://mycob-ref/amr/card.json /home/cromwell/mycob-ref/amr/CARD/ |& tee -a /home/cromwell/aws.log

echo $FILES | jq --raw-output '.[] | .[]' | xargs -L1 -I'{}' aws --endpoint-url=https://storage.yandexcloud.net s3 cp {} /home/cromwell/fastq/ |& tee -a /home/cromwell/aws.log
jq --arg key0 'task-id' --arg value0 ${TASK_ID} --arg key1 'attempt' --arg value1 ${ATTEMPT} '. | .[$key0]=$value0 | .[$key1]=$value1' <<<'{}' > /home/cromwell/labels.json
export LOCAL_FILES=$(echo $FILES | jq 'map(. | map(. | sub("s3:\/\/mycob-userdata"; "\/home\/cromwell")))' | jq -c 'map(. | map(. | sub(env.USER_ID; "fastq")))')
jq --arg key0 'processing.Files' --argjson value0 ${LOCAL_FILES} \
	--arg key1 'processing.TaskID' --arg value1 ${TASK_ID} \
	--arg key2 'processing.Attempt' --arg value2 ${ATTEMPT} \
	'. | .[$key0]=$value0 | .[$key1]=$value1 | .[$key2]=$value2' <<<'{}' > /home/cromwell/local_inputs.json
jq --slurp '.[0] * .[1]' /home/cromwell/local_inputs.json /home/cromwell/yandex_inputs.json > /home/cromwell/inputs.json

sudo systemctl start env-cromwell.service

sudo sysctl --system
cd /home/cromwell

while :
do
   VERSION=$(curl --silent -X GET "http://127.0.0.1:8000/engine/v1/version" -H "accept: application/json" | jq --raw-output '.cromwell')
   DATE=$(date)
   echo ${DATE}   ${VERSION} |& tee -a /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
   sleep 10
   if [ "${VERSION}" = "65" ]
   then
	break
   fi
done

SUBMIT=$(curl --silent -X POST "http://127.0.0.1:8000/api/workflows/v1" -H "accept: application/json" -H "Content-Type: multipart/form-data" -F "workflowSource=@/home/cromwell/run.wdl" -F "workflowInputs=@/home/cromwell/inputs.json;type=application/json" -F "workflowOptions=@/home/cromwell/options.json;type=application/json" -F "labels=@/home/cromwell/labels.json;type=application/json" -F "workflowDependencies=@/home/cromwell/imports.zip;type=application/x-zip-compressed")
echo "Submiting:" ${SUBMIT} |& tee -a /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
WORKFLOW_ID=$(echo ${SUBMIT} | jq --raw-output ".id")
echo "Workflow:" ${WORKFLOW_ID} |& tee -a /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo "VM name:" ${VM_NAME} |& tee -a /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
echo ${WORKFLOW_ID} > /home/cromwell/workflow.${TASK_ID}.${ATTEMPT}
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/workflow.${TASK_ID}.${ATTEMPT} s3://mycob-cromwell-outputs/workflow/workflow.${TASK_ID}.${ATTEMPT} |& tee -a /home/cromwell/aws.log
#CALLBACK=$(curl --silent -X GET "https://seq24.ru/api/yandex_cromwell/id_set/${VM_NAME}/${WORKFLOW_ID}")
#echo "Callback:" ${CALLBACK} |& tee -a /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log

while :
do
   STATUS=$(curl --silent -X GET "http://127.0.0.1:8000/api/workflows/v1/${WORKFLOW_ID}/metadata?includeKey=status" -H "accept: application/json" | jq --raw-output '.status')
   DATE=$(date)
   echo ${DATE}   ${STATUS} |& tee -a /s3/mycob-cromwell-logs/task_logs/task.${TASK_ID}.${ATTEMPT}.log
   sleep 30
   if [ "${STATUS}" = "Running" ] && [ "${WORKFLOW_TAILED}" != 1 ]
   then
      tail -f /home/cromwell/cromwell-workflow-logs/workflow.${WORKFLOW_ID}.log |& sudo tee -a /dev/ttyS2 > /dev/null &
      WORKFLOW_TAILED=1
   fi
   if [ "${STATUS}" != "" ]
   then
      echo ${STATUS} > /home/cromwell/status.${WORKFLOW_ID}
   fi
   aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/status.${WORKFLOW_ID} s3://mycob-cromwell-outputs/status/status.${WORKFLOW_ID} |& tee -a /home/cromwell/aws.log
   if [ "${STATUS}" = "Succeeded" ] || [ "${STATUS}" = "Failed" ] || [ "${STATUS}" = "Aborted" ]
   then
	echo ${STATUS} > /home/cromwell/status.${WORKFLOW_ID}
	break
   fi
done


curl --silent -o /home/cromwell/metadata.${WORKFLOW_ID} -X GET "http://127.0.0.1:8000/api/workflows/v1/${WORKFLOW_ID}/metadata?expandSubWorkflows=true" -H "accept: application/json"
curl --silent -o /home/cromwell/outputs.${WORKFLOW_ID} -X GET "http://127.0.0.1:8000/api/workflows/v1/${WORKFLOW_ID}/outputs" -H "accept: application/json"
curl --silent -o /home/cromwell/timing.${WORKFLOW_ID}.html -X GET "http://127.0.0.1:8000/api/workflows/v1/${WORKFLOW_ID}/timing" -H "accept: application/json"
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/metadata.${WORKFLOW_ID} s3://mycob-cromwell-metadata/metadata.${WORKFLOW_ID} |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/outputs.${WORKFLOW_ID} s3://mycob-cromwell-outputs/outputs/outputs.${WORKFLOW_ID} |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/timing.${WORKFLOW_ID}.html s3://mycob-cromwell-outputs/timing/timing.${WORKFLOW_ID}.html |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp $(cat "/home/cromwell/outputs.${WORKFLOW_ID}" | jq --raw-output '.outputs."human.full_tsv"') s3://mycob-temp/${WORKFLOW_ID}.tsv.gz |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp $(cat "/home/cromwell/outputs.${WORKFLOW_ID}" | jq --raw-output '.outputs."human.region_coverage"') s3://mycob-temp/${WORKFLOW_ID}_coverage.regions.bed.gz |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp $(cat "/home/cromwell/outputs.${WORKFLOW_ID}" | jq --raw-output '.outputs."human.bp_coverage"') s3://mycob-temp/${WORKFLOW_ID}_coverage.bp.bed.gz |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/mycob-cromwell-outputs/ s3://mycob-cromwell-outputs/ --recursive |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/mycob-cromwell-logs/wf_logs/ s3://mycob-cromwell-logs/wf_logs/ --recursive |& tee -a /home/cromwell/aws.log
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/mycob-cromwell-logs/call_logs/ s3://mycob-cromwell-logs/call_logs/ --recursive |& tee -a /home/cromwell/aws.log

tar cfv - \
	--hard-dereference \
	--exclude='GRCh38_ERCC.bowtie2.tar' \
	--exclude='GRCh38_transcriptome.bowtie2.tar' \
	--exclude='k2_standard_08gb_20231009.tar.gz' \
	--exclude='k2_viral_20231009.tar.gz' \
    --exclude='1000G_phase1.snps.high_confidence.hg38.vcf.gz' \
	--exclude='1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi' \
	--exclude='1000g_pon.hg38.vcf.gz' \
	--exclude='1000g_pon.hg38.vcf.gz.tbi' \
	--exclude='af-only-gnomad.hg38.vcf.gz' \
	--exclude='af-only-gnomad.hg38.vcf.gz.tbi' \
	--exclude='dbsnp_138.b38.vcf.gz' \
	--exclude='dbsnp_138.b38.vcf.gz.tbi' \
	--exclude='hapmap_3.3.hg38.vcf.gz' \
	--exclude='hapmap_3.3.hg38.vcf.gz.tbi' \
	--exclude='Homo_sapiens_assembly38.bed.gz' \
	--exclude='Homo_sapiens_assembly38.bed.gz.tbi' \
	--exclude='Homo_sapiens_assembly38.dict' \
	--exclude='Homo_sapiens_assembly38.fasta' \
	--exclude='Homo_sapiens_assembly38.fasta.64.alt' \
	--exclude='Homo_sapiens_assembly38.fasta.64.amb' \
	--exclude='Homo_sapiens_assembly38.fasta.64.ann' \
	--exclude='Homo_sapiens_assembly38.fasta.64.bwt' \
	--exclude='Homo_sapiens_assembly38.fasta.64.pac' \
	--exclude='Homo_sapiens_assembly38.fasta.64.sa' \
	--exclude='Homo_sapiens_assembly38.fasta.fai' \
	--exclude='Homo_sapiens_assembly38.known_indels.vcf.gz' \
	--exclude='Homo_sapiens_assembly38.known_indels.vcf.gz.tbi' \
	--exclude='Homo_sapiens_assembly38.primary_contigs.list' \
	--exclude='Homo_sapiens_assembly38_repeats.interval_list' \
	--exclude='Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' \
	--exclude='Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi' \
	--exclude='small_exac_common_3.hg38.vcf.gz' \
	--exclude='small_exac_common_3.hg38.vcf.gz.tbi' \
	--exclude='clinvar.vcf.gz' \
	--exclude='clinvar.vcf.gz.tbi' \
	--exclude='hg38_avsnp150.txt' \
	--exclude='hg38_avsnp150.txt.idx' \
	--exclude='hg38_dbnsfp35a.txt' \
	--exclude='hg38_dbnsfp35a.txt.idx' \
	--exclude='hg38_dbscsnv11.txt' \
	--exclude='hg38_dbscsnv11.txt.idx' \
	--exclude='hg38_gnomad211_exome.txt' \
	--exclude='hg38_gnomad211_exome.txt.idx' \
	--exclude='hg38_gnomad211_genome.txt' \
	--exclude='hg38_gnomad211_genome.txt.idx' \
	--exclude='hg38_intervar_20180118.txt' \
	--exclude='hg38_intervar_20180118.txt.idx' \
	--exclude='*.bam' \
    --exclude='*.sam' \
	--exclude='*.bam.bai' \
	--exclude='*.fastq' \
	--exclude='*.fastq.gz' \
	/home/cromwell/cromwell-executions/ | pigz -p 16 | aws --endpoint-url=https://storage.yandexcloud.net s3 cp - s3://mycob-cromwell-execution/${WORKFLOW_ID}.tar.gz
aws --endpoint-url=https://storage.yandexcloud.net s3 cp /home/cromwell/status.${WORKFLOW_ID} s3://mycob-cromwell-outputs/status/status.${WORKFLOW_ID} |& tee -a /home/cromwell/aws.log

cat /home/cromwell/aws.log | pigz -p 8 | aws --endpoint-url=https://storage.yandexcloud.net s3 cp - s3://mycob-cromwell-logs/task_logs/aws.${TASK_ID}.${ATTEMPT}.log.gz
#CROMWELL_ID=$(cat metadata.json | jq --raw-output '.id')
echo "Shuthdown VM in 60 sec..."
sleep 60
sudo shutdown -h now

