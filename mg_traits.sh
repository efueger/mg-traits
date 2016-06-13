#!/bin/bash
set -x
set -o pipefail

START_TIME=`date +%s.%N`

echo "Environment variables:"

source ~/.bashrc
source /bioinf/home/epereira/workspace/mg-traits/resources/config.bash
source /bioinf/home/epereira/workspace/mg-traits/resources/config.proxy

echo -e "\tJob ID: ${JOB_ID}"
echo -e "\tTarget database: ${target_db_user}@${target_db_host}:${target_db_port}/${target_db_name}"
echo -e "\tCD-HIT-DUP: ${cd_hit_dup}"
echo -e "\tCD-HIT-EST: ${cd_hit_est}"
echo -e "\tCD-HIT-MMS: ${cd_hit_mms}"
echo -e "\tFragGeneScan: ${frag_gene_scan}"
echo -e "\tUPro: ${upro}"
echo -e "\tR: ${r_interpreter}"
echo -e "\tTemp dir: ${temp_dir}"
echo -e "\tMG traits dir: ${mg_traits_dir}"
echo -e "\tJob out dir: ${job_out_dir}"
echo -e "\tMT admin mail: ${mt_admin_mail}"


###########################################################################################################
# 0 - Parse parameters
###########################################################################################################

# urldecode input
string=$(echo $1 | sed -e 's/&/|/g' -e 's/\+/ /g' -e 's/%25/%/g' -e 's/%20/ /g' -e 's/%09/ /g' -e 's/%21/!/g' -e 's/%22/"/g' -e 's/%23/#/g' -e 's/%24/\$/g' -e 's/%26/\&/g' -e 's/%27/'\''/g' -e 's/%28/(/g' -e 's/%29/)/g' -e 's/%2a/\*/g' -e 's/%2b/+/g' -e 's/%2c/,/g' -e 's/%2d/-/g' -e 's/%2e/\./g' -e 's/%2f/\//g' -e 's/%3a/:/g' -e 's/%3b/;/g' -e 's/%3d/=/g' -e 's/%3e//g' -e 's/%3f/?/g' -e 's/%40/@/g' -e 's/%5b/\[/g' -e 's/%5c/\\/g' -e 's/%5d/\]/g' -e 's/%5e/\^/g' -e 's/%5f/_/g' -e 's/%60/`/g' -e 's/%7b/{/g' -e 's/%7c/|/g' -e 's/%7d/}/g' -e 's/%7e/~/g' -e 's/%09/      /g')

# set delimiter
IFS="|"

# parse input
echo "Input parameters:"
for pair in $string; do
key=${pair%%=*}
value=${pair#*=}

printf "\t${key}=${value}\n";

if [ "${key}" = "sample_label" ]; then
	SAMPLE_LABEL="${value}";
fi

if [ "${key}" = "mg_url" ]; then
	MG_URL="${value}";
fi

if [ "${key}" = "customer" ]; then
	CUSTOMER="${value}";
fi

if [ "${key}" = "sample_environment" ]; then
	SAMPLE_ENVIRONMENT="${value}";
fi

if [ "${key}" = "time_submitted" ]; then
	SUBMIT_TIME="${value}";
fi

if [ "${key}" = "make_public" ]; then
	MAKE_PUBLIC="${value}";
fi

if [ "${key}" = "keep_data" ]; then
	KEEP_DATA="${value}";
fi

if [ "${key}" = "id" ]; then
	MG_ID="${value}";
fi

done



#####################################################################################################
# Define database communication function
#####################################################################################################

function db_error_comm() {
  echo "UPDATE mg_traits.mg_traits_jobs SET time_finished = now(), return_code = 1, error_message = '${1}' \
  WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
}


DB_RESULT=$( echo "UPDATE mg_traits.mg_traits_jobs SET time_started = now(), job_id = ${JOB_ID}, cluster_node = '${HOSTNAME}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = ${MG_ID};" \
| psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" )


cd "${THIS_JOB_TMP_DIR}"



###########################################################################################################
# 4 - Finish Jobs
###########################################################################################################


echo frag_gene_scan=$frag_gene_scan >> 00-environment
echo sina=$sina >> 00-environment
echo sina_arb_pt_server=$sina_arb_pt_server >> 00-environment
echo sina_version=$sina_version >> 00-environment
echo sina_seed=$sina_seed >> 00-environment
echo sina_seed_version=$sina_seed_version >> 00-environment
echo sina_ref=$sina_ref >> 00-environment
echo sina_ref_version=$sina_ref_version >> 00-environment
echo SINA_SOCKET=$SINA_SOCKET >> 00-environment

echo uproc=$uproc >> 00-environment
echo r_interpreter=$r_interpreter >> 00-environment
echo target_db_host=$target_db_host >> 00-environment
echo target_db_port=$target_db_port >> 00-environment
echo target_db_user=$target_db_user >> 00-environment
echo target_db_name=$target_db_name >> 00-environment
echo mt_admin_mail=$mt_admin_mail >> 00-environment
echo THIS_JOB_TMP_DIR=$THIS_JOB_TMP_DIR >> 00-environment
echo TFFILE=$TFFILE >> 00-environment
echo PFAM_ACCESSIONS=$PFAM_ACCESSIONS >> 00-environment
echo SAMPLE_LABEL=$SAMPLE_LABEL >> 00-environment
echo RAW_FASTA=$RAW_FASTA >> 00-environment
echo GC=$GC >> 00-environment
echo VARGC=$VARGC >> 00-environment
echo NUM_BASES=$NUM_BASES >> 00-environment
echo NUM_READS=$NUM_READS >> 00-environment
echo THIS_JOB_ID=$JOB_ID >> 00-environment
echo temp_dir=$temp_dir >> 00-environment
echo FAILED_JOBS_DIR=$FAILED_JOBS_DIR >> 00-environment
echo RUNNING_JOBS_DIR=$RUNNING_JOBS_DIR >> 00-environment
echo FINISHED_JOBS_DIR=$FINISHED_JOBS_DIR >> 00-environment
echo SUBJOBS=$SUBJOBS >> 00-environment
echo uproc_pfam=$uproc_pfam >> 00-environment
echo uproc_model=$uproc_model >> 00-environment
echo MG_ID=$MG_ID >> 00-environment
echo SINA_LOG_DIR=$SINA_LOG_DIR >> 00-environment
echo SLV_FILE=$SLV_FILE >> 00-environment
echo ARBHOME=$ARBHOME >> 00-environment
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH >> 00-environment
echo NSLOTS=$NSLOTS >> 00-environment



qsub -pe threaded 8 -l h=\!mg32 -N $FINISHJOBID -o $THIS_JOB_TMP_DIR -e $THIS_JOB_TMP_DIR -l ga -j y -terse -P megx.p -R y -m sa -M $mt_admin_mail -hold_jid $FGS_JOBARRAYID,$SINA_JOBARRAYID  /bioinf/home/epereira/workspace/mg-traits/resources/finish_runner.dev.sh $THIS_JOB_TMP_DIR



END_TIME=`date +%s.%N`
