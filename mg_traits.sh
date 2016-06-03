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


###########################################################################################################
# 1 - Check database connection
###########################################################################################################


DB_RESULT=$( echo "UPDATE mg_traits.mg_traits_jobs SET time_started = now(), job_id = ${JOB_ID}, cluster_node = '${HOSTNAME}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = ${MG_ID};" \
| psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" )


if [[ "$?" -ne "0" ]]; then
  email_comm "Cannot connect to database. Output:${DB_RESULT}"
  exit 2
fi

if [[ "${DB_RESULT}" != "UPDATE 1" ]]; then
  email_comm "sample name ${SAMPLE_LABEL} is not in database Result:${DB_RESULT}"  
  exit 2
fi


###########################################################################################################
# 2 - Create job directory
###########################################################################################################

echo "This job tmp dir: ${THIS_JOB_TMP_DIR}"; 

rm -r ${THIS_JOB_TMP_DIR}  # CHANGE THIS FOR REAL DATA!!!!!!!!!!
mkdir "${THIS_JOB_TMP_DIR}" && cd "${THIS_JOB_TMP_DIR}"
mkdir "${THIS_JOB_TMP_DIR_DATA}" && mkdir "${SINA_LOG_DIR}"


echo "Logs, data and temp files will be written to:$(pwd)"
if [[ "$(pwd)" != "${THIS_JOB_TMP_DIR}" ]]; then 
 email_comm  "Could not access job temp dir ${THIS_JOB_TMP_DIR}"
 db_error_comm "Could not access job temp dir ${THIS_JOB_TMP_DIR}"
 
 exit 2; 
fi


###########################################################################################################
# 3 - Download file from MG URL
###########################################################################################################

# validate MG URL 
echo "${MG_URL}"
REGEX='(https?|ftp|file)://[-A-Za-z0-9\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\+&@#/%=~_|]'

if [[ ! ${MG_URL} =~ ${REGEX} ]]; then
  
  DB_COM=$( db_error_com "Not a valid URL" ) 
  email_comm "Invalid URL ${MG_URL} output db: $DB_COM"
  
  exit 1
fi 



# check if it already exist on our DB
# URLDB=$(psql -t -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" -c \
# "SELECT count(*) FROM mg_traits.mg_traits_jobs where mg_url = '${MG_URL}' AND sample_label NOT ILIKE '%test% AND return_code = 0'")
#     
# if [[ "${URLDB}" -gt 1 ]]; then 
#   email_comm "The URL "${MG_URL}" has been already succesfully crunched. If the file is different please change the file name"
#   db_error_comm "The URL ${MG_URL} has been already succesfully crunched. If the file is different please change the file name."
#   exit 1
# fi


# download MG_URL
printf "Downloading ${MG_URL} to ${RAW_DOWNLOAD}..."
curl -s "${MG_URL}" > "${RAW_DOWNLOAD}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "Could not retrieve ${MG_URL}"
  db_error_comm  "Could not retrieve ${MG_URL}"
  exit 1
fi

# compress data
gunzip -qc "${RAW_DOWNLOAD}" > "${RAW_FASTA}"
if [[ "$?" -ne "0" ]]; then 
  echo "File was uncompressed"
  rm "${RAW_FASTA}"; mv "${RAW_DOWNLOAD}" "${RAW_FASTA}" # NEEDS REVIEW: TRAP AND DB COMMUNICATION?
fi
  
###########################################################################################################
# 4 -  Validate file
###########################################################################################################

printf "Validating file..."
${fasta_file_check} "${RAW_FASTA}" "${FASTA_BAD}"
FASTA_ERROR_CODE="$?"

if [[ "${FASTA_ERROR_CODE}" -ne "0" ]]; then
  FASTA_BAD_HEADER=$(grep '>' "${FASTA_BAD}" | tr -d '>'); 
  
  email_comm "${MG_URL} is not a valid FASTA file. FASTA validation failed at sequence ${FASTA_BAD_HEADER}, error: ${FASTA_ERROR_CODE}. See ${FASTA_BAD}. ${fasta_file_check} ${RAW_FASTA} ${FASTA_BAD}"
  db_error_comm  "${MG_URL} is not a valid FASTA file. Sequence validation failed. Error: ${FASTA_ERROR_CODE}. See ${FASTA_BAD}"
  
  exit 1
fi


###########################################################################################################
# 5 - Check for utilities and directories
###########################################################################################################

ERROR_MESSAGE=$(\
check_required_writable_directories "${temp_dir}" "${RUNNING_JOBS_DIR}" "${FAILED_JOBS_DIR}" "${job_out_dir}";
check_required_readable_directories "${mg_traits_dir}"; 
check_required_programs "${cd_hit_dup}" "${cd_hit_est}" "${cd_hit_mms}" "${uproc}" "${r_interpreter}";\
)


if [[ -n "${ERROR_MESSAGE}" ]]; then
  email_comm "Not found ${ERROR_MESSAGE}"
  db_error_comm "Not found ${ERROR_MESSAGE}"
  exit 2; 
fi

###########################################################################################################
# 6 - Download data files from SVN
###########################################################################################################

# pfam downlaod
echo "${PFAM_ACCESSIONS_URL}"
curl -s "${PFAM_ACCESSIONS_URL}" > "${PFAM_ACCESSIONS}"

if [[ "$?" -ne "0" ]]; then
  email_comm "Could not retrieve ${PFAM_ACCESSIONS_URL} $http_proxy ${PFAM_ACCESSIONS}"
  db_error_comm "Could not retrieve ${PFAM_ACCESSIONS_URL}"
  exit 1; 
fi


# tf_file download
echo "${TFFILE_URL}"
curl -s "${TFFILE_URL}" > "${TFFILE}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "Could not retrieve ${TFFILE_URL}"
  db_error_comm "Could not retrieve ${TFFILE_URL}"
  exit 1;
fi
  

# slv download
echo "${SLV_TAX_URL}"
curl -s "${SLV_TAX_URL}" > "${SLV_FILE}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "Could not retrieve ${SLV_TAX_URL}"
  db_error_comm "Could not retrieve ${SLV_TAX_URL}"
  exit 1; 
fi

###########################################################################################################
# 7 - Check for duplicates
###########################################################################################################

printf "Removing duplicated sequences..."
"${cd_hit_dup}" -i "${RAW_FASTA}" -o "${UNIQUE}" > "${UNIQUE_LOG}"
CD_HIT_ERROR_CODE="$?"

if [[ "${CD_HIT_ERROR_CODE}" -ne "0" ]]; then 
  email_comm "${cd_hit_dup} -i ${RAW_FASTA} -o /dev/null > ${UNIQUE_LOG}
exited with RC ${CD_HIT_ERROR_CODE} in job ${JOB_ID}\nFiles are at: ${FAILED_JOBS_DIR}/job-${JOB_ID}"
  db_error_comm "${MG_URL} could not be processed by cd-hit-dup"
  exit 2;
fi

NUM_READS=$(grep 'Total number of sequences:'  "${UNIQUE_LOG}" | awk '{print $(NF)}')
NUM_UNIQUE=$(grep 'Number of clusters found:'  "${UNIQUE_LOG}" | awk '{print $(NF)}')

echo "Number of sequences: ${NUM_READS}"
echo "Number of unique sequences: ${NUM_UNIQUE}"

if [[ "$NUM_READS" -ne "$NUM_UNIQUE" ]]; then 
  email_comm "${MG_URL} contains duplicates. Please provide a pre-processed metagenome."
  db_error_comm "${MG_URL} contains duplicates"
  exit 1; 
fi

###########################################################################################################
# 8 - Calculate sequence statistics
###########################################################################################################

printf "Calculating sequence statistics..."

# infoseq
infoseq "${RAW_FASTA}" -only -pgc -length -noheading -auto > "${INFOSEQ_TMPFILE}"

if [[ "$?" -ne "0" ]]; then  
  email_comm "infoseq ${RAW_FASTA} -only -pgc -length -noheading -auto > ${INFOSEQ_TMPFILE}
exited with RC $? in job ${JOB_ID} Files are at: ${FAILED_JOBS_DIR}/job-${JOB_ID}"
  db_error_comm "Cannot calculate sequence statistics. Please contact adminitrator."
  exit 2; 
fi


# seq stats
"${r_interpreter}" --vanilla --slave seq_stats.R "${INFOSEQ_MGSTATS}" "${INFOSEQ_TMPFILE}"
SEQ_STATS_ERROR_CODE="$?"

if [[ "${SEQ_STATS_ERROR_CODE}" -ne "0" ]]; then 
  email_comm "${r_interpreter} --vanilla --slave
exited with RC ${SEQ_STATS_ERROR_CODE} in job ${JOB_ID}. Infoseq script. Files are at: ${FAILED_JOBS_DIR}/job-${JOB_ID}"
  db_error_comm "Cannot process sequence statistics. Please contact adminitrator."
  exit 2; 
fi

NUM_BASES=$(cut -f1 "${INFOSEQ_MGSTATS}" -d ' '); 
GC=$(cut -f2 "${INFOSEQ_MGSTATS}" -d ' '); 
VARGC=$(cut -f3 "${INFOSEQ_MGSTATS}" -d ' ')
printf "Number of bases: %d\nGC content: %f\nGC variance: %f\n" "${NUM_BASES}" "${GC}" "${VARGC}"



###########################################################################################################
# 1 - run fgs
###########################################################################################################

mkdir split_qc && cd split_qc
#Split original
printf "Splitting file ("${NSEQ}" seqs file)..."
awk -vn="${NSEQ}" 'BEGIN {n_seq=0;partid=1;} /^>/ {if(n_seq%n==0){file=sprintf("05-part-%d.fasta",partid);partid++;} print >> file; n_seq++; next;} { print >> file; }' < ../"${RAW_FASTA}"

"${fgs_runner}" "${RAW_FASTA}" "${NSLOTS}" "${NSEQ}"
ERROR_FGS=$?

if [[ "${ERROR_FGS}" -ne "0" ]]; then
  email_comm  "${frag_gene_scan} -genome=${IN_FASTA_FILE} -out=${IN_FASTA_FILE}.genes10 -complete=0 -train=sanger_10
exited with RC ${ERROR_FGS} in job ${JOB_ID}"
  db_error_comm "FragGeneScan failed. Please contact adminitrator."
  exit 2
fi

cd ../


###########################################################################################################
# 2 - run sortmerna
###########################################################################################################


"${sortmerna_runner}" "${NAM}" "${NSLOTS}" "${RES}"
ERROR_SORTMERNA=$?

if [[ "${ERROR_SORTMERNA}" -ne "0" ]]; then
  email_comm "${sortmerna} --reads ${SE} -a ${NSLOTS} --ref ${DB}/rRNA_databases/silva-bac-16s-id90.fasta ...
exited with RC ${ERROR_SORTMERNA} in job ${JOB_ID}."
  db_error_comm "sormerna failed. Please contact adminitrator"
  exit 2
fi


###########################################################################################################
# 3 - run SINA
###########################################################################################################

mkdir split_smr && cd split_smr

"${sina_runner}" "${NAM}" "${NSLOTS}" "${nSEQ}" "${RES}"
ERROR_SINA=$?

if [[ "${ERROR_SINA}" -ne "0" ]]; then
  email_comm "${sina} -i ${RES}/split_smr/spout_\${SGE_TASK_ID} -o ${RES}/split_smr/\${SGE_TASK_ID}.16S.align.fasta ...
exited with RC ${ERROR_SINA} in job ${JOB_ID}."
  db_error_comm "sina failed. Please contact adminitrator"
  exit 2
fi
cd ../


END_TIME=`date +%s.%N`