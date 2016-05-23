#!/bin/bash
set -x
set -o pipefail

START_TIME=`date +%s.%N`

echo "Environment variables:"

source ~/.bashrc
source config.bash

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
# Parse parameters
###########################################################################################################

# urldecode input
url_test="https://owncloud.mpi-bremen.de/index.php/s/RDB4Jo0PAayg3qx/download?path=%2F2014%2Fdatasets%2Fworkable%2Fmetagenomes%2Fnon-merged&files=OSD1_R1_shotgun_workable.fastq.gz"

string=$(echo $1 | sed -e 's/&/|/g' -e 's/\+/ /g' -e 's/%25/%/g' -e 's/%20/ /g' -e 's/%09/ /g' -e 's/%21/!/g' -e 's/%22/"/g' -e 's/%23/#/g' -e 's/%24/\$/g' -e 's/%26/\&/g' -e 's/%27/'\''/g' -e 's/%28/(/g' -e 's/%29/)/g' -e 's/%2a/\*/g' -e 's/%2b/+/g' -e 's/%2c/,/g' -e 's/%2d/-/g' -e 's/%2e/\./g' -e 's/%2f/\//g' -e 's/%3a/:/g' -e 's/%3b/;/g' -e 's/%3d/=/g' -e 's/%3e//g' -e 's/%3f/?/g' -e 's/%40/@/g' -e 's/%5b/\[/g' -e 's/%5c/\\/g' -e 's/%5d/\]/g' -e 's/%5e/\^/g' -e 's/%5f/_/g' -e 's/%60/`/g' -e 's/%7b/{/g' -e 's/%7c/|/g' -e 's/%7d/}/g' -e 's/%7e/~/g' -e 's/%09/      /g')

string=$(echo $url_test | sed -e 's/&/|/g' -e 's/\+/ /g' -e 's/%25/%/g' -e 's/%20/ /g' -e 's/%09/ /g' -e 's/%21/!/g' -e 's/%22/"/g' -e 's/%23/#/g' -e 's/%24/\$/g' -e 's/%26/\&/g' -e 's/%27/'\''/g' -e 's/%28/(/g' -e 's/%29/)/g' -e 's/%2a/\*/g' -e 's/%2b/+/g' -e 's/%2c/,/g' -e 's/%2d/-/g' -e 's/%2e/\./g' -e 's/%2f/\//g' -e 's/%3a/:/g' -e 's/%3b/;/g' -e 's/%3d/=/g' -e 's/%3e//g' -e 's/%3f/?/g' -e 's/%40/@/g' -e 's/%5b/\[/g' -e 's/%5c/\\/g' -e 's/%5d/\]/g' -e 's/%5e/\^/g' -e 's/%5f/_/g' -e 's/%60/`/g' -e 's/%7b/{/g' -e 's/%7c/|/g' -e 's/%7d/}/g' -e 's/%7e/~/g' -e 's/%09/      /g')

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

echo -e "SAMPLE_LABEL=${SAMPLE_LABEL}\nMG_URL=${MG_URL}\nCUSTOMER=${CUSTOMER}\nSAMPLE_ENVIRONMENT=${SAMPLE_ENVIRONMENT}\nSUBMIT_TIME=${SUBMIT_TIME}\nMAKE_PUBLIC=${MAKE_PUBLIC}\nMAKE_PUBLIC=${MAKE_PUBLIC}\nKEEP_DATA=${KEEP_DATA}\nMG_ID=${MG_ID}" > tmp.vars 

###########################################################################################################
# Download file from MG URL
###########################################################################################################

# validate URL and check if it already exist on our DB
echo "${MG_URL}"
REGX='(https?|ftp|file)://[-A-Za-z0-9\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\+&@#/%=~_|]'
[[ ! "${MG_URL}" =~ $REGX && "${DB}"==1 ]] && echo "ERROR_MG_URL=1" >> tmp.vars && "${BIN}"/db_commands.bash && cd .. &&  mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}" && exit 1
[[ "${DB}"==1 ]] && URLDB=$(psql -t -U $target_db_user -h $target_db_host -p $target_db_port -d $target_db_name -c "SELECT count(*) FROM mg_traits.mg_traits_jobs where mg_url = '${MG_URL}' AND sample_label NOT ILIKE '%test% AND return_code = 0'")
[[ "${URLDB}" -gt 1 && "${DB}"==1 ]] && echo "ERROR_URLDB=1" >> tmp.vars && "${BIN}"/db_commands.bash && cd .. && mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"


printf "Downloading "${MG_URL} to "${RAW_DOWNLOAD}..."
curl -s "${MG_URL}" > "${RAW_DOWNLOAD}"
[[ "$?" -ne "0"  && "${DB}"==1 ]] && echo "ERROR_MG_DOWNLOAD=1" >> tmp.vars && "${BIN}"/db_commands.bash && cd .. &&  mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}" && exit 1

gunzip -qc "${RAW_DOWNLOAD}" > "${RAW_FASTA}"
[[ "$?" -ne "0" ]] && echo "File was uncompressed" && rm "${RAW_FASTA}" && mv "${RAW_DOWNLOAD}" "${RAW_FASTA}" # NEEDS REVIEW: TRAP AND DB COMMUNICATION?

###########################################################################################################
# Validate file
###########################################################################################################

printf "Validating file..."
"${BIN}"/fasta_file_check.pl "${RAW_FASTA}" "${BAD_FASTA}"
FASTA_ERROR_CODE="$?"

if [[ "$FASTA_ERROR_CODE" -ne "0" && "${DB}"==1 ]]; then

  FASTA_BAD_HEADER=$(grep '>' ${BAD_FASTA} | tr -d '>'); FASTA_BAD_SEQ=$( cat ${BAD_FASTA} );
  echo -e "FASTA_ERROR=${FASTA_ERROR_CODE}\nFASTA_BAD_HEADER=${FASTA_BAD_HEADER}\nBAD_FASTA=${BAD_FASTA}\nFASTA_BAD_SEQ=${FASTA_BAD_SEQ}" >> tmp.vars
  "${BIN}"/db_commands.bash; cd ..; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 1

fi

###########################################################################################################
# Create job directory
###########################################################################################################

echo "This job tmp dir: ${THIS_JOB_TMP_DIR}"; 

mkdir "${THIS_JOB_TMP_DIR}" && cd "${THIS_JOB_TMP_DIR}"
mkdir "${THIS_JOB_TMP_DIR_DATA}" && mkdir "${SINA_LOG_DIR}"

echo "Logs, data and temp files will be written to:$(pwd)"
if [[ "$(pwd)" != "$THIS_JOB_TMP_DIR" ]]; then exit 2; fi

###########################################################################################################
# Check for utilities and directories
###########################################################################################################

check_required_writeble_directories "${temp_dir}" "${RUNNING_JOBS_DIR}" "${FAILED_JOBS_DIR}" "${job_out_dir}"
echo -e "ERROR_WDIRECTORIES=${ERROR_WDIRECTORIES}" >> tmp.vars
if [[ -n "${ERROR_WDIRECTORIES}" ]]; then echo "non writable directories ${ERROR_WDIRECTORIES}"; exit 2; fi

check_required_readable_directories "${mg_traits_dir}"
echo -e "ERROR_RDIRECTORIES=${ERROR_RDIRECTORIES}" >> tmp.vars
if [[ -n "${ERROR_RDIRECTORIES}" ]]; then echo "non writable directories ${ERROR_RDIRECTORIES}"; exit 2; fi

check_required_programs "${cd_hit_dup}"	"${cd_hit_est}" "${cd_hit_mms}" "${uproc}" "${r_interpreter}"
echo -e "ERROR_UTILITIES=${ERROR_UTILITIES}" >> tmp.vars
if [[ -n "${ERROR_UTILITIES}" ]]; then echo "not found ${ERROR_UTILITIES}"; exit 2; fi

###########################################################################################################
# Download data files from SVN
###########################################################################################################

echo "${PFAM_ACCESSIONS_URL}"
curl -s "${PFAM_ACCESSIONS_URL}" > "${PFAM_ACCESSIONS}"
if [[ "$?" -ne "0" ]]; then echo -e "ERROR_PFAM=1" >> tmp.vars; cd ..; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 1; fi


echo "${TFFILE_URL}"
curl -s "${TFFILE_URL}" > "${TFFILE}"
if [[ "$?" -ne "0" ]]; then echo -e "ERROR_TF=1" >> tmp.vars; cd ..;  mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 1; fi
  

echo "${SLV_TAX_URL}"
curl -s "${SLV_TAX_URL}" > "${SLV_FILE}"
[[ "$?" -ne "0" ]] then echo "ERROR_SLV=1" >> tmp.vars; cd ..;  mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 1; fi

###########################################################################################################
# Check for duplicates
###########################################################################################################

printf "Removing duplicated sequences..."
"${cd_hit_dup} -i "${RAW_FASTA}" -o "${UNIQUE} > "${UNIQUE_LOG}"
CD_HIT_ERROR_CODE="$?"

if [[ "$?" -ne "0" ]]; then echo -e "ERROR_CD_HIT=1\nCD_HIT_ERROR_CODE=${CD_HIT_ERROR_CODE}" >> tmp.vars; cd ..; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 2; fi

NUM_READS=$(grep 'Total number of sequences:'  "${UNIQUE_LOG}" | awk '{print $(NF)}')
NUM_UNIQUE=$(grep 'Number of clusters found:'  "${UNIQUE_LOG}" | awk '{print $(NF)}')

echo "Number of sequences: ${NUM_READS}"
echo "Number of unique sequences: ${NUM_UNIQUE}"
if [[ "$NUM_READS" -ne "$NUM_UNIQUE" ]]; then echo "ERROR_NUM_UNIQUE=1" >> tmp.vars; cd ..; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 1; fi

###########################################################################################################
# Calculate sequence statistics
###########################################################################################################

printf "Calculating sequence statistics..."
infoseq "${RAW_FASTA}" -only -pgc -length -noheading -auto > "${INFOSEQ_TMPFILE}"
if [[ "$?" -ne "0" ]]; then  echo "ERROR_INFOSEQ=1" >> tmp.vars; cd ..; mv mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 2; fi

"${r_interpreter}" --vanilla --slave seq_stats.R "${INFOSEQ_MGSTATS}" "${INFOSEQ_TMPFILE}"
SEQ_STATS_ERROR_CODE="$?"
if [[ "$?" -ne "0" ]]; then echo -e "ERROR_SEQ_STATS=1\nSEQ_STATS_ERROR_CODE=${SEQ_STATS_ERROR_CODE}" cd ..; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"; exit 2; fi

NUM_BASES=$(cut -f1 "${INFOSEQ_MGSTATS}" -d ' '); GC=$(cut -f2 "${INFOSEQ_MGSTATS}" -d ' '); VARGC=$(cut -f3 "${INFOSEQ_MGSTATS}" -d ' ')
printf "Number of bases: %d\nGC content: %f\nGC variance: %f\n" "${NUM_BASES}" "${GC}" "${VARGC}"

###########################################################################################################
# Parition the data
###########################################################################################################

#Split original
printf "Splitting file ("${NSEQ}" seqs file)..."
awk -vn="${NSEQ}" 'BEGIN {n_seq=0;partid=1;} /^>/ {if(n_seq%n==0){file=sprintf("05-part-%d.fasta",partid);partid++;} print >> file; n_seq++; next;} { print >> file; }' < "${RAW_FASTA}"

SUBJOBS=$(ls -1 05-part*.fasta | wc -l)
echo "Split into ${SUBJOBS} sub jobs..."



