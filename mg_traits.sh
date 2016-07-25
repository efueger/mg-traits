#!/bin/bash
set -x
set -o pipefail

START_TIME=$(date +%s.%N)

echo "Environment variables:"

source ~/.bashrc
source /bioinf/projects/megx/mg-traits/resources/config_files/config.bash

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

if [[ "${key}" = "sample_label" ]]; then
	SAMPLE_LABEL="${value}";
fi

if [[ "${key}" = "mg_url" ]]; then
	MG_URL="${value}";
fi

if [[ "${key}" = "customer" ]]; then
	CUSTOMER="${value}";
fi

if [[ "${key}" = "sample_environment" ]]; then
	SAMPLE_ENVIRONMENT="${value}";
fi

if [[ "${key}" = "time_submitted" ]]; then
	SUBMIT_TIME="${value}";
fi

if [[ "${key}" = "make_public" ]]; then
	MAKE_PUBLIC="${value}";
fi

if [[ "${key}" = "keep_data" ]]; then
	KEEP_DATA="${value}";
fi

if [[ "${key}" = "id" ]]; then
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

###########################################################################################################
# 1 - Check database connection
###########################################################################################################

DB_RESULT=$( echo "UPDATE mg_traits.mg_traits_jobs SET time_started = now(), job_id = ${JOB_ID}, cluster_node = '${HOSTNAME}' WHERE sample_label = '${SAMPLE_LABEL}' AND id = ${MG_ID};" \
| psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" )

if [[ "$?" -ne "0" ]]; then
  email_comm "Cannot connect to database. Output:${DB_RESULT}"
  cleanup && exit 2
fi

if [[ "${DB_RESULT}" != "UPDATE 1" ]]; then
  email_comm "sample name ${SAMPLE_LABEL} is not in database Result:${DB_RESULT}"  
  cleanup && exit 2
fi

###########################################################################################################
# 2 - Create job directory
###########################################################################################################

echo "This job tmp dir: ${THIS_JOB_TMP_DIR}";  

# rm -r ${THIS_JOB_TMP_DIR}  # CHANGE THIS FOR REAL DATA!!!!!!!!!! 
# qdel -u megxnet  # CHANGE THIS FOR REAL DATA!!!!!!!!!! 
# echo "UPDATE mg_tratis.mg_traits_jobs  SET return_code = 130 WHERE return_code = -1;" \
#  | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"

# rm -r /bioinf/projects/megx/scratch/mg-traits/running_jobs/job-83*  # CHANGE THIS FOR REAL DATA

mkdir "${THIS_JOB_TMP_DIR}" && cd "${THIS_JOB_TMP_DIR}"
mkdir "${THIS_JOB_TMP_DIR_DATA}" && mkdir "${SINA_LOG_DIR}"

echo "Logs, data and temp files will be written to:$(pwd)"
if [[ "$(pwd)" != "${THIS_JOB_TMP_DIR}" ]]; then 
 email_comm  "Could not access job temp dir ${THIS_JOB_TMP_DIR} in $(pwd)"
 db_error_comm "Could not access job temp dir ${THIS_JOB_TMP_DIR}"
 
 cleanup && exit 2; 
fi

###########################################################################################################
# 3 - Download file from MG URL
###########################################################################################################

# validate MG URL 
echo "${MG_URL}"
REGEX='(https?|ftp|file)://[-A-Za-z0-9\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\+&@#/%=~_|]'
  
if [[ ! ${MG_URL} =~ ${REGEX} ]]; then
  email_comm "Invalid URL ${MG_URL} output db: ${DB_COM} ${SAMPLE_LABEL}"
  db_error_comm "Not valid URL: ${MG_URL}";
  cleanup && exit 1
fi 

# check if it already exist on our DB            

if [[ "${SAMPLE_LABEL}" != "test_label" ]]; then
  URLDB=$(psql -t -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" -c \
  "SELECT count(*) FROM mg_traits.mg_traits_jobs where mg_url = '${MG_URL}' AND sample_label NOT ILIKE 'test_label AND return_code = 0'")
    
  if [[ "${URLDB}" -gt 1 ]]; then 
     email_comm "The URL "${MG_URL}" has been already succesfully crunched. If the file is different please change the file name"
     db_error_comm "The URL ${MG_URL} has been already succesfully crunched. If the file is different please change the file name."
     cleanup && exit 1
  fi 
fi

# # download MG_URL
printf "Downloading ${MG_URL} to ${RAW_DOWNLOAD}..."
curl -s "${MG_URL}" > "${RAW_DOWNLOAD}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "Could not retrieve ${MG_URL}"
  db_error_comm  "Could not retrieve ${MG_URL}"
  cleanup && exit 1
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
  
  cleanup && exit 1
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
  cleanup && exit 2; 
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
  cleanup && exit 1; 
fi


# tf_file download
echo "${TFFILE_URL}"
curl -s "${TFFILE_URL}" > "${TFFILE}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "Could not retrieve ${TFFILE_URL}"
  db_error_comm "Could not retrieve ${TFFILE_URL}"
  cleanup && exit 1;
fi
  

# slv download
echo "${SLV_TAX_URL}"
curl -s "${SLV_TAX_URL}" > "${SLV_FILE}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "Could not retrieve ${SLV_TAX_URL}"
  db_error_comm "Could not retrieve ${SLV_TAX_URL}"
  cleanup && exit 1; 
fi

###########################################################################################################
# 7 - Check for duplicates
###########################################################################################################

#printf "Removing duplicated sequences..."
#qsub  -l h="mg9.mpi-bremen.de|mg10.mpi-bremen.de|mg11.mpi-bremen.de|mg12.mpi-bremen.de|mg13.mpi-bremen.de|mg14.mpi-bremen.de|mg15.mpi-bremen.de|mg16.mpi-bremen.de,exclusive" \
#-sync y -pe threaded $NSLOTS  "${vsearch_runner}" "${RAW_FASTA}" "${UNIQUE}" "${UNIQUE_LOG}" "${NSLOTS}"

#VSEARCH_ERROR_CODE="$?"

#rm -r /bioinf/projects/megx/scratch/mg-traits/failed_jobs/job*
#rm -r /bioinf/projects/megx/scratch/mg-traits/running_jobs/job*

#if [[ "${CD_HIT_ERROR_CODE}" -ne "0" ]]; then 
#  email_comm "${cd_hit_dup} -i ${RAW_FASTA} -o /dev/null > ${UNIQUE_LOG}
#exited with RC ${CD_HIT_ERROR_CODE} in job ${JOB_ID}\nFiles are at: ${FAILED_JOBS_DIR}/job-${JOB_ID}"
#  db_error_comm "${MG_URL} could not be processed by cd-hit-dup"
#  cleanup && exit 2;
#fi

#### ONLY FOR TARA!!!! ######
MG_URL_LOG=$( echo "${MG_URL}" | sed 's/pre-process.SR.*.fasta/pre-process.SR_vsearch.log/')
curl -s "${MG_URL_LOG}" > pre-process.SR_vsearch.log
NUM_READS=$( sed -n 3p pre-process.SR_vsearch.log | cut -f10 -d" " )
#### ONLY FOR TARA!!!! ######


#NUM_READS=$(grep 'Total number of sequences:'  "${UNIQUE_LOG}" | awk '{print $(NF)}')
#NUM_UNIQUE=$(grep 'Number of clusters found:'  "${UNIQUE_LOG}" | awk '{print $(NF)}')

#echo "Number of sequences: ${NUM_READS}"
#echo "Number of unique sequences: ${NUM_UNIQUE}"

#if [[ "$NUM_READS" -ne "$NUM_UNIQUE" ]]; then 
#  email_comm "${MG_URL} contains duplicates. Please provide a pre-processed metagenome."
#  db_error_comm "${MG_URL} contains duplicates"
#  cleanup && exit 1; 
#fi

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
  cleanup && exit 2; 
fi

# seq stats
Rscript --vanilla "${seq_stats}" "${INFOSEQ_TMPFILE}" "${INFOSEQ_MGSTATS}"
SEQ_STATS_ERROR_CODE="$?"

if [[ "${SEQ_STATS_ERROR_CODE}" -ne "0" ]]; then 
  email_comm "Rscript --vanilla ${seq_stats} ${INFOSEQ_TMPFILE} ${INFOSEQ_MGSTATS}
exited with RC ${SEQ_STATS_ERROR_CODE} in job ${JOB_ID}. Infoseq script. Files are at: ${FAILED_JOBS_DIR}/job-${JOB_ID}"
  db_error_comm "Cannot process sequence statistics. Please contact adminitrator."
  cleanup && exit 2; 
fi

NUM_BASES=$(cut -f1 "${INFOSEQ_MGSTATS}" -d ' '); 
GC=$(cut -f2 "${INFOSEQ_MGSTATS}" -d ' '); 
VARGC=$(cut -f3 "${INFOSEQ_MGSTATS}" -d ' ')
printf "Number of bases: %d\nGC content: %f\nGC variance: %f\n" "${NUM_BASES}" "${GC}" "${VARGC}"


###########################################################################################################
# define environment for sub jobs
###########################################################################################################

cat > 01-job_env << EOF 
JOB_ID="${JOB_ID}"
SAMPLE_LABEL="${SAMPLE_LABEL}"
MG_ID="${MG_ID}"
RAW_FASTA="${RAW_FASTA}"
GC="${GC}"
VARGC="${VARGC}"
NUM_BASES="${NUM_BASES}"
NUM_READS="${NUM_READS}"
RUNNING_JOBS_DIR="${RUNNING_JOBS_DIR}"
MG_URL="${MG_URL}"
EOF

###########################################################################################################
# 1 - run fgs
###########################################################################################################

#Split original
awk -vn="${NSEQ}" 'BEGIN {n_seq=0;partid=1;} /^>/ {if(n_seq%n==0){file=sprintf("05-part-%d.fasta",partid);partid++;} print >> file; n_seq++; next;} { print >> file; }' < "${RAW_FASTA}"
NFILES=$(ls -1 05-part*.fasta | wc -l)

qsub  -t 1-"${NFILES}" -pe threaded "${NSLOTS}" -N "${FGS_JOBARRAYID}" ${fgs_runner}
# "${fgs_runner}" "${NSLOTS}" "${NFILES}" "${FGS_JOBARRAYID}"

ERROR_FGS=$?

if [[ "${ERROR_FGS}" -ne "0" ]]; then
  email_comm  "${frag_gene_scan} -genome=${IN_FASTA_FILE} -out=${IN_FASTA_FILE}.genes10 -complete=0 -train=sanger_10
exited with RC ${ERROR_FGS} in job ${JOB_ID}"
  db_error_comm "FragGeneScan failed. Please contact adminitrator."
  cleanup && exit 2
fi

############################################################################################################
## 2 - run sortmerna
############################################################################################################

#MEM=$(free -m | grep Mem | awk '{printf "%d",$2/3}')
MEM=4000
"${sortmerna}" --reads "${RAW_FASTA}" -a "${NSLOTS}" --ref \
"${DB}"/rRNA_databases/silva-bac-16s-id90.fasta,\
"${DB}"/index/silva-bac-16s-db:"${DB}"/rRNA_databases/silva-arc-16s-id95.fasta,\
"${DB}"/index/silva-arc-16s-db:"${DB}"/rRNA_databases/silva-euk-18s-id95.fasta,\
"${DB}"/index/silva-euk-18s-db --blast 1 --fastx --aligned "${SORTMERNA_OUT}" -v --log -m "${MEM}" --best 1 > sortmerna.log

if [[ "${ERROR_SORTMERNA}" -ne "0" ]]; then
  email_comm "${sortmerna} --reads ${RAW_FASTA} -a ${NSLOTS} --ref ${DB}/rRNA_databases/silva-bac-16s-id90.fasta ...
exited with RC ${ERROR_SORTMERNA} in job ${JOB_ID}."
  db_error_comm "sortmerna failed. Please contact adminitrator"
  cleanup && exit 2
fi

NUM_RNA=$(egrep -c ">" "${SORTMERNA_OUT}".fasta)

if [[ "${NUM_RNA}" -eq "0" ]]; then
  email_comm "not RNA sequence found by sortmerna"
  db_error_comm "no RNA sequence found by sortmerna"
  cleanup && exit 2
fi  


###########################################################################################################
# 3 - run SINA
###########################################################################################################

awk -vn="${nSEQ}" 'BEGIN {n_seq=0;partid=1;} /^>/ {if(n_seq%n==0){file=sprintf("06-part-%d.fasta",partid);partid++;} print >> file; n_seq++; next;} { print >> file; }' < "${SORTMERNA_OUT}".fasta
nFILES=$(ls -1 06-part*.fasta | wc -l)

qsub -pe threaded "${NSLOTS}" -t 1-"${nFILES}" -N "${SINA_JOBARRAYID}" "${sina_runner}"
# "${sina_runner}" "${NSLOTS}" "${nFILES}" "${SINA_JOBARRAYID}"

ERROR_SINA=$?

if [[ "${ERROR_SINA}" -ne "0" ]]; then
  email_comm "qsub -pe threaded ${NSLOTS} -t 1-${nFILES} -N ${SINA_JOBARRAYID} ${sina_runner} failed
exited with RC ${ERROR_SINA} in job ${JOB_ID}."
  db_error_comm "sina failed. Please contact adminitrator"
  cleanup && exit 2
fi

###########################################################################################################
# 4 - run finish traits
###########################################################################################################

qsub -sync y -pe threaded "${NSLOTS}" -l h=\!mg32 -N "${FINISHJOBID}" -o "${THIS_JOB_TMP_DIR}" -e "${THIS_JOB_TMP_DIR}" -l ga -j y -terse -P megx.p -R y -m sa -M "${mt_admin_mail}" \
-hold_jid "${FGS_JOBARRAYID}","${SINA_JOBARRAYID}"  /bioinf/projects/megx/mg-traits/resources/bin/finish_runner.sh "${THIS_JOB_TMP_DIR}"

if [[ "$?" -ne "0" ]]; then
  email_comm "qsub finish_runner.sh failed"
  db_error_comm "qsub finish_runner.sh failed"
fi 

###########################################################################################################
# 5 - finished job communication: update mg_traits_jobs
###########################################################################################################

END_TIME=`date +%s.%N`
RUN_TIME=`echo "${END_TIME}"-"${START_TIME}" | bc -l`
# mv $THIS_JOB_TMP_DIR $FINISHED_JOBS_DIR


echo "UPDATE mg_traits.mg_traits_jobs SET total_run_time = total_run_time + "${RUN_TIME}", time_protocol = time_protocol \
|| ('${JOB_ID}', 'mg_traits', ${RUN_TIME})::mg_traits.time_log_entry WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
| psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"

###########################################################################################################
# 6 remove preprcess data 
###########################################################################################################

# FILE=$( echo "SELECT mg_url FROM mg_traits.mg_traits_jobs WHERE label=' ${SAMPLE_LABEL}' AND id = '${MG_ID}'" \ |
# | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" )
# FILE=$(echo $FILE | sed 's/file:\/\///')
# rm $FILE


