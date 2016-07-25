#!/bin/bash
#$ -j y
#$ -cwd

set -x
set -o pipefail

#####################################
## define variables 
#####################################
START_TIME=`date +%s.%N`

SAMPLE_LABEL="${1}"
NSLOTS="${2}"
PREPROCESS_DIR="/bioinf/projects/megx/mg-traits/TARA_crunch/preprocess_data"
WORKING_TMP_DIR="${PREPROCESS_DIR}/${SAMPLE_LABEL}"

#NSLOTS=12
# pear
pear="/bioinf/software/pear/pear-0.9.8/bin/pear"
# bbduk
bbduk="/bioinf/software/bbmap/bbmap-35.14/bbduk.sh"
# vsearch
vsearch_runner="/bioinf/projects/megx/mg-traits/resources/bin/vsearch_runner.sh"
# db communication
target_db_user=epereira
target_db_host=antares
target_db_port=5434
target_db_name=megdb_r8

# email
mt_admin_mail=epereira@mpi-bremen.de

# declare AFILE="/bioinf/projects/megx/TARA/rDNAs/compute/input/assem-file-name2url.txt"
declare AFILE="/bioinf/projects/megx/mg-traits/resources/assem-file-name2url.txt"

######################################
# Functions
######################################

function email_comm() {
mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${1}"
EOF
}


function db_error_comm() {
  echo "UPDATE epereira.preprocess_jobs SET time_finished = now(), return_code = 1, error_message = '${1}' \
  WHERE sample_label = '${SAMPLE_LABEL}' AND job_id = '${JOB_ID}';" | psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
}

function cleanup() {
mv "${WORKING_TMP_DIR}" /bioinf/projects/megx/mg-traits/TARA_crunch/failed_preprocess_data/
}
trap cleanup SIGINT SIGKILL SIGTERM

#####################################################
# insert job in table preprocess_jobs
#####################################################
#JOB_ID=1
DB_RESULT=$( echo  "INSERT INTO epereira.preprocess_jobs (customer, mg_url, sample_label, sample_env_ontology, time_started, job_id, cluster_node ) \
VALUES ('anonymous', 'file://no_file', '${SAMPLE_LABEL}', 'marine', now(),'${JOB_ID}', '${HOSTNAME}' );" \
| psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}" );

if [[ "$?" -ne "0" ]]; then
  email_comm "Cannot connect to database. Output:${DB_RESULT}"
  cleanup && exit 2
fi

if [[ "${DB_RESULT}" != "INSERT 0 1" ]]; then
  email_comm "sample name ${SAMPLE_LABEL} is not in database Result:${DB_RESULT}"  
  cleanup && exit 2
fi


#####################################################
# Map and link the TARA samples name with their file.
#####################################################

declare -r MINOV='10'
RFILES=$(grep "${SAMPLE_LABEL}" "${AFILE}")
NRFILE=$(echo "${RFILES}" | wc -l) #How many files do we have

echo "${RFILES}" | \
    cut -f 2 -d ' ' | \
    while read LINE
    do
        N=$(basename "${LINE}")
        ln -s "/bioinf/projects/megx/TARA/assemblies/FASTQ/${N}" .
        #ln -s "/bioinf/home/epereira/workspace/mg-traits/tara_prepross/data/toyFASTQ/${N}" .
    done

echo "UPDATE epereira.preprocess_jobs SET mg_url = '${RFILES}' WHERE sample_label = '${SAMPLE_LABEL}' AND  job_id = '${JOB_ID}';" | \
psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"


######################################################    
# Combine all TARA files from the same sample name.
######################################################
R1_raw_all=R1_raw_all.fastq.gz
R2_raw_all=R2_raw_all.fastq.gz
cat *_1.fastq.gz > "${R1_raw_all}"
cat *_2.fastq.gz > "${R2_raw_all}"

rm *_1.fastq.gz #Remove original file links
rm *_2.fastq.gz 

########################################################
# Remove adapters
########################################################
PE1_rmadapt="R1_rmadapt.fastq"
PE2_rmadapt="R2_rmadapt.fastq"
SE_rmadapt="SR_rmadapt.fastq"
TARAADAP="/bioinf/software/bbmap/bbmap-35.14/resources/tara.fa.gz"

${bbduk} in="${R1_raw_all}" in2="${R2_raw_all}" out1="${PE1_rmadapt}" out2="${PE2_rmadapt}" \
outs="${SE_rmadapt}" qin=33 minlen=45 ktrim=r k=25 mink=11 ref="${TARAADAP}" hdist=1 tbo tpe maxns=0 threads="${NSLOTS}"

if [[ "$?" -ne "0" ]]; then 
  email_comm "remove adaptares bbduk.sh failed"
  db_error_comm "remove adaptares bbduk.sh failed"
  cleanup && exit 2
fi

rm R1_raw_all.fastq.gz #Remove raw data
rm R2_raw_all.fastq.gz 

#######################################################
# Define memory usage. 
#######################################################
declare MEM=$(free -g | grep Mem | awk '{printf "%dG",$2/3}') #We use one third of the memory available
echo "${MEM}"

########################################################
# Run PEAR to merge the data
########################################################
${pear} -j "${NSLOTS}" -y "${MEM}" -v "${MINOV}" -f "${PE1_rmadapt}" -r "${PE2_rmadapt}" -o "pear"

if [[ "$?" -ne "0" ]]; then 
  email_comm "pear merge sequences failed"
  db_error_comm "pear merge sequences failed"
  cleanup && exit 2
fi

NDISCARD=$(echo $(wc *.discarded.fastq -l | cut -f1 -d" " ) / 4 | bc)
echo "${NDISCARD}"

rm "${PE1_rmadapt}" "${PE2_rmadapt}" #Remove the unmerged reads and discarded data
[[ -s pear.discarded.fastq ]] && rm *discarded.fastq

#########################################################
# Quality trim non merged
#########################################################

if [[ -s "pear.unassembled.forward.fastq" ]]; then
  PE1_qc_nonmerged="R1_qc_nonmerged.fasta"
  PE2_qc_nonmerged="R2_qc_nonmerged.fasta"
  
  ${bbduk}  in="pear.unassembled.forward.fastq" in2="pear.unassembled.reverse.fastq" \
  out1="${PE1_qc_nonmerged}" out2="${PE2_qc_nonmerged}" outs="tmp.SR.fasta" qin=33 minlen=45 qtrim=rl \
  trimq=20 ktrim=r k=25 mink=11 ref="${TARAADAP}" hdist=1 tbo tpe maxns=0 threads="${NSLOTS}"

  
  if [[ "$?" -ne "0" ]]; then  
    email_comm "bbduk.sh quality trim non merged failed"
    db_error_comm "bbduk.sh quality trim non merged failed"
    cleanup && exit 2
  fi  
 
  NUNMERGED=$(echo $(wc "pear.unassembled.forward.fastq" -l | cut -f1 -d" " ) / 4 | bc)
  echo "${NUNMERGED}"
  rm pear.unassembled*.fastq

fi

#########################################################
# Quality trim merged + SE_rmadapt
#########################################################

if [[ -s "pear.assembled.fastq" || -s "${SE_rmadapt}" ]]; then

  NMERGED=$(echo $(wc "pear.assembled.fastq" -l | cut -f1 -d" " ) / 4 | bc)
  echo "${NMERGED}"

  cat *assembled.fastq >> "${SE_rmadapt}"
  SR_qc="SR.qc.fasta"
  
  ${bbduk} in="${SE_rmadapt}" out1="${SR_qc}" qin=33 minlen=45 qtrim=rl trimq=20 ktrim=r k=25 mink=11 \
  ref="${TARAADAP}" hdist=1 tbo tpe maxns=0 threads="${NSLOTS}"

  if [[ "$?" -ne "0" ]]; then 
    email_comm "bbduk.sh quality trim merged + SE failed"
    db_error_comm "bbduk.sh quality trim merged + SE failed"
    cleanup && exit 2
  fi   
fi


##########################################################
# Concatenate all results
##########################################################

cat "${PE1_qc_nonmerged}" "${PE2_qc_nonmerged}" "tmp.SR.fasta" >> "${SR_qc}"
rm "tmp.SR.fasta" "${SE_rmadapt}" "${PE1_qc_nonmerged}" "${PE2_qc_nonmerged}" *assembled.fastq pear.discarded.fastq

#####################################################################
# Remove duplicates
#####################################################################

SE_nodup=pre-process.SR.fasta # Same name as PROCESS_FASTA in config.bash. It is defined here so it can be run independelty.
SE_log=pre-process.SR_vsearch.log

qsub -l h="mg9.mpi-bremen.de|mg10.mpi-bremen.de|mg11.mpi-bremen.de|mg12.mpi-bremen.de|mg13.mpi-bremen.de|mg14.mpi-bremen.de|mg15.mpi-bremen.de|mg16.mpi-bremen.de,exclusive" \
-sync y -pe threaded "${NSLOTS}" "${vsearch_runner}"  "${SR_qc}" "${SE_nodup}" "${SE_log}" "${NSLOTS}"

 if [[ "$?" -ne "0" ]]; then  
    email_comm "vsearch failed"
    db_error_comm  "vsearch failed"
    cleanup && exit 2
  fi

rm "${SR_qc}" *.clstr

END_TIME=`date +%s.%N`

#########################################################################
# time registration
#########################################################################

 echo "UPDATE mg_traits.mg_traits_jobs SET total_run_time = total_run_time + "${RUN_TIME}", time_protocol = time_protocol \
|| ('${JOB_ID}', 'mg_traits', ${RUN_TIME})::mg_traits.time_log_entry WHERE sample_label = '${SAMPLE_LABEL}' AND id = '${MG_ID}';" \
| psql -U "${target_db_user}" -h "${target_db_host}" -p "${target_db_port}" -d "${target_db_name}"
