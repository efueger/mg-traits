#!/bin/bash
set -x
set -o pipefail

START_TIME=`date +%s.%N`
source ~/.bashrc
source config.bash


###########################################################################################################
# 1 - run fgs
###########################################################################################################

mkdir split_qc && cd split_qc

bash ${BIN}/fgs_runner.sh $NAM $NSLOTS $NSEQ
ERROR_FGS=$?


if [[ "${ERROR_FGS}" -ne "0" ]]; then
  echo "ERROR_FGS=${ERROR_FGS}" >> tmp.vars; 
  "${BIN}"/db_commands2.2.bash;
  cd ../..; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"
  err "fgs_runner.sh failed"
  exit 2
fi

cd ../


###########################################################################################################
# 2 - run sortmerna
###########################################################################################################


bash ${BIN}/sortmerna_runner.sh $NAM $NSLOTS $RES
ERROR_SORTMERNA=$?

if [[ "${ERROR_SORTMERNA}" -ne "0" ]]; then
  echo "ERROR_SORTMERNA=${ERROR_SORTMERNA}" >> tmp.vars; 
  "${BIN}"/db_commands2.2.bash;
  cd ../; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"
  err "sortmerna_runner.sh failed"
  exit 2
fi


###########################################################################################################
# 3 - run SINA
###########################################################################################################

mkdir split_smr && cd split_smr

bash ${BIN}/sina_runner.sh $NAM $NSLOTS $nSEQ $RES
ERROR_SINA=$?

if [[ "${ERROR_SINA}" -ne "0" ]]; then
   echo "ERROR_SINA=${ERROR_SINA}" >> tmp.vars; 
  "${BIN}"/db_commands2.2.bash;
  cd ../../; mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"
  err "sina_runner.sh failed"
  exit 2
fi
cd ../



