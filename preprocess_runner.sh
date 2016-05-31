#!/bin/bash -l
set -x
set -e
set -o pipefail
#set -o errexit
set -o errtrace
set -o nounset

module load java
module load bbmap
module load pear
module load sortmerna
module load coreutils

#####################################
## define functions and variables 
#####################################

err() {
  echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: $@" >&2
}

NAM="${1}"
NSLOTS="${2}"
RES="${3}"


declare AFILE="/bioinf/projects/megx/TARA/rDNAs/compute/input/assem-file-name2url.txt"
declare -r MINOV='10'
RFILES=$(grep "${NAM}" "${AFILE}")
NRFILE=$(echo "${RFILES}" | wc -l) #How many files do we have

#####################################################
# Map and link the TARA samples name with their file.
#####################################################
echo "${RFILES}" | \
    cut -f 2 -d ' ' | \
    while read LINE
    do
        N=$(basename "${LINE}")
        #ln -s "/bioinf/projects/megx/TARA/assemblies/FASTQ/${N}" .
        ln -s "/bioinf/home/epereira/workspace/mg-traits/tara_prepross/data/toyFASTQ/${N}"
    done

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
PE1_rmadapt="${RES}/${NAM}.R1_rmadapt.fastq"
PE2_rmadapt="${RES}/${NAM}.R2_rmadapt.fastq"
SE_rmadapt="${RES}/${NAM}.SR_rmadapt.fastq"
TARAADAP=/bioinf/software/bbmap/bbmap-35.14/resources/tara.fa.gz
bbduk.sh in="${R1_raw_all}" in2="${R2_raw_all}" out1="${PE1_rmadapt}" out2="${PE2_rmadapt}" \
outs="${SE_rmadapt}" qin=33 minlen=45 ktrim=r k=25 mink=11 ref="${TARAADAP}" hdist=1 tbo tpe maxns=0 threads="${NSLOTS}"

if [[ "$?" -ne "0" ]]; then 
  err "remove adaptares bbduk.sh failed"
  exit 2
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
/bioinf/software/pear/pear-0.9.8/bin/pear -j "${NSLOTS}" -y "${MEM}" -v "${MINOV}" -f "${PE1_rmadapt}" -r "${PE2_rmadapt}" -o "${RES}/${NAM}.pear"

if [[ "$?" -ne "0" ]]; then 
  err "pear merge sequences failed"
  exit 2
fi

NDISCARD=$(echo $(wc "${RES}"/*.discarded.fastq -l | cut -f1 -d" " ) / 4 | bc)
echo "${NDISCARD}"

rm "${PE1_rmadapt}" "${PE2_rmadapt}" #Remove the unmerged reads and discarded data
[[ -s "${RES}"/${NAM}.pear.discarded.fastq ]] && rm "${RES}"/*discarded.fastq

#########################################################
# Quality trim non merged
#########################################################

if [[ -s "${RES}/${NAM}.pear.unassembled.forward.fastq" ]]; then
  PE1_qc_nonmerged="${RES}/${NAM}.R1_qc_nonmerged.fasta"
  PE2_qc_nonmerged="${RES}/${NAM}.R2_qc_nonmerged.fasta"
  bbduk.sh in="${RES}/${NAM}.pear.unassembled.forward.fastq" in2="${RES}/${NAM}.pear.unassembled.reverse.fastq" \
  out1="${PE1_qc_nonmerged}" out2="${PE2_qc_nonmerged}" outs="${RES}/tmp.SR.fasta" qin=33 minlen=45 qtrim=rl \
  trimq=20 ktrim=r k=25 mink=11 ref="${TARAADAP}" hdist=1 tbo tpe maxns=0 threads="${NSLOTS}"

  if [[ "$?" -ne "0" ]]; then 
    err "bbduk.sh quality trim unmerged sequences failed"
    exit 2
  fi  
 
  NUNMERGED=$(echo $(wc "${RES}"/"${NAM}.pear.unassembled.forward.fastq" -l | cut -f1 -d" " ) / 4 | bc)
  echo "${NUNMERGED}"
  rm "${RES}"/"${NAM}".pear.unassembled*.fastq

fi

#########################################################
# Quality trim merged + SE_rmadapt
#########################################################

if [[ -s "${RES}/${NAM}.pear.assembled.fastq" || -s "${SE_rmadapt}" ]]; then

  NMERGED=$(echo $(wc "${RES}"/"${NAM}.pear.assembled.fastq" -l | cut -f1 -d" " ) / 4 | bc)
  echo "${NMERGED}"

  cat "${RES}"/*assembled.fastq >> "${SE_rmadapt}"
  SR_qc="${RES}/${NAM}.SR.qc.fasta"
  bbduk.sh in="${SE_rmadapt}" out1="${SR_qc}" qin=33 minlen=45 qtrim=rl trimq=20 ktrim=r k=25 mink=11 \
  ref="${TARAADAP}" hdist=1 tbo tpe maxns=0 threads="${NSLOTS}"

  if [[ "$?" -ne "0" ]]; then 
    err "bbduk.sh quality trim merged sequences failed"
    exit 2
  fi   
fi

##########################################################
# Concatenate all results
##########################################################
cat "${PE1_qc_nonmerged}" "${PE2_qc_nonmerged}" "${RES}/tmp.SR.fasta" >> "${SR_qc}"
rm "${RES}/tmp.SR.fasta" "${SE_rmadapt}" "${PE1_qc_nonmerged}" "${PE2_qc_nonmerged}" "${RES}"/*assembled.fastq "${RES}"/"${NAM}".pear.discarded.fastq


