#!/bin/bash -l
set -x
set -e
set -o pipefail
set -o errexit
set -o errtrace
set -o nounset

module load sortmerna

NAM=$1
NSLOTS=$2
RES=$3

SE=${NAM}.SR.qc.fasta

sortmerna="/bioinf/software/sortmerna/sortmerna-2.0/bin/sortmerna"
DB="/bioinf/software/sortmerna/sortmerna-2.0/"
#MEM=$(free -m | grep Mem | awk '{printf "%d",$2/3}')
MEM=4000
$sortmerna --reads "${SE}" -a ${NSLOTS} --ref ${DB}/rRNA_databases/silva-bac-16s-id90.fasta,${DB}/index/silva-bac-16s-db:${DB}/rRNA_databases/silva-arc-16s-id95.fasta,${DB}/index/silva-arc-16s-db:${DB}/rRNA_databases/silva-euk-18s-id95.fasta,${DB}/index/silva-euk-18s-db --blast 1 --fastx --aligned "${RES}/${NAM}.sortmerna.rDNA" -v --log -m ${MEM} --best 1

  