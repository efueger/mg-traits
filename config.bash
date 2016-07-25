NSLOTS=4

## folders
mg_traits_dir="/bioinf/projects/megx/mg-traits/bin"
temp_dir="/bioinf/projects/megx/scratch/mg-traits"
job_out_dir="/vol/tmp/megx/"

RUNNING_JOBS_DIR="${temp_dir}/running_jobs/"
FAILED_JOBS_DIR="${temp_dir}/failed_jobs/"
FINISHED_JOBS_DIR="${temp_dir}/finished_jobs/"
THIS_JOB_TMP_DIR=$(readlink -m "${RUNNING_JOBS_DIR}/job-${JOB_ID}")
THIS_JOB_TMP_DIR_DATA="${THIS_JOB_TMP_DIR}/data/"
SINA_LOG_DIR="${THIS_JOB_TMP_DIR}/sina_log"

## mails
mt_admin_mail="epereira@mpi-bremen.de"

## cdhit
cd_hit_dup="/bioinf/projects/megx/mg-traits/resources/bin/cd-hit-otu-illumina-0.0.1/cd-hit-dup-0.0.1-2011-09-30/cd-hit-dup"
#cd_hit_dup="/bioinf/projects/megx/mg-traits/resources/bin/cdhit-master/cd-hit-auxtools/cd-hit-dup"cd_hit_dup_version="4.6"
cd_hit_est="/bioinf/projects/megx/mg-traits/resources/bin/cdhit-master/cd-hit-est"
cd_hit_mms="/bioinf/projects/megx/mg-traits/resources/bin/cdhit-master/make_multi_seq.pl"
cd_hit_version="4.6"

## fgs
frag_gene_scan="/bioinf/software/fraggenescan/fraggenescan-1.19/run_FragGeneScan.pl"
frag_gene_scan_version="1.19"

## uproc
uproc_version="1.2"
uproc="/bioinf/software/uproc/uproc-1.2/bin/uproc-dna"
uproc_pfam="/local/biodb/uproc/pfam28"
uproc_pfam_version="pfam28"
uproc_model="/local/biodb/uproc/model"

## sortmerna
sortmerna="/bioinf/software/sortmerna/sortmerna-2.0/bin/sortmerna"
DB="/bioinf/software/sortmerna/sortmerna-2.0/"

## sina                                
sina="/bioinf/software/sina/sina-1.3.0rc/sina"
sina_arb_pt_server="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib/arb_pt_server"
sina_version="1.2.11"

sina_seed="/local/biodb/mg-traits/sina/ssu_seed_50_26_05_13_cut_t.arb"  # SINA SEED HAS TO BE UPDATED????!!!!!
sina_seed_version="ssu_seed_50_26_05_13_cut_t"

#sina_ref="/local/biodb/mg-traits/sina/ssuref_silva_nr99_115_20_07_13.arb"
#sina_ref="/bioinf/projects/megx/mg-traits/resources/sina/SSURef_Nr99_123.1_SILVA_03_03_16_opt.arb" NEEDS PERMISSIONS TO RUN IN /bioinf/projects/megx !!!!
sina_ref="/bioinf/home/epereira/workspace/mg-traits/resources/sina/SSURef_Nr99_123.1_SILVA_03_03_16_opt.arb"
sina_ref_version="ssuref_silva_nr99_123_03_03_16.arb"

SINA_SOCKET=":/tmp/mg_traits_pt_"$(tr -cd '[:alnum:]' < /dev/urandom | fold -w 32 | head -n 1)

ARBHOME="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/"
LD_LIBRARY_PATH="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib:/bioinf/software/gcc/gcc-4.9/lib64:/usr/lib/libgomp.so.1:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH

## R
r_interpreter="/bioinf/software/R/R-3.2.3/bin/R"
r_interpreter_version="3.1.3"

## scripts
fasta_file_check="/bioinf/projects/megx/mg-traits/resources/bin/fasta_file_check.pl"
fgs_runner="/bioinf/projects/megx/mg-traits/resources/bin/fgs_runner2.sh"
sortmerna_runner="/bioinf/projects/megx/mg-traits/resources/bin/sortmerna_runner.sh"
sina_runner="/bioinf/projects/megx/mg-traits/resources/bin/sina_runner2.sh"
seq_stats="/bioinf/projects/megx/mg-traits/resources/bin/seq_stats.R"
cd_hit_dup_runner="/bioinf/projects/megx/mg-traits/resources/bin/cd_hit_dup_runner.sh"
vsearch_runner="/bioinf/projects/megx/mg-traits/resources/bin/vsearch_runner.sh"

## URLs: pfam, silva, tf_file
#PFAM_ACCESSIONS_URL="https://colab.mpi-bremen.de/micro-b3/svn/analysis-scripts/trunk/mg-traits/data/pfam27_acc.txt"
PFAM_ACCESSIONS_URL="file:///bioinf/projects/megx/mg-traits/resources/pfam/pfam28_acc.txt"
TFFILE_URL="https://colab.mpi-bremen.de/micro-b3/svn/analysis-scripts/trunk/mg-traits/data/TF.txt"
SLV_TAX_URL="https://colab.mpi-bremen.de/micro-b3/svn/analysis-scripts/trunk/mg-traits/data/silva_tax_order_115.txt"
PFAM_ACCESSIONS=$THIS_JOB_TMP_DIR/data/pfam28_acc.txt
TFFILE=$THIS_JOB_TMP_DIR/data/TF.txt
SLV_FILE=$THIS_JOB_TMP_DIR/data/silva_tax_order_115.txt

## names
RAW_DOWNLOAD="01-raw-download"
RAW_FASTA="01-raw-fasta"
FASTA_BAD="01-bad-fasta"
UNIQUE="02-unique-sequences"
UNIQUE_LOG="02-unique-sequences.log"
CLUST95="03-clustered-sequences"
CLUST95_LOG=$CLUST95".log"
CLUST95_CLSTR=$CLUST95".clstr"
INFOSEQ_TMPFILE="04-stats-tempfile"
INFOSEQ_MGSTATS="04-mg_stats"
SORTMERNA_OUT="06-sortmerna"
FGS_JOBARRAYID="mt-$JOB_ID-fgs"
SINA_JOBARRAYID="mt-$JOB_ID-sina"
FINISHJOBID="mt-$JOB_ID-finish"
TMP_VOL_FILE=/vol/tmp/megx/${JOB_NAME}.${JOB_ID}
NSEQ=2000000 # set to 2000 000 
nSEQ=2000  # set to 1000

##########################
### functions
#########################

function cleanup {
if [[ -f ${TMP_VOL_FILE} ]];then
mv "${TMP_VOL_FILE}" "${THIS_JOB_TMP_DIR}"
fi
mv "${THIS_JOB_TMP_DIR}" "${FAILED_JOBS_DIR}"
}


function check_required_programs() {

    req_progs=("$@"); 
    ERROR_UTILITIES=""
   
    for p in ${req_progs[@]}; do
        hash "${p}" 2>&- || \
             if  [[ -z "${ERROR_UTILITIES}" ]]; then 
	        ERROR_UTILITIES="${p}"
	      else
	        ERROR_UTILITIES="${ERROR_UTILITIES} ${p}"   
	     fi   	      
    done
    echo "${ERROR_UTILITIES}"
}


function check_required_writable_directories() {

    req_dirs=("$@"); 
    ERROR_WDIRECTORIES=""
   
    for d in ${req_dirs[@]}; do
        if [[ ! -d "${d}" && ! -w "${d}" ]]; then 
          if [[ -z "${ERROR_WDIRECTORIES}" ]]; then 
            ERROR_WDIRECTORIES="${d}"
	  else
	    ERROR_WDIRECTORIES="${ERROR_WDIRECTORIES} ${d}"   
 	  fi 	  
        fi         
    done 
    echo "${ERROR_WDIRECTORIES}"
}


function check_required_readable_directories() {

    req_dirs=("$@"); 
    ERROR_RDIRECTORIES=""
   
    for d in ${req_dirs[@]}; do
        if [[ ! -d "${d}" && ! -r "${d}" ]]; then 
          if [[ -z "${ERROR_RDIRECTORIES}" ]]; then 
            ERROR_RDIRECTORIES="${d}"
	  else
	    ERROR_RDIRECTORIES="${ERROR_RDIRECTORIES} ${d}"   
 	  fi
        fi 
    done
    echo "${ERROR_RDIRECTORIES}"                                                                 
}


function email_comm() {
mail -s "mg_traits:${JOB_ID} failed" "${mt_admin_mail}" <<EOF
"${1}"
EOF
}

trap cleanup SIGINT SIGKILL SIGTERM

