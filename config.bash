## folders

TEST_DIR=/bioinf/home/epereira/workspace/mg-traits/mg_traits_tests/

mg_traits_dir="${TEST_DIR}/bioinf/projects/megx/mg-traits/bin"
temp_dir="${TEST_DIR}/bioinf/projects/megx/scratch/mg-traits"
job_out_dir="${TEST_DIR}/vol/tmp/megx/"

RUNNING_JOBS_DIR="${temp_dir}/running_jobs/"
FAILED_JOBS_DIR="${temp_dir}/failed_jobs/"
FINISHED_JOBS_DIR="${temp_dir}/finished_jobs/"

THIS_JOB_TMP_DIR=$(readlink -m "$RUNNING_JOBS_DIR/job-$JOB_ID")
THIS_JOB_TMP_DIR_DATA=$THIS_JOB_TMP_DIR/data/
SINA_LOG_DIR=$THIS_JOB_TMP_DIR/sina_log

## mails
mt_admin_mail="epereira@mpi-bremen.de"

## proxy
export http_proxy=http://webproxy.mpi-bremen.de:3128
export https_proxy=https://webproxy.mpi-bremen.de:3128

## cdhit
cd_hit_dup="/bioinf/software/cd-hit/cd-hit-4.6/cd-hit-dup"
cd_hit_dup_version="0.5"
cd_hit_est="/bioinf/software/cd-hit/cd-hit-4.6/cd-hit-est"
cd_hit_mms="/bioinf/software/cd-hit/cd-hit-4.6/make_multi_seq.pl"
cd_hit_version="4.6.1"

## fgs
frag_gene_scan="/bioinf/software/fraggenescan/fraggenescan-1.19/run_FragGeneScan.pl"
frag_gene_scan_version="1.19"

## uproc
uproc_version="1.1.2"
uproc=/bioinf/software/uproc/uproc-1.1/bin/uproc-dna
uproc_pfam="/vol/biodb/uproc/pfam27"
uproc_pfam_version="pfam27"
uproc_model="/vol/biodb/uproc/model"

## sina
sina="/arb/software/arb_ubuntu_1004/latest/sina-1.3.0/sina"
sina_arb_pt_server="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib/arb_pt_server"
sina_version="1.2.11"
sina_seed="/local/biodb/mg-traits/sina/ssu_seed_50_26_05_13_cut_t.arb"
sina_seed_version="ssu_seed_50_26_05_13_cut_t"
sina_ref="/local/biodb/mg-traits/sina/ssuref_silva_nr99_115_20_07_13.arb"
sina_ref_version="ssuref_silva_nr99_115_20_07_13.arb"
ARBHOME="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/"
LD_LIBRARY_PATH="/bioinf/projects/megx/mg-traits/bin/sina-1.2.13/lib:/bioinf/software/gcc/gcc-4.9/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH

## R
r_interpreter="/bioinf/software/R/R-3.1.2/bin/R"
r_interpreter_version="3.1.2"

## URLs: pfam, silva, tf_file
PFAM_ACCESSIONS_URL="https://colab.mpi-bremen.de/micro-b3/svn/analysis-scripts/trunk/mg-traits/data/pfam27_acc.txt"
TFFILE_URL="https://colab.mpi-bremen.de/micro-b3/svn/analysis-scripts/trunk/mg-traits/data/TF.txt"
SLV_TAX_URL="https://colab.mpi-bremen.de/micro-b3/svn/analysis-scripts/trunk/mg-traits/data/silva_tax_order_115.txt"
PFAM_ACCESSIONS=$THIS_JOB_TMP_DIR/data/pfam27_acc.txt
TFFILE=$THIS_JOB_TMP_DIR/data/TF.txt
SLV_FILE=$THIS_JOB_TMP_DIR/data/silva_tax_order_115.txt

## names
RAW_DOWNLOAD="01-raw-download"
RAW_FASTA="01-raw-fasta"
BAD_FASTA="01-bad-fasta"
UNIQUE="02-unique-sequences"
UNIQUE_LOG="02-unique-sequences.log"
CLUST95="03-clustered-sequences"
CLUST95_LOG=$CLUST95".log"
CLUST95_CLSTR=$CLUST95".clstr"
INFOSEQ_TMPFILE="04-stats-tempfile"
INFOSEQ_MGSTATS="04-mg_stats"
FGS_JOBARRAYID="mt-$JOB_ID-fgs"
SINA_JOBARRAYID="mt-$JOB_ID-sina"
FINISHJOBID="mt-$JOB_ID-finish"
TMP_VOL_FILE=/vol/tmp/megx/${JOB_NAME}.${JOB_ID}
NSEQ=500

##########################
### functions
#########################

function cleanup {
if [ -f ${TMP_VOL_FILE} ];then
mv ${TMP_VOL_FILE} $THIS_JOB_TMP_DIR
fi
}

trap cleanup EXIT ERR SIGINT SIGKILL SIGTERM


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
    echo $ERROR_UTILITIES
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
    echo $ERROR_WDIRECTORIES                                                                  
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
    echo $ERROR_RDIRECTORIES                                                                   
}
