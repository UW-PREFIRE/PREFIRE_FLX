#!/usr/bin/env bash 

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./run.sh    OR    bash run.sh

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run.sh; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

absfpath() {
  # Generate absolute filepath from a relative (or even an absolute) filepath.
  #
  # Based on (circa Oct 2023) https://stackoverflow.com/questions/3915040/how-to-obtain-the-absolute-path-of-a-file-via-shell-bash-zsh-sh
  # 
  # $1     : a relative (or even an absolute) filepath
  # Returns the corresponding absolute filepath.
  if [ -d "$1" ]; then
    # dir
    (cd "$1"; pwd)
  elif [ -f "$1" ]; then
    # file
    if [[ $1 = /* ]]; then
      echo "$1"
    elif [[ $1 == */* ]]; then
      echo "$(cd "${1%/*}"; pwd)/${1##*/}"
    else
      echo "$(pwd)/$1"
    fi
  fi
}

activate_conda_env () {
  . "$1"/bin/activate;
}

deactivate_conda_env () {
  . "$1"/bin/deactivate;
}

#set -ve;  # Exit on the first error, and print out commands as we execute them
set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
readonly base_dir="$(absfpath ".")";

hn=`hostname -s`;  # Hostname

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

non_SDPS_hostname="longwave";

input_dir="${base_dir}/inputs";

# PROC_MODE=1
#cfg_str1="ATRACK_IDXRANGE_0BASED_INCLUSIVE:5800:5900,${input_dir}/PREFIRE_SAT1_1B-RAD_R01_P00_20241007075724_01877.nc||${input_dir}/PREFIRE_SAT1_AUX-MET_R01_P00_20241007075724_01877.nc||${input_dir}/PREFIRE_SAT1_2B-MSK_R01_P00_20241007075724_01877.nc||${input_dir}/PREFIRE_SAT1_2B-CLD_R01_P00_20241007075724_01877.nc";
cfg_str2="ATRACK_IDXRANGE_0BASED_INCLUSIVE:5400:6000,${input_dir}/PREFIRE_SAT2_1B-RAD_R01_P00_20241007071543_02040.nc||${input_dir}/PREFIRE_SAT2_AUX-MET_R01_P00_20241007071543_02040.nc||${input_dir}/PREFIRE_SAT2_2B-MSK_R01_P00_20241007071543_02040.nc||${input_dir}/PREFIRE_SAT2_2B-CLD_R01_P00_20241007071543_02040.nc";


# Specify that numpy, scipy, et cetera should not use more than one thread or
#  process):
MKL_NUM_THREADS=1;
NUMEXPR_NUM_THREADS=1;
OMP_NUM_THREADS=1;
VECLIB_MAXIMUM_THREADS=1;
OPENBLAS_NUM_THREADS=1;
export MKL_NUM_THREADS NUMEXPR_NUM_THREADS OMP_NUM_THREADS;
export VECLIB_MAXIMUM_THREADS OPENBLAS_NUM_THREADS;


# Some environment vars that convey configuration info to the algorithm:

this_top_dir="$(absfpath "${base_dir}/..")";

PACKAGE_TOP_DIR="${this_top_dir}";
ANCILLARY_DATA_DIR="${this_top_dir}/dist/ancillary";

#====== For certain types of development and testing ======
#SRF_DISAMBIG_STR="SRF_v12_2023-08-09";

#export SRF_DISAMBIG_STR;
#==========================================================

output_dir_base="${base_dir}/outputs";

# Processing mode #1 (PROC_MODE=1): Create part of a 2B-FLX granule using
#                                   AUX-MET preliminary merged surface type

  # * Only increment 'Rxx' when the resulting products will be DAAC-ingested
PRODUCT_FULLVER="R01_P00";
  # Special form ('R00_Syy') when processing simulated observations:
#PRODUCT_FULLVER="R00_S07";

# Make required environment vars available:
export PACKAGE_TOP_DIR ANCILLARY_DATA_DIR;
export PRODUCT_FULLVER;

# Check if output file directory exists; if not, bail:
tmpdir="${output_dir_base}";
test -d "${tmpdir}" || { echo "Output directory does not exist: ${tmpdir}"; exit 1; }

# If custom conda environment files exist, activate that conda environment:
conda_env_dir="${this_top_dir}/dist/c_env_for_PREFIRE_FLX";
if [ -d "${conda_env_dir}" ]; then
   activate_conda_env "${conda_env_dir}";
fi

# Execute any necessary machine setup instructions ((un)loading modules, etc.):
if [ "x$hn" = "x$non_SDPS_hostname" ]; then
   . "${this_top_dir}/dist/perform_machine_setup.sh";
fi

OUTPUT_DIR=${output_dir_base}/m1;
PROC_MODE=1;
export OUTPUT_DIR PROC_MODE;
   # Execute script that writes a new 'prdgit_version.txt', which contains
   #  product moniker(s) and current (latest) git hash(es) that are part of the
   #  provenance of this package's product(s).
   # *** This step should not be done within the SDPS, since that file is
   #     created just before delivery to the SDPS.
if [ ! -f "${this_top_dir}/dist/for_SDPS_delivery.txt" ]; then
   python "${this_top_dir}/dist/determine_prdgit.py";
fi
for cfg_str in ${cfg_str1} ${cfg_str2}
do
   ATRACK_IDX_RANGE_0BI=${cfg_str%,*};

   tmp_str=${cfg_str##*,};  # L1B, AUX-MET, 2B-MSK, 2B-CLD, [AUX-SAT]
   tmpA_str=${tmp_str%||*};  # L1B, AUX-MET, 2B-MSK, [2B-CLD]
   tmpB_str=${tmp_str##*||};  # 2B-CLD? or AUX-SAT?
   if [[ "$tmpB_str" == *"AUX-SAT"* ]]; then
      AUX_SAT_FILE=$tmpB_str;
      export AUX_SAT_FILE;
      L2B_CLD_FILE=${tmpA_str##*||};
      tmp_str=${tmpA_str%||*};  # L1B, AUX-MET, 2B-MSK
   else
      L2B_CLD_FILE=$tmpB_str;
      tmp_str=$tmpA_str;  # L1B, AUX-MET, 2B-MSK
   fi
   tmpA_str=${tmp_str%||*};  # L1B, AUX-MET
   L2B_MSK_FILE=${tmp_str##*||};
   L1B_RAD_FILE=${tmpA_str%||*};
   AUX_MET_FILE=${tmpA_str##*||};
   export ATRACK_IDX_RANGE_0BI L1B_RAD_FILE AUX_MET_FILE L2B_MSK_FILE;
   export L2B_CLD_FILE;

   if [ "x$1" = "x-i" ]; then
      python "${this_top_dir}"/dist/produce_L2.py -i;
   else
      if [ "x$hn" = "x$non_SDPS_hostname" ]; then
         #!!! The MATLAB module on longwave does not seem to set these needed
         #     paths (for standalone M-executable execution) in any way, so do
         #     it here:
         MATLAB_BASEDIR=/opt/matlab/r2021b;
         LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MATLAB_BASEDIR/bin/glnxa64:$MATLAB_BASEDIR/runtime/glnxa64";
         export LD_LIBRARY_PATH;
      fi

      python "${this_top_dir}"/dist/produce_L2.py $@;
   fi
done

# If custom conda environment files exist, DEactivate that conda environment:
if [ -d "${conda_env_dir}" ]; then
   deactivate_conda_env "${conda_env_dir}";
fi
