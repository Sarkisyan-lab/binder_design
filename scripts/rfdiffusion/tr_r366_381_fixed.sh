# PARAMS
JOB_NAME="tr_r366_381_fixed_v0"
OUTPUT_DIR="../../outputs/rfdiffusion/${JOB_NAME}"
INPUT_PDB="../../inputs/mrp_2_complex.pdb"
CONTIGS="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-365/15-15/B381-600]"
NUM_DESIGNS=100
DIFF_STEPS=50
HOTSPOT_RES="[A318,A320,A322,A324]"

# Function to find the next available directory name
function next_dir {
    suffix=0
    while [ -d "${1}_${suffix}" ]; do
        let suffix+=1
    done
    echo "${1}_${suffix}"
}


# Initialize parameters with default values
RF_DIR=""

# Parse command-line arguments
while (( "$#" )); do
  case "$1" in
    --rf-dir)
      RF_DIR=$2
      shift 2
      ;;
    # Other flags...
    *)
      shift
      ;;
  esac
done


# Check if OUTPUT_DIR exists and create a new versioned directory
if [ -d "$OUTPUT_DIR" ]; then
    echo "$OUTPUT_DIR exists."
    OUTPUT_DIR=$(next_dir $OUTPUT_DIR)
    echo "Creating $OUTPUT_DIR."
    mkdir -p $OUTPUT_DIR
else
    echo "Creating $OUTPUT_DIR."
    mkdir -p $OUTPUT_DIR
fi


# cd to the RFdiffusion directory
cd $RF_DIR

# Ensure conda is initialized
# CONDA_BASE=$(conda info --base)
# source "$CONDA_BASE/etc/profile.d/conda.sh"

# conda activate SE3nv

# cd /home/harsh/code/RFdiffusion/
# scripts/run_inference.py inference.output_prefix=$OUTPUT_DIR/$JOB_NAME inference.input_pdb=$INPUT_PDB "contigmap.contigs=$CONTIGS" inference.num_designs=$NUM_DESIGNS "ppi.hotspot_res=$HOTSPOT_RES"
