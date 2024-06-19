# PARAMS
JOB_NAME="tr_r366_381_flexible"
OUTPUT_DIR="/home/harsh/code/RFdiffusion_jobs/outputs/fresh_outs/$JOB_NAME"
INPUT_PDB="/home/harsh/code/structures_out/mrp_2_complex.pdb"
CONTIGS="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-365/5-25/B381-600]"
NUM_DESIGNS=100
DIFF_STEPS=50
HOTSPOT_RES="[A318,A320,A322,A324]"

if [ -d "$OUTPUT_DIR" ]; then
    echo "$OUTPUT_DIR exists."
else
    echo "Creating $OUTPUT_DIR."
    mkdir $OUTPUT_DIR
fi

# Ensure conda is initialized
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate SE3nv

cd /home/harsh/code/RFdiffusion/
scripts/run_inference.py inference.output_prefix=$OUTPUT_DIR/$JOB_NAME inference.input_pdb=$INPUT_PDB "contigmap.contigs=$CONTIGS" inference.num_designs=$NUM_DESIGNS "ppi.hotspot_res=$HOTSPOT_RES"