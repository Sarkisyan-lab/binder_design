import os
import argparse
import subprocess
from dotenv import load_dotenv

# load .env
load_dotenv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/.env")

# Create the parser
parser = argparse.ArgumentParser()

# Add an argument
parser.add_argument(
    "--num_designs", type=int, help="Number of RF designs to generate", default=100
)
parser.add_argument("--job_name", type=str, help="job_name", required=True)
parser.add_argument("--pdb_file_path", type=str, help="pdb file path", required=True)
parser.add_argument(
    "--hotspot_residues", type=str, help="Hotspot residues", required=True
)
parser.add_argument("--contigs", type=str, help="Contigs", required=True)
parser.add_argument("--cuda_device_id", type=int, help="Cuda device id", default=0)
parser.add_argument(
    "--logs_dir",
    type=str,
    help="Logs dir to save logs. Logs won't be saved if not provided",
    default=None,
)

# Parse the arguments
args = parser.parse_args()

# get the script file path
FILE_PATH = os.path.abspath(__file__)
BASE_PATH = os.path.dirname(os.path.dirname(FILE_PATH))


def make_latest_dir(job_name):
    dir_path = os.path.join(os.getenv("OUTPUT_DIR_PATH"), "rfdiffusion")
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    dir_contents = os.listdir(dir_path)

    created_versions = [
        dir_name for dir_name in dir_contents if dir_name.startswith(f"{job_name}_v")
    ]
    if len(created_versions) == 0:
        folder_path = os.path.join(dir_path, f"{job_name}_v0")
    else:
        latest_repo = max(
            [int(file.split("_")[-1].replace("v", "")) for file in created_versions]
        )
        folder_path = os.path.join(dir_path, f"{job_name}_v{latest_repo + 1}")

    print(f"\nMaking folder: {folder_path}\n")
    os.makedirs(folder_path)
    return folder_path


if __name__ == "__main__":

    # Load the .parser file
    num_designs = args.num_designs
    job_name = args.job_name
    hotspot_residues = args.hotspot_residues
    contigs = args.contigs
    input_pdb_path = args.pdb_file_path

    rf_diff_repo_path = os.getenv("RF_DIFFUSION_REPO_PATH")
    print(rf_diff_repo_path)

    # make output dir what's needed
    folder_path = make_latest_dir(job_name)

    commands = [
        f"source {os.getenv('CONDA_PATH')}",
        f"conda activate {os.getenv('RFDIFF_ENV')}",
        f"cd {os.getenv('RF_DIFFUSION_REPO_PATH')}",
        f"CUDA_VISIBLE_DEVICES={args.cuda_device_id}",
        f'scripts/run_inference.py inference.output_prefix={folder_path} inference.input_pdb={input_pdb_path} "contigmap.contigs={args.contigs}" inference.num_designs={args.num_designs} "ppi.hotspot_res={args.hotspot_residues}"',
    ]

    if args.logs_dir:
        logs_dir = os.path.join(
            args.logs_dir, f"rfdiffusion/{folder_path.split('/')[-1]}"
        )
        print(f"\nSaving logs to: {logs_dir}\n")
        os.makedirs(logs_dir, exist_ok=True)

        with open(os.path.join(logs_dir, "stdout.txt"), "w") as stdout_file, open(
            os.path.join(logs_dir, "stderr.txt"), "w"
        ) as stderr_file:

            # Save outputs to log file
            subprocess.run(
                ["bash", "-c", "\n".join(commands)],
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
            )

    else:
        subprocess.run(["bash", "-c", "\n".join(commands)], text=True)

    print("\nFINISHED!\n")
