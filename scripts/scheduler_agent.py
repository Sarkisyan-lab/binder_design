import os
import glob
import json
import time
import argparse
import subprocess
import wandb
from dotenv import load_dotenv
import numpy as np


# load .env
load_dotenv()

CHAIN_B_DESC = ">ChainB|common_chain|MRP-2"
CHAIN_B_SEQ = "MFPLVRELVCGDKTELLEPGWKNRSSIPHVTKCGQHTDFSTIPTLFLVVFSFIVIYQLLHTRTAQLRCFSPISFRIILGCLLVVDLIATVIYDLYLYVSQSPNFEVVHFYGDLVQFGGICLALILTVACKNKGIITSGVITLYWLLVVVCGIPEFRFYLSGFIYNEYALEGIRATLYIIAFTFSALELFLCCFADVPSDMYKSESSCPEYTASFINRLTFQWFTGLAYLGNKKSLENEDLWDLNEIDKAENLIPSFMQNLKPRIDEYHQNIKKDPSAALPKNHPSFVIPIFKTYKYTLLAGFFYKLCFDMLQFLAPQLLKQLIGFIEDKNQPVWIGCSIVGIMFFSSFLQSMFLHQYYHSMFRLGMHVRSVLTSAVYSKALNLSNEARKGKTIGAIVNLMSVDIQKIQDMAPTIMLFWSAPLQIFLSIYFLWKFLGVAALAGLVVLILALPVNGLIAIQMRKCQTEQMKLKDERIKMMSEILNGMKVLKLYSWERSMENMVLKIRERELHILKKLSYFMAAIVFSWICAPFLASVISFVVYVYLDPENNVLTPEITFVALSLFDILRMPLAMVAMVYGEAVQCSVSNTRLKEFFAAEEMSPQTSISHGETDSAIEVENGLFSWSSDEDPTLREISFKIQKGQLVAIVGKVGSGKSSLLHALLGEMNKLSGSVQINGNIAYVPQQAWIQNMSLRNNILFNKPYDLENYEDVVKNCALKEDLANLPAGDRTEIGEKGINLSGGQKQRVSLARAVYQNPDIILLDDPLSAVDSHVGKHIFENVISSSTGCLASKTRVLVTHGLTYLKHCDQLIVLKEGTISELGTYQELLNNSGAFAEFLEEFLIEESKTRGRVASIGDGSGEVDEILRDLGQVKPGILKRLESHLSQESDKEDTSARAIEYSRDSSRRSVLLHSPRSQHEENEALLGAISEDVPAQENTQLIEKETVETGKVKFEVYIAYFQAISIPITLLFFFLYVGSSGLGILSNFYLAKLSDHAKSGNRTSSDAKMELGIYAVLGMGQSFVVLIASIILTIGVLRASRILHAGLLGNIMRSPMAFFDVTPIGRILNRIGKDIEAIDRTLPDVIRHMSMTIFNVVATLVVIMWATPWAGIAFAILSVIYFIVLRFYISTSRQLKRLESASRSPIYSHFQESIQGASSIRAFGVVDNFIKQSQQRVDDHLIAYYPSIVANRWLAVRLEMVGNLIVLSAAGAAVYFRDSPGLSAGLVGLSVSYALNITQTLNWAVRMTSELETNIVSVERIKEYTVTPTEGNNSRRLAAKSWPEKGEISIKNFSVRYRPGLDLVLHGISAHIAPSEKVGIVGRTGAGKSSLTLALFRIIEADGGSIEIDGINIANLQLEQLRSCLTIVPQDPVLFSGTMKMNLDPFSAYSDSQVWEALENAHLKPFVKSLQDGLEHKISEGGENLSVGQRQLICLARALLRKTKVLVLDEAAAAVDVETDSLIQKTIREQFKECTVLTIAHRLNTVMDSDRLLVLDKGRVAEFDSPKNLLANPDGIFYSMAKDANVV"
TABLE_NAME = "parafold_jobs_table"
TABLE_COLUMNS = [
    "file_name",
    "msa_status",
    "msa_device",
    "folding_status",
    "folding_device",
]
ARTIFACT_NAME = "parafold_jobs"
STATUS_COMPLETED = "completed"
STATUS_NOT_STARTED = "not_started"
STATUS_IN_PROGRESS = "in_progress"


class Scheduler:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.wandb = wandb.init(project=os.getenv("WANDB_PROJECT_NAME"))
        self.output_file_dir = os.getenv("OUTPUT_DIR_PATH")

    def add_common_chain_to_fasta(self, fasta_file_path: str):
        with open(fasta_file_path, "r") as f:
            fasta_content = f.readlines()

        filtered_fasta_content = [
            line.strip() for line in fasta_content if line != "\n"
        ]

        if len(filtered_fasta_content) == 2:
            filtered_fasta_content += [CHAIN_B_DESC, CHAIN_B_SEQ]
            with open(fasta_file_path, "w") as f:
                f.write("\n".join(filtered_fasta_content))

        elif len(filtered_fasta_content) == 4:
            print("Common Chain already exists")

        else:
            print("Fasta file has more than 2 sequences")

    def get_fasta_file_paths(self):
        # check both .fasta and .fa files
        prot_mpnn_out_files = os.path.join(
            os.getenv("OUTPUT_DIR_PATH"), "ligand_mpnn/processed"
        )
        fasta_files = glob.glob(f"{prot_mpnn_out_files}/*.fasta") + glob.glob(
            f"{prot_mpnn_out_files}/*.fa"
        )
        return fasta_files

    @staticmethod
    def obtain_parafold_out_path(fasta_file_dir: str):
        if fasta_file_dir.endswith(".fa"):
            trimmed_fasta_dirname = fasta_file_dir.replace(".fa", "")
        elif fasta_file_dir.endswith(".fasta"):
            trimmed_fasta_dirname = fasta_file_dir.replace(".fasta", "")
        else:
            print("Invalid fasta file")
        # TODO: Remove outmulti once sto is done
        out_path = trimmed_fasta_dirname.replace("ligand_mpnn/processed/", "parafold/")
        return out_path

    def submit_job_lilibet(self, fasta_file_dir, job_type="folding"):
        final_output_dir = self.obtain_parafold_out_path(fasta_file_dir)
        common_chain_msa_path = os.path.join(
            os.getenv("OUTPUT_DIR_PATH"), "mrp_2_msas/msas"
        )

        if job_type == "folding":
            af_multimer_command = f"bash run_alphafold.sh -d {os.getenv('AFDB_PATH')} -o {os.path.dirname(final_output_dir)} -m model_1_multimer_v3 -p multimer -i {fasta_file_dir} -t 2020-12-01 -Y common_chain -Z {common_chain_msa_path} -u {self.args.gpu_id} -r best -N 1"
        else:
            af_multimer_command = f"bash run_alphafold.sh -d {os.getenv('AFDB_PATH')} -o {final_output_dir} -m model_1_multimer_v3 -p multimer -i {fasta_file_dir} -t 2020-12-01"

        commands = [
            f"source {os.getenv('CONDA_PATH')}",
            f"cd {os.getenv('PARAFOLD_REPO_PATH')}",
            f"conda activate {os.getenv('PARAFOLD_ENV')}",
            af_multimer_command,
        ]

        try:
            logs_dir = os.path.join(
                args.logs_dir, f"parafold/{final_output_dir.split('/')[-1]}_{job_type}"
            )
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
            return True

        except Exception as e:
            print(f"Error submitting job: {e}")
            return False

    def retrieve_wandb_jobs(self, shuffle=True):
        # Retrieve the created table using the artifact you created.
        artifact = self.wandb.use_artifact(f"{ARTIFACT_NAME}:latest")
        table = artifact.get(TABLE_NAME).get_dataframe().values
        dict_table = [
            {
                "file_name": row[0],
                "msa_status": row[1],
                "msa_device": row[2],
                "folding_status": row[3],
                "folding_device": row[4],
            }
            for row in table
        ]
        if shuffle:
            np.random.shuffle(dict_table)
        return dict_table

    def retrieve_local_jobs(self):
        fasta_files = self.get_fasta_file_paths()
        table = []
        for fasta_file in fasta_files:
            out_file_dir = self.obtain_parafold_out_path(fasta_file)

            # check msa_gen status
            msa_status = STATUS_NOT_STARTED
            msa_device = "none"

            chain_id_map_path = os.path.join(out_file_dir, "msas/chain_id_map.json")
            if os.path.isfile(chain_id_map_path):
                with open(chain_id_map_path, "r") as f:
                    chain_id_map = json.load(f)
                if len(chain_id_map) == 2:
                    msa_status = STATUS_COMPLETED
                    msa_device = self.args.device

            # check GPU status
            ranking_debug_file_path = os.path.join(out_file_dir, "ranking_debug.json")
            if os.path.isfile(ranking_debug_file_path):
                folding_status = STATUS_COMPLETED
                folding_device = self.args.device
            else:
                folding_status = STATUS_NOT_STARTED
                folding_device = "none"

            # file_name, msa_status, msa_device, folding_status, folding_device
            table.append(
                {
                    "file_name": out_file_dir.split("/")[-1],
                    "msa_status": msa_status,
                    "msa_device": msa_device,
                    "folding_status": folding_status,
                    "folding_device": folding_device,
                }
            )
        return table

    def upload_table_wandb(self, table, table_name=TABLE_NAME):
        artifact = wandb.Artifact(ARTIFACT_NAME, type="dataset")
        table_to_upload = wandb.Table(
            columns=TABLE_COLUMNS,
            data=[list(item.values()) for item in table],
        )
        artifact.add(table_to_upload, table_name)
        wandb.log_artifact(artifact)

    def sync_table(self, force=False):
        # Retrieve the table from wandb
        wandb_table = self.retrieve_wandb_jobs()
        print(f"Total Wandb jobs: {len(wandb_table)}")

        # Generate the table from the local directory
        local_table = self.retrieve_local_jobs()
        print(f"Total Local jobs: {len(local_table)}")

        wandb_job_names = [job["file_name"] for job in wandb_table]
        local_job_names = [job["file_name"] for job in local_table]

        for job in local_job_names:
            # If job found in wandb_table, check if everything is updated
            local_job_idx = local_job_names.index(job)
            if job in wandb_job_names:
                wandb_job_idx = wandb_job_names.index(job)
                for j_type in ["msa", "folding"]:
                    if (
                        wandb_table[wandb_job_idx][f"{j_type}_status"]
                        != STATUS_COMPLETED
                        and local_table[local_job_idx][f"{j_type}_status"]
                        == STATUS_COMPLETED
                    ):
                        print(f"Making {j_type} {job} as complete in wandb table...")
                        wandb_table[wandb_job_idx][f"{j_type}_status"] = local_table[
                            local_job_idx
                        ][f"{j_type}_status"]
                        wandb_table[wandb_job_idx][
                            f"{j_type}_device"
                        ] = self.args.device

                    if (
                        force == True
                        and wandb_table[wandb_job_idx][f"{j_type}_status"]
                        == STATUS_IN_PROGRESS
                        and local_table[local_job_idx][f"{j_type}_status"]
                        == STATUS_NOT_STARTED
                    ):
                        print(f"Making {job} as not_started in wandb table...")
                        wandb_table[wandb_job_idx][
                            f"{j_type}_status"
                        ] = STATUS_NOT_STARTED
                        wandb_table[wandb_job_idx][f"{j_type}_device"] = "none"

            # If job not found in wandb_table, add it
            else:
                print(f"Adding new job: {job} to wandb table...")
                wandb_table.append(local_table[local_job_idx])
        print("Finished Syncing...")
        self.upload_table_wandb(wandb_table)

    def update_table(self, table: np.array, file_name: str, status: str, job_type: str):
        print(f"Updating {file_name} {job_type} status to {status}...")
        # update row in the table
        for row_idx, row in enumerate(table):
            if row["file_name"] == file_name:
                table[row_idx][f"{job_type}_status"] = status
                table[row_idx][f"{job_type}_device"] = self.args.device
        self.upload_table_wandb(table)

    def check_if_msas_exist(self, out_file_dir: str):
        # check msa_gen status
        feat_pkl_file = os.path.join(out_file_dir, "features.pkl")
        uniprot_sto_file = os.path.join(out_file_dir, "msas/A/uniprot_hits.sto")
        chain_id_map_path = os.path.join(out_file_dir, "msas/chain_id_map.json")
        if os.path.isfile(chain_id_map_path):
            with open(chain_id_map_path, "r") as f:
                chain_id_map = json.load(f)
            if (
                len(chain_id_map) == 2
                and os.path.isfile(feat_pkl_file)
                and os.path.isfile(uniprot_sto_file)
            ):
                return True
        return False

    def sync_msas(self, msa_folder_path: str, from_device: str, to_device: str):
        if from_device == to_device:
            print("From and to device are same. No need to sync")

        if to_device == "lilibet":
            if from_device == "hpc":
                print("Syncing from hpc to lilibet")

                if os.path.exists(msa_folder_path):
                    print(f"Deleting existing output directory: {msa_folder_path}")
                    subprocess.run(f"rm -r {msa_folder_path}", shell=True)

                hpc_path = msa_folder_path.replace(
                    os.getenv("LILIBET_BASE"), os.getenv("HPC_BASE")
                )

                subprocess.run(
                    f"sshpass -p {os.getenv('PASSWORD_HPC')} scp -r {os.getenv('HPC_SERVER')}:{hpc_path} {os.path.dirname(msa_folder_path)}",
                    shell=True,
                )

    def main(self, job_type: str, max_jobs=1, sleep_duration_sec=60):
        # OUTPUT_DIR_PATH = os.getenv("OUTPUT_DIR_PATH")
        fasta_files_local = self.get_fasta_file_paths()
        fasta_file_names = [
            file.split("/")[-1].replace(".fasta", "").replace(".fa", "")
            for file in fasta_files_local
        ]

        job_idx = 0
        sleep_counter = 0
        while job_idx <= max_jobs:
            print(f"\nJob Idx: {job_idx}: Checking for jobs to submit...\n")

            # sync and retrieve current jobs
            self.sync_table()
            curr_wandb_jobs = self.retrieve_wandb_jobs()

            job_found = False
            for job in curr_wandb_jobs:
                # folding job
                if (
                    job_type == "folding"
                    and job["msa_status"]
                    == STATUS_COMPLETED  # completed msa generation
                    and job["folding_status"]
                    == STATUS_NOT_STARTED  # not started folding
                ) or (
                    job_type == "msa"
                    and job["msa_status"]
                    == STATUS_NOT_STARTED  # incomplete msa generation
                    and job["folding_status"]
                    == STATUS_NOT_STARTED  # not started folding
                ):
                    job_found = True
                    self.update_table(
                        table=curr_wandb_jobs,
                        file_name=job["file_name"],
                        status=STATUS_IN_PROGRESS,
                        job_type=job_type,
                    )
                    full_inp_path = fasta_files_local[
                        fasta_file_names.index(job["file_name"])
                    ]
                    if not self.check_if_msas_exist(
                        self.obtain_parafold_out_path(full_inp_path),
                    ):
                        self.sync_msas(
                            msa_folder_path=self.obtain_parafold_out_path(
                                full_inp_path
                            ),
                            from_device=job["msa_device"],
                            to_device=self.args.device,
                        )
                    # TODO: Submit Job
                    print(f"\nJob Idx: {job_idx}: Submitting {job['file_name']} \n")

                    self.add_common_chain_to_fasta(full_inp_path)
                    self.submit_job_lilibet(
                        fasta_file_dir=full_inp_path,
                        job_type=job_type,
                    )
                    print(f"\nFinished Job: {job['file_name']} \n")
                    self.update_table(
                        table=curr_wandb_jobs,
                        file_name=job["file_name"],
                        status=STATUS_COMPLETED,
                        job_type=job_type,
                    )
                    self.sync_table()
                    job_idx += 1
                    break
            if sleep_counter == 5:
                self.sync_table()
                return None
            if job_found == False:
                sleep_counter += 1
                print(
                    f"Found no relevant jobs... Sleeping for {sleep_duration_sec} seconds..."
                )
                time.sleep(sleep_duration_sec)
        print("\nAll Jobs Finished!")


if __name__ == "__main__":
    # Parse the arguments
    # Create the parser
    parser = argparse.ArgumentParser()

    # Add an argument
    parser.add_argument(
        "--device",
        type=str,
        help="Option between hpc / lilibet / jex",
        default="lilibet",
    )
    parser.add_argument(
        "--logs_dir",
        type=str,
        help="Path to logs_dir.",
        required=True,
    )

    parser.add_argument(
        "--max_jobs",
        type=int,
        help="Number of jobs to submit",
        required=True,
    )

    parser.add_argument(
        "--job_type",
        type=str,
        help="Type of job to submit: msa / folding",
        required=True,
    )
    parser.add_argument(
        "--sleep_duration_sec",
        type=int,
        help="Duration to sleep if job not found",
        default=60,
    )
    parser.add_argument(
        "--gpu_id",
        type=int,
        help="GPU ID",
        default=0,
    )

    args = parser.parse_args()
    scheduler = Scheduler(args)

    try:
        scheduler.main(
            job_type=args.job_type,
            max_jobs=args.max_jobs,
            sleep_duration_sec=args.sleep_duration_sec,
        )

    except KeyboardInterrupt:
        print("Interrupted by user")
        scheduler.sync_table()
        print("Exiting...")
