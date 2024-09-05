import os
import subprocess
import firebase_admin
from firebase_admin import db
import glob
import pandas as pd
from argparse import Namespace
from click import command, option
import json

# constants
STATUS_UNASSIGNED = "unassigned"
STATUS_NOT_STARTED = "not_started"
STATUS_COMPLETED = "completed"
DEVICE_LILIBET = "lilibet"
DEVICE_HPC_CX3 = "hpc_cx3"
DEVICE_HPC_BASE = "hpc_base"
DEVICE_JEX = "jex"
JOB_TYPE_MSA = "msa"
JOB_TYPE_FOLDING = "folding"


class TaskScheduler:
    def __init__(self, config: dict, device=DEVICE_LILIBET):
        self.config = Namespace(**config)
        self.device = device
        cred_obj = firebase_admin.credentials.Certificate(
            os.path.join(
                os.path.dirname(os.path.dirname(__file__)), "src", "login_key.json"
            )
        )
        self.default_app = firebase_admin.initialize_app(
            cred_obj, {"databaseURL": self.config.firebase_db_url}
        )
        self.ref = db.reference("/")

    def empty_db(self):
        objs = self.ref.get()
        for key, _ in objs.items():
            self.ref.child(key).set({})

    def fasta_file_to_fasta_str(self, fasta_file_path: str) -> str:
        with open(fasta_file_path, "r") as f:
            fasta_lines = f.readlines()
        fasta_content = fasta_lines[0] + "".join(fasta_lines[1:]).replace("\n", "")
        return fasta_content

    def save_fasta_from_fasta_str(self, fasta_str: str, file_name: str) -> str:
        if not os.path.exists(os.path.dirname(file_name)):
            os.makedirs(os.path.dirname(file_name), exist_ok=True)
        with open(file_name, "w") as f:
            for line in fasta_str.split():
                f.write(line + "\n")

    def inflate_tasks_from_fasta_dir(self, fasta_dir):
        for fasta_file in glob.glob(os.path.join(fasta_dir, "*.fasta")) + glob.glob(
            os.path.join(fasta_dir, "*.fa")
        ):
            if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
                fasta_content = self.fasta_file_to_fasta_str(fasta_file)
                print("Adding task for", os.path.basename(fasta_file))
                self.ref.push(
                    {
                        "file_name": os.path.basename(fasta_file).split(".")[0],
                        "seq": fasta_content,
                        "msa_status": STATUS_NOT_STARTED,
                        "folding_status": STATUS_NOT_STARTED,
                        "msa_device": STATUS_UNASSIGNED,
                        "folding_device": STATUS_UNASSIGNED,
                    }  # type: ignore
                )

    def get_tasks_by_filter(self, filter_criteria: dict) -> list:
        objs = self.ref.get()
        filtered_tasks = []
        for key, value in objs.items():
            if all(value[k] == v for k, v in filter_criteria.items()):
                filtered_tasks.append(value)
        return filtered_tasks

    def get_local_tasks_status(self, fasta_dir: str) -> list:
        local_task_status = []
        for fasta_file in glob.glob(os.path.join(fasta_dir, "*.fasta")) + glob.glob(
            os.path.join(fasta_dir, "*.fa")
        ):
            if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
                file_name = os.path.basename(fasta_file).split(".")[0]
                task_status = {
                    "file_name": file_name,
                    "seq": self.fasta_file_to_fasta_str(fasta_file),
                    "msa_status": STATUS_NOT_STARTED,
                    "folding_status": STATUS_NOT_STARTED,
                    "msa_device": STATUS_UNASSIGNED,
                    "folding_device": STATUS_UNASSIGNED,
                }

                # If folding / msa is completed
                if self.device == DEVICE_LILIBET:
                    base_output_dir = self.config.lilibet_output_dir
                elif self.device == DEVICE_HPC_CX3 or self.device == DEVICE_HPC_BASE:
                    base_output_dir = self.config.hpc_output_dir
                elif self.device == DEVICE_JEX:
                    base_output_dir = self.config.jex_output_dir
                else:
                    raise ValueError(f"Unknown device: {self.device}")

                msas_dir = os.path.join(base_output_dir, "msas")
                predictions_dir = os.path.join(
                    base_output_dir, "predictions", file_name
                )
                # Check if MSAs exist
                if os.path.exists(msas_dir):
                    msas_folder_contents = os.listdir(msas_dir)
                    for file in msas_folder_contents:
                        if file_name in file and file.endswith(".a3m"):
                            task_status["msa_status"] = STATUS_COMPLETED
                            task_status["msa_device"] = self.device
                            break

                # Check if predictions exist
                if os.path.exists(predictions_dir):
                    predictions_folder_contents = os.listdir(predictions_dir)
                    for file in predictions_folder_contents:
                        if "rank_001_alphafold2_multimer_v3" in file:
                            task_status["folding_status"] = STATUS_COMPLETED
                            task_status["folding_device"] = self.device
                            break

                local_task_status.append(task_status)
        return local_task_status

    def find_db_id_for_file_name(self, file_name: str) -> dict:
        db_tasks = self.ref.get("/")[0]
        for db_task_key, db_task_value in db_tasks.items():
            if db_task_value["file_name"] == file_name:
                return db_task_key
        return None

    def update_db_with_completed_local_tasks(self, local_tasks_status: list):
        # Assuming that the db has the entire list of tasks already
        db_tasks = self.ref.get()
        for local_task in local_tasks_status:
            # Find the task in the db for a given file name
            resp_db_task_key, resp_db_task_value = self.find_db_task_for_file_name(
                local_task["file_name"], db_tasks
            )
            for job_type in ["msa", "folding"]:
                if (
                    local_task[f"{job_type}_status"] == STATUS_COMPLETED
                    and resp_db_task_value[f"{job_type}_status"] != STATUS_COMPLETED
                ):
                    print(
                        f"Updating {job_type} status for {local_task['file_name']} to {STATUS_COMPLETED}"
                    )
                    # Update the db with the completed task
                    self.ref.child(resp_db_task_key).child(f"{job_type}_status").set(
                        STATUS_COMPLETED
                    )
                    self.ref.child(resp_db_task_key).child(f"{job_type}_device").set(
                        self.device
                    )

    def submit_job_msa_local_hpc_cx3_local(self, job_details: dict):
        fasta_file_name = job_details["file_name"]
        # Create directories locally
        os.makedirs(
            os.path.join(self.config.hpc_output_dir, fasta_file_name, "input"),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(self.config.hpc_output_dir, fasta_file_name, "output"),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(self.config.hpc_output_dir, fasta_file_name, "logs"),
            exist_ok=True,
        )

        # Save the fasta file locally
        self.save_fasta_from_fasta_str(
            fasta_str=job_details["seq"],
            file_name=os.path.join(
                self.config.hpc_output_dir,
                fasta_file_name,
                "input",
                f"{fasta_file_name}.fasta",
            ),
        )

        commands = [
            "#!/bin/bash",
            f"#PBS -l select=1:ncpus={self.config.hpc_msa_job_num_cpus_local}:mem={self.config.hpc_msa_job_mem_gb_local}gb",
            f"#PBS -l walltime={self.config.hpc_msa_job_time_local}",
            f"#PBS -N {fasta_file_name}",
            f"#PBS -e {self.config.hpc_output_dir}/{fasta_file_name}/logs/$PBS_JOBID_msa_local.err",
            f"#PBS -o {self.config.hpc_output_dir}/{fasta_file_name}/logs/$PBS_JOBID_msa_local.out",
            f"exec > >(tee -a {self.config.hpc_output_dir}/{fasta_file_name}/logs/$PBS_JOBNAME_msa_local.out)",
            "cd $PBS_O_WORKDIR",
            'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
            f"conda activate {self.config.hpc_colabfold_conda_env}",
            f"colabfold_search {self.config.hpc_output_dir}/{fasta_file_name}/input/{fasta_file_name}.fasta {self.config.hpc_dbs_loc} {self.config.hpc_output_dir}/{fasta_file_name}/output/msas",
            f"scp -P 10002 -r {self.config.hpc_output_dir}/{fasta_file_name}/output {self.config.lilibet_host}:{self.config.lilibet_output_dir}/{fasta_file_name}/",
        ]

        job_file_path = os.path.join(
            self.config.hpc_output_dir,
            fasta_file_name,
            "input",
            f"{fasta_file_name}_msa_local.pbs",
        )
        with open(job_file_path, "w") as f:
            for command in commands:
                f.write(command + "\n")

        # submit_job
        # try:
        #     subprocess.run(f"qsub {job_file_path}", shell=True)
        # except subprocess.CalledProcessError as e:
        #     raise Exception(f"Error submitting job on HPC: {e}")

    def submit_job_msa_local_hpc_from_fasta_dir(
        self, fasta_dir: str, use_templates=True, copy_to_lilibet=True
    ):
        """
        Use this function to submit a batch MSA generation job to the HPC
        which uses the high throughput mmseqs2 implementation of colabfold
        which is only setup on cx3 cluster

        Creates an intermediate directory structure
        Args:
            fasta_dir (str): Directory containing fasta files to run MSA generation on.
                This assumes every fasta has one protein (could be multimer).
            use_templates (bool): Whether to use templates for MSA generation
            copy_to_lilibet (bool): Whether to copy the results to lilibet
        """
        # Check how many files are in the fasta_dir
        fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta")) + glob.glob(
            os.path.join(fasta_dir, "*.fa")
        )
        if len(fasta_files) == 0:
            raise ValueError(f"No fasta files found in {fasta_dir}")
        else:
            print(f"Found {len(fasta_files)} fasta files in {fasta_dir}.")

        # Convert the fasta files into a csv file for colabfold_search
        fasta_data = []
        for fasta_file in fasta_files:
            with open(fasta_file, "r") as f:
                lines = f.readlines()
                description = lines[0].strip().lstrip(">")
                sequence = "".join(line.strip() for line in lines[1:])
                fasta_data.append({"id": description, "seq": sequence})

        df = pd.DataFrame(fasta_data)
        csv_output_path = os.path.join(fasta_dir, "fasta_contents.csv")
        df.to_csv(csv_output_path, index=False)

        # mkdir logs and msas
        logs_dir = os.path.join(self.config.hpc_output_dir, "logs", "msas")
        msas_dir = os.path.join(self.config.hpc_output_dir, "msas")
        os.makedirs(logs_dir, exist_ok=True)
        os.makedirs(msas_dir, exist_ok=True)

        # Check if templates are used
        if use_templates:
            use_templates_idx = "1"
        else:
            use_templates_idx = "0"

        # Run colabfold_search
        commands = [
            "#!/bin/bash",
            f"#PBS -l select=1:ncpus={self.config.hpc_clf_search_num_cpus}:mem={self.config.hpc_clf_search_mem_gb}gb",
            f"#PBS -l walltime={self.config.hpc_clf_search_time}",
            f"#PBS -e {logs_dir}/",
            f"#PBS -o {logs_dir}/",
            f"exec > >(tee -a {logs_dir}/"
            + "$colabfold_search_batch_${PBS_JOBID}.out)",
            "cd $PBS_O_WORKDIR",
            'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
            f"conda activate {self.config.hpc_colabfold_conda_env}",
            f"colabfold_search {csv_output_path} {self.config.hpc_dbs_loc} {msas_dir} --threads {self.config.hpc_clf_search_num_cpus-4} --use-templates {use_templates_idx} --db-load-mode 2 --use-env 1 --db2 pdb100_230517",
        ]

        if copy_to_lilibet:
            commands.append(
                f"scp -P 10002 -r {msas_dir} {self.config.lilibet_host}:{self.config.lilibet_output_dir}"
            )

        job_file_path = os.path.join(logs_dir, f"msa_search.pbs")
        with open(job_file_path, "w") as f:
            for command in commands:
                f.write(command + "\n")

        # submit_job
        try:
            subprocess.run(f"qsub {job_file_path}", shell=True)
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error submitting job on HPC: {e}")

    def job_housekeeping(self, job_details: dict, device: str):
        if device == DEVICE_HPC_CX3 or device == DEVICE_HPC_BASE:
            base_path = self.config.hpc_output_dir
        elif device == DEVICE_LILIBET:
            base_path = self.config.lilibet_output_dir
        elif device == DEVICE_JEX:
            base_path = self.config.jex_output_dir
        else:
            raise ValueError(f"Unknown device: {device}")

        fasta_file_name = job_details["file_name"]

        # Check if fasta exists locally
        fasta_file_full_path = os.path.join(
            base_path,
            "input",
            f"{fasta_file_name}.fasta",
        )
        if not os.path.exists(fasta_file_full_path):
            raise FileNotFoundError(f"Fasta file not found: {fasta_file_full_path}")

        # Check if MSAs exist
        MSA_FOUND = False
        msa_output_dir = os.path.join(base_path, "msas")
        if not os.path.exists(msa_output_dir) or len(os.listdir(msa_output_dir)) == 0:
            MSA_FOUND = False
        else:
            msas_list = os.listdir(msa_output_dir)
            for msa_file in msas_list:
                if msa_file.endswith(".a3m") and fasta_file_name in msa_file:
                    MSA_FOUND = True
                    break
        if not MSA_FOUND:
            raise FileNotFoundError(f"MSAs do not exist for {fasta_file_name}!")

        # Make logs directory if doesn't exist
        logs_dir = os.path.join(base_path, "logs", fasta_file_name)
        os.makedirs(logs_dir, exist_ok=True)

        # Make predictions directory
        predictions_dir = os.path.join(base_path, "predictions", fasta_file_name)
        os.makedirs(predictions_dir, exist_ok=True)
        return logs_dir, predictions_dir, msa_output_dir, fasta_file_full_path

    def submit_job_folding_hpc_cx3(
        self, job_details: dict, copy_to_lilibet: bool = True
    ):
        fasta_file_name = job_details["file_name"]
        # Check if fasta exists locally
        logs_dir, predictions_dir, msa_output_dir, fasta_file_full_path = (
            self.job_housekeeping(job_details, device=DEVICE_HPC_CX3)
        )
        drop_out_str = " --use-dropout" if self.config.colabfold_dropout else ""
        # Folding commands
        commands = [
            "#!/bin/bash",
            f"#PBS -l select=1:ncpus={self.config.hpc_folding_job_num_cpus}:mem={self.config.hpc_folding_job_mem_gb}gb:ngpus=1",
            f"#PBS -l walltime={self.config.hpc_folding_job_time}",
            f"#PBS -N {fasta_file_name}",
            f"#PBS -e {self.config.hpc_output_dir}/logs/{fasta_file_name}/",
            f"#PBS -o {self.config.hpc_output_dir}/logs/{fasta_file_name}/",
            f"exec > >(tee -a {self.config.hpc_output_dir}/logs/{fasta_file_name}/"
            + "${PBS_JOBNAME}_${PBS_JOBID}_folding.out)",
            "cd $PBS_O_WORKDIR",
            'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
            f"conda activate {self.config.hpc_colabfold_conda_env}",
            f"colabfold_batch {fasta_file_full_path} {predictions_dir} --templates --amber --overwrite-existing-results --num-models {self.config.colabfold_num_models} --num-recycle {self.config.colabfold_num_recycle} {drop_out_str}",
        ]
        if copy_to_lilibet:
            commands.append(
                f"scp -P 10002 -r {predictions_dir} {self.config.lilibet_host}:{self.config.lilibet_output_dir}/predictions/",
            )

        job_file_path = os.path.join(logs_dir, f"{fasta_file_name}_folding.pbs")
        with open(job_file_path, "w") as f:
            for command in commands:
                f.write(command + "\n")

        # submit_job
        # try:
        #     subprocess.run(f"qsub {job_file_path}", shell=True)
        # except subprocess.CalledProcessError as e:
        #     raise Exception(f"Error submitting job on HPC: {e}")

    def submit_job_folding_lilibet(self, job_details_list: list):
        print(
            f"Running Folding for {len(job_details_list)} sequences on lilibet parallely"
        )
        fold_cmd = []
        # scp_cmd = []
        for job_details in job_details_list:
            fasta_file_name = job_details["file_name"]
            logs_dir, predictions_dir, msa_output_dir, fasta_file_full_path = (
                self.job_housekeeping(job_details, device=DEVICE_LILIBET)
            )
            msa_file_path = os.path.join(msa_output_dir, f"{fasta_file_name}.a3m")
            template_file_path = os.path.join(
                msa_output_dir, f"{fasta_file_name}_pdb100_230517.m8"
            )
            fold_cmd.append(
                f"colabfold_batch {msa_file_path} {predictions_dir} --num-recycle {self.config.colabfold_num_recycle} --num-models {self.config.colabfold_num_models}"
            )
            # scp_cmd.append(
            #     f"scp -r -P 10002 {predictions_dir} {self.config.lilibet_host}:{self.config.lilibet_output_dir}/predictions/"
            # )
        fold_cmd_str = " & ".join(fold_cmd)
        # scp_cmd_str = " & ".join(scp_cmd)

        # define a common set of commands
        commands = [
            f"source ~/anaconda3/etc/profile.d/conda.sh",
            f"conda activate {self.config.lilibet_colabfold_conda_env}",
            fold_cmd_str,
            "wait",
            # scp_cmd_str,
        ]

        # return commands
        with open(os.path.join(logs_dir, "folding.txt"), "w") as stdout_file, open(
            os.path.join(logs_dir, "folding.err"),
            "w",
        ) as stderr_file:
            print("Submitting folding job on lilibet")
            subprocess.run(
                ["bash", "-c", "\n".join(commands)],
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
            )
            print("Finished folding on lilibet")
        return True

    def run_job(self, job_details_list: list):
        # Set the job status to queued
        for job_details in job_details_list:
            print(f"Setting {job_details['file_name']} to queued")
            self.ref.child(
                self.find_db_id_for_file_name(job_details["file_name"])
            ).child(f"folding_status").set("queued")
            self.ref.child(
                self.find_db_id_for_file_name(job_details["file_name"])
            ).child(f"folding_device").set(self.device)

        # Submit the job
        if self.device == DEVICE_LILIBET:
            job_status = self.submit_job_folding_lilibet(job_details_list)
        else:
            raise ValueError(
                f"Folding on other devices not implemented yet: {self.device}"
            )

        # Update the job status if the job finished successfully
        for job_details in job_details_list:
            if job_status:
                print(f"Setting {job_details['file_name']} to completed")
                self.ref.child(
                    self.find_db_id_for_file_name(job_details["file_name"])
                ).child(f"folding_status").set(STATUS_COMPLETED)
            else:
                print(f"Setting {job_details['file_name']} to not started")
                self.ref.child(
                    self.find_db_id_for_file_name(job_details["file_name"])
                ).child(f"folding_status").set(STATUS_NOT_STARTED)
                self.ref.child(
                    self.find_db_id_for_file_name(job_details["file_name"])
                ).child(f"folding_device").set(STATUS_UNASSIGNED)


@command()
@option("--config_file_path", help="Path to the config file")
@option("--device_name", default="lilibet", help="Device to run the scheduler on")
@option("--max_jobs", default=100, help="Max number of jobs to run")
@option("--num_jobs_per_gpu", default=2, help="Number of jobs to run per GPU")
def run_scheduler_lilibet(
    config_file_path: str,
    device_name: str,
    max_jobs: int,
    num_jobs_per_gpu: int = 2,
):
    # load the config file
    with open(config_file_path, "r") as f:
        config = json.load(f)

    # Initialize the TaskScheduler
    ts = TaskScheduler(config, device_name)

    print("Starting Scheduler...")
    for job_idx in range(1, max_jobs + 1):

        db_tasks = ts.get_tasks_by_filter(
            {
                "folding_status": STATUS_NOT_STARTED,
                "msa_status": STATUS_COMPLETED,
            }
        )
        if len(db_tasks) == 0:
            print("No tasks to run")
            break

        # Run the task
        curr_jobs = db_tasks[:num_jobs_per_gpu]
        print(f"{job_idx}. Running jobs: {[a['file_name'] for a in curr_jobs]}")
        ts.run_job(curr_jobs)
        print(f"Finished running all jobs for job idx: {job_idx}")


if __name__ == "__main__":
    run_scheduler_lilibet()
