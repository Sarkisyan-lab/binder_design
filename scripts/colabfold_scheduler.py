import os
import subprocess
import firebase_admin
from firebase_admin import db
import glob
import pandas as pd
from argparse import Namespace
from click import command, option

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

    def find_db_task_for_file_name(self, file_name: str, db_tasks: dict) -> dict:
        for db_task_key, db_task_value in db_tasks.items():
            if db_task_value["file_name"] == file_name:
                return db_task_key, db_task_value
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

    def submit_job_folding_hpc_cx3(
        self, job_details: dict, copy_to_lilibet: bool = True, hpc_type: str = "cx3"
    ):
        fasta_file_name = job_details["file_name"]

        # Save the fasta file locally
        fasta_file_full_path = os.path.join(
            self.config.hpc_output_dir,
            "input",
            f"{fasta_file_name}.fasta",
        )
        if not os.path.exists(fasta_file_full_path):
            raise FileNotFoundError(f"Fasta file not found: {fasta_file_full_path}")

        drop_out_str = " --use-dropout" if self.config.colabfold_dropout else ""

        # Check if MSAs exist
        msa_output_dir = os.path.join(self.config.hpc_output_dir, "msas")
        if not os.path.exists(msa_output_dir) or len(os.listdir(msa_output_dir)) == 0:
            raise FileNotFoundError(f"MSAs do not exist for {fasta_file_name}!")

        # Make logs directory if doesn't exist
        logs_dir = os.path.join(self.config.hpc_output_dir, "logs", fasta_file_name)
        os.makedirs(logs_dir, exist_ok=True)

        # Make predictions directory
        predictions_dir = os.path.join(
            self.config.hpc_output_dir, "predictions", fasta_file_name
        )
        os.makedirs(predictions_dir, exist_ok=True)

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

    def submit_job_msa_server_lilibet(self, job_details: dict):
        fasta_file_name = job_details["file_name"]
        # Create directories locally
        os.makedirs(
            os.path.join(self.config.lilibet_output_dir, fasta_file_name, "input"),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(self.config.lilibet_output_dir, fasta_file_name, "output"),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(self.config.lilibet_output_dir, fasta_file_name, "logs"),
            exist_ok=True,
        )

        # Save the fasta file locally
        fasta_file_full_path = os.path.join(
            self.config.lilibet_output_dir,
            fasta_file_name,
            "input",
            f"{fasta_file_name}.fasta",
        )
        self.save_fasta_from_fasta_str(
            fasta_str=job_details["seq"],
            file_name=fasta_file_full_path,
        )
        commands = [
            f"source ~/anaconda3/etc/profile.d/conda.sh",
            f"conda activate {self.config.lilibet_colabfold_conda_env}",
            "export JAX_PLATFORMS=cpu",
            f"colabfold_batch --templates --msa-only --overwrite-existing-results {fasta_file_full_path} {self.config.lilibet_output_dir}/{fasta_file_name}/output/",
        ]

        logs_dir = os.path.join(self.config.lilibet_output_dir, fasta_file_name, "logs")
        with open(
            os.path.join(logs_dir, "msa_gen_server.txt"), "w"
        ) as stdout_file, open(
            os.path.join(logs_dir, "msa_gen_server.err"),
            "w",
        ) as stderr_file:
            print(f"Submitting MSA generation job on lilibet for: {fasta_file_name}")
            subprocess.run(
                ["bash", "-c", "\n".join(commands)],
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
            )
            print(f"Finished MSA generation job on lilibet for: {fasta_file_name}")

    def submit_job_folding_lilibet(self, job_details: dict):
        fasta_file_name = job_details["file_name"]
        # Create directories locally
        os.makedirs(
            os.path.join(self.config.lilibet_output_dir, fasta_file_name, "input"),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(self.config.lilibet_output_dir, fasta_file_name, "output"),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(self.config.lilibet_output_dir, fasta_file_name, "logs"),
            exist_ok=True,
        )

        # Save the fasta file locally
        fasta_file_full_path = os.path.join(
            self.config.lilibet_output_dir,
            fasta_file_name,
            "input",
            f"{fasta_file_name}.fasta",
        )
        self.save_fasta_from_fasta_str(
            fasta_str=job_details["seq"],
            file_name=fasta_file_full_path,
        )
        drop_out_str = " --use-dropout " if self.config.colabfold_dropout else " "
        commands = [
            f"source ~/anaconda3/etc/profile.d/conda.sh",
            f"conda activate {self.config.lilibet_colabfold_conda_env}",
            f"colabfold_batch --templates --amber --overwrite-existing-results --num-models {self.config.colabfold_num_models}{drop_out_str}--num-recycle {self.config.colabfold_num_recycle} {fasta_file_full_path} {self.config.lilibet_output_dir}/{fasta_file_name}/output/",
        ]

        logs_dir = os.path.join(self.config.lilibet_output_dir, fasta_file_name, "logs")
        with open(
            os.path.join(logs_dir, "folding_gen_server.txt"), "w"
        ) as stdout_file, open(
            os.path.join(logs_dir, "folding_gen_server.err"),
            "w",
        ) as stderr_file:
            print(f"Submitting folding job on lilibet for: {fasta_file_name}")
            subprocess.run(
                ["bash", "-c", "\n".join(commands)],
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
            )
            print(f"Finished folding job on lilibet for: {fasta_file_name}")

    def run_job(self, job_details: dict, job_type: str):
        # Set the job status to queued
        self.ref.child(self.find_db_task_for_file_name(job_details["file_name"])).child(
            f"{job_type}_status"
        ).set("queued")
        self.ref.child(self.find_db_task_for_file_name(job_details["file_name"])).child(
            f"{job_type}_device"
        ).set(self.device)

        # Submit the job
        if job_type == JOB_TYPE_MSA:
            if self.device == DEVICE_LILIBET:
                job_status = self.submit_job_msa_server_lilibet(job_details)
            elif self.device == DEVICE_HPC_CX3 or self.device == DEVICE_HPC_BASE:
                job_status = self.submit_job_msa_local_hpc_cx3_local(job_details)
            else:
                raise ValueError(
                    f"MSAs on other devices not implemented yet: {self.device}"
                )

        elif job_type == "folding":
            if self.device == DEVICE_LILIBET:
                job_status = self.submit_job_folding_lilibet(job_details)
            else:
                raise ValueError(
                    f"Folding on other devices not implemented yet: {self.device}"
                )
        if job_status is True:
            self.ref.child(
                self.find_db_task_for_file_name(job_details["file_name"])
            ).child(f"{job_type}_status").set(STATUS_COMPLETED)
        else:
            self.ref.child(
                self.find_db_task_for_file_name(job_details["file_name"])
            ).child(f"{job_type}_status").set(STATUS_NOT_STARTED)
            self.ref.child(
                self.find_db_task_for_file_name(job_details["file_name"])
            ).child(f"{job_type}_device").set(STATUS_UNASSIGNED)


if __name__ == "__main__":

    @command()
    @option("--config_file_path", help="Path to the config file")
    @option("--device", default="lilibet", help="Device to run the scheduler on")
    @option("--max_jobs", default=100, help="Max number of jobs to run")
    @option(
        "--task_type",
        default="folding",
        help="Type of task to run. Choose between 'folding' and 'msa'",
    )
    def run_scheduler(
        config_file_path: str, device_name: str, max_jobs: int, task_type: str
    ):
        # Initialize the TaskScheduler
        ts = TaskScheduler(config_file_path, device_name)

        print("Starting Scheduler...")
        for job_idx in range(1, max_jobs + 1):
            # Obtain the task which isn't started yet.
            if task_type == JOB_TYPE_MSA:
                filter_criteria = {
                    "msa_status": STATUS_NOT_STARTED,
                }
            elif task_type == JOB_TYPE_FOLDING:
                filter_criteria = {
                    "folding_status": STATUS_NOT_STARTED,
                    "msa_status": STATUS_COMPLETED,
                }
            else:
                raise ValueError(f"Unknown task type: {task_type}")

            db_tasks = ts.get_tasks_by_filter(filter_criteria)
            if len(db_tasks) == 0:
                print("No tasks to run")
                break
            # Run the task
            print(f"{job_idx}. Running task: {db_tasks[0]}")
            ts.run_job(db_tasks[0], task_type)
