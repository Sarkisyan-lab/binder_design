import os
from os import path
import subprocess
import firebase_admin
from firebase_admin import db
import glob
import json
from argparse import Namespace
from click import command, option


class TaskScheduler:
    def __init__(self, config: dict, device="lilibet"):
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
                        "msa_status": "not_started",
                        "folding_status": "not_started",
                        "msa_device": "unassigned",
                        "folding_device": "unassigned",
                    }  # type: ignore
                )

    def get_tasks_by_filter(self, filter_criteria: dict) -> list:
        objs = self.ref.get()
        filtered_tasks = []
        for key, value in objs.items():
            if all(value[k] == v for k, v in filter_criteria.items()):
                filtered_tasks.append(value)
        return filtered_tasks

    def get_local_tasks_status(self, fasta_dir: str) -> dict:
        local_task_status = []
        for fasta_file in glob.glob(os.path.join(fasta_dir, "*.fasta")) + glob.glob(
            os.path.join(fasta_dir, "*.fa")
        ):
            if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
                file_name = os.path.basename(fasta_file).split(".")[0]
                task_status = {
                    "file_name": file_name,
                    "seq": self.fasta_file_to_fasta_str(fasta_file),
                    "msa_status": "not_started",
                    "folding_status": "not_started",
                    "msa_device": "unassigned",
                    "folding_device": "unassigned",
                }

                # If folding / msa is completed
                ind_output_dir = os.path.join(
                    self.config.lilibet_output_dir, file_name, "output"
                )
                if os.path.exists(ind_output_dir):
                    folder_contents = os.listdir(ind_output_dir)
                    for file in folder_contents:
                        if file.endswith(".a3m") or file.endswith(".pickle"):
                            task_status["msa_status"] = "completed"
                            task_status["msa_device"] = self.device
                        elif "rank_001_alphafold2_multimer_v3" in file:
                            task_status["folding_status"] = "completed"
                            task_status["folding_device"] = self.device
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
                    local_task[f"{job_type}_status"] == "completed"
                    and resp_db_task_value[f"{job_type}_status"] != "completed"
                ):
                    # Update the db with the completed task
                    self.ref.child(resp_db_task_key).child(f"{job_type}_status").set(
                        "completed"
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
            "exec > >(tee -a $PBS_O_WORKDIR/logs/$PBS_JOBNAME_msa_local.out)",
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
        try:
            subprocess.run(f"qsub {job_file_path}", shell=True)
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error submitting job on HPC: {e}")

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
        if job_type == "msa":
            if self.device == "lilibet":
                self.submit_job_msa_server_lilibet(job_details)
            elif self.device == "hpc_cx3" or self.device == "hpc_base":
                self.submit_job_msa_local_hpc_cx3_local(job_details)
            else:
                raise ValueError(
                    f"MSAs on other devices not implemented yet: {self.device}"
                )

        elif job_type == "folding":
            if self.device == "lilibet":
                self.submit_job_folding_lilibet(job_details)
            else:
                raise ValueError(
                    f"Folding on other devices not implemented yet: {self.device}"
                )


if __name__ == "__main__":

    @command()
    @option("--config_file_path", help="Path to the config file")
    @option("--device", default="lilibet", help="Device to run the scheduler on")
    @option("--max_jobs", default=100, help="Max number of jobs to run")
    @option("--task_type", default="folding", help="Type of task to run")
    def run_scheduler(
        config_file_path: str, device_name: str, max_jobs: int, task_type: str
    ):
        # Initialize the TaskScheduler
        ts = TaskScheduler(config_file_path, device_name)

        print("Starting Scheduler")
        for job_idx in range(1, max_jobs + 1):
            # Obtain the task which isn't started yet.
            db_tasks = ts.get_tasks_by_filter({f"{task_type}_status": "not_started"})
            if len(db_tasks) == 0:
                print("No tasks to run")
                break
            # Run the task
            print(f"{job_idx}. Running task: {db_tasks[0]}")
            ts.run_job(db_tasks[0], task_type)
