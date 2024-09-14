import os
import subprocess
import firebase_admin
from firebase_admin import db
import glob
import pandas as pd
from argparse import Namespace, ArgumentParser
import json
import time
import logging
from dotenv import load_dotenv

import sys

sys.path.append(os.path.abspath(os.path.join(__file__, "../")))
import utils

# constants
STATUS_UNASSIGNED = "unassigned"
STATUS_QUEUED = "queued"
STATUS_RUNNING = "running"
STATUS_NOT_STARTED = "not_started"
STATUS_COMPLETED = "completed"
DEVICE_LILIBET = "lilibet"
DEVICE_HPC_CX3 = "hpc_cx3"
DEVICE_HPC_BASE = "hpc_base"
DEVICE_HPC_HX1 = "hpc_hx1"
DEVICE_JEX = "jex"
JOB_TYPE_MSA = "msa"
JOB_TYPE_FOLDING = "folding"

# load dotenv
load_dotenv(os.path.abspath(os.path.join(__file__, "../../.env")))


class TaskScheduler:
    def __init__(self, args: Namespace = None):
        """
        Args:
            args (Namespace): Arguments from the command line. Not required if not submitting jobs.
            logging_level (str): Logging level (info/debug). Defaults to "info".
        """
        self.logger = logging.getLogger("scheduler_body")
        self.logger.setLevel(logging.INFO)

        # If no args were provided
        self.args = args
        if self.args is None:
            self.logger.warning("No arguments provided.")

        cred_obj = firebase_admin.credentials.Certificate(
            os.path.abspath(os.path.join(__file__, "../src/login_key.json"))
        )
        self.default_app = firebase_admin.initialize_app(
            cred_obj, {"databaseURL": os.environ["FIREBASE_DB_URL"]}
        )
        self.ref = db.reference("/")

        # Args into self vars
        self.device = utils.obtain_device_name()
        self.output_dir = args.output_dir

    def empty_db(self):
        """
        Empties the entire database
        """
        objs = self.ref.get()
        for key, _ in objs.items():
            self.ref.child(key).set({})
        self.logger.info(f"Emptied the database containing {len(objs)} tasks.")

    # def inflate_tasks_from_fasta_dir(
    #     self, fasta_dir, set_default_msas_status=STATUS_COMPLETED
    # ):
    #     """
    #     Inflates the database with tasks from a fasta directory

    #     Args:
    #         fasta_dir (str): Path to the fasta directory
    #         set_default_msas_status (str): Status to set for the MSAs. Defaults to STATUS_COMPLETED.
    #     """
    #     for fasta_file in glob.glob(
    #         os.path.join(self.output_dir, "input", "*.fasta")
    #     ) + glob.glob(os.path.join(self.output_dir, "input", "*.fa")):
    #         fasta_content = self.fasta_file_to_fasta_str(fasta_file)
    #         self.logger.info(f"Adding task: {os.path.basename(fasta_file)}")
    #         self.ref.push(
    #             {
    #                 "file_name": os.path.basename(fasta_file).split(".")[0],
    #                 "seq": fasta_content,
    #                 "msa_status": set_default_msas_status,
    #                 "folding_status": STATUS_NOT_STARTED,
    #                 "msa_device": STATUS_UNASSIGNED,
    #                 "folding_device": STATUS_UNASSIGNED,
    #             }
    #         )

    def get_tasks_by_filter(self, filter_criteria: dict) -> list:
        """
        Gets tasks from the database based on the filter criteria

        Args:
            filter_criteria (dict): Filter criteria

        Returns:
            list: List of tasks
        """
        objs = self.ref.get()
        filtered_tasks = []
        for key, value in objs.items():
            if all(value[k] == v for k, v in filter_criteria.items()):
                filtered_tasks.append(value)
        self.logger.info(
            f"Found {len(filtered_tasks)} tasks matching the filter criteria"
        )
        return filtered_tasks

    def get_local_tasks_status(self) -> list:
        """
        Get the status of the tasks on the local machine

        Returns:
            list: List of tasks with their status
        """
        local_task_status = []
        msas_dir = os.path.join(self.output_dir, "msas")
        fasta_dir = os.path.join(self.output_dir, "input")

        if not os.path.exists(msas_dir):
            msas_folder_contents = []
        else:
            msas_folder_contents = os.listdir(msas_dir)

        for fasta_file in glob.glob(os.path.join(fasta_dir, "*.fasta")) + glob.glob(
            os.path.join(fasta_dir, "*.fa")
        ):

            file_name = os.path.basename(fasta_file).split(".")[0]
            task_status = {
                "file_name": file_name,
                "seq": utils.fasta_file_to_fasta_str(fasta_file),
                "msa_status": STATUS_NOT_STARTED,
                "folding_status": STATUS_NOT_STARTED,
                "msa_device": STATUS_UNASSIGNED,
                "folding_device": STATUS_UNASSIGNED,
            }
            predictions_dir = os.path.join(self.output_dir, "predictions", file_name)
            # Check if MSAs exist
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

    def find_db_id_for_file_name(self, file_name: str, db_tasks: dict = None) -> dict:
        """
        Find the database id for a given file name

        Args:
            file_name (str): File name
            db_tasks (dict): Database tasks. If None, will fetch from the database.

        Returns:
            dict: Database id
        """
        if db_tasks is None:
            db_tasks = self.ref.get("/")[0]
        for db_task_key, db_task_value in db_tasks.items():
            if db_task_value["file_name"] == file_name:
                return db_task_key
        return None

    def overwrite_db_with_local_tasks(self, local_tasks_status: list):
        """
        Update the database with the local tasks

        Args:
            local_tasks_status (list): List of tasks with their status
        """

        db_tasks = self.ref.get()
        for local_task in local_tasks_status:
            # Find the task in the db for a given file name
            resp_db_task_key = self.find_db_id_for_file_name(
                local_task["file_name"], db_tasks
            )
            db_task_value = db_tasks[resp_db_task_key]
            if db_task_value is None:
                self.logger.info(
                    f"Task {local_task['file_name']} not found in db. Inserting."
                )
                self.ref.push(local_task)
                continue

            for job_type in ["msa", "folding"]:
                if (
                    local_task[f"{job_type}_status"]
                    != db_task_value[f"{job_type}_status"]
                ):
                    self.logger.info(
                        f"Updating {job_type} status for {local_task['file_name']} from {db_task_value[f'{job_type}_status']} to {local_task[f'{job_type}_status']}"
                    )
                    # Update the db with the completed task
                    self.ref.child(resp_db_task_key).child(f"{job_type}_status").set(
                        local_task[f"{job_type}_status"]
                    )
                    self.ref.child(resp_db_task_key).child(f"{job_type}_device").set(
                        self.device
                    )

    def submit_job_msa_local_hpc_cx3_local(self, job_details: dict):
        """
        Submit a job to the HPC to run MSA generation on the local machine

        Args:
            job_details (dict): Job details
        """
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
                fasta_data.append({"id": description, "sequence": sequence})

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
            f"exec > >(tee -a {logs_dir}/" + "colabfold_search_batch_${PBS_JOBID}.out)",
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

    def job_housekeeping(self, job_details: dict):
        """
        Perform housekeeping tasks for a job.
        
        Args:
            job_details (dict): Job details.
        """
        fasta_file_name = job_details["file_name"]

        # Check if fasta exists locally
        fasta_file_path = os.path.join(
            self.output_dir,
            "input",
            f"{fasta_file_name}.fasta",
        )
        msas_file_path = os.path.join(self.output_dir, "msas")
        predictions_dir = os.path.join(self.output_dir, "predictions", fasta_file_name)

        assert utils.check_if_fasta_file_exists(fasta_file_path), f"Fasta file not found: {fasta_file_name}"
        assert utils.check_if_msas_exist(msas_file_path, fasta_file_name), f"MSAs do not exist for {fasta_file_name}!"

        # Make predictions directory
        os.makedirs(predictions_dir, exist_ok=True)

    # def submit_job_folding_hpc_cx3(
    #     self, job_details_list: list, copy_to_lilibet: bool = True
    # ):
    #     print(
    #         f"Running Folding for {len(job_details_list)} sequences on {self.device} parallely"
    #     )
    #     fold_cmd = []
    #     scp_cmd = [] if copy_to_lilibet else None
    #     for job_details in job_details_list:
    #         fasta_file_name = job_details["file_name"]
    #         logs_dir, predictions_dir, msa_output_dir, fasta_file_full_path = (
    #             self.job_housekeeping(job_details, device=DEVICE_HPC_CX3)
    #         )
    #         msa_file_path = os.path.join(msa_output_dir, f"{fasta_file_name}.a3m")
    #         template_file_path = os.path.join(
    #             msa_output_dir, f"{fasta_file_name}_pdb100_230517.m8"
    #         )
    #         drop_out_str = "--use-dropout" if self.config.colabfold_dropout else ""
    #         fold_cmd.append(
    #             f"colabfold_batch {msa_file_path} {predictions_dir} --num-recycle {self.config.colabfold_num_recycle} --num-models {self.config.colabfold_num_models} {drop_out_str}"
    #         )
    #         if copy_to_lilibet:
    #             scp_cmd.append(
    #                 f"scp -r -P 10002 {predictions_dir} {self.config.lilibet_host}:{self.config.lilibet_output_dir}/predictions/"
    #             )
    #     fold_cmd_str = " & ".join(fold_cmd)
    #     scp_cmd_str = " & ".join(scp_cmd) if copy_to_lilibet else None

    #     # Folding commands
    #     commands = [
    #         "#!/bin/bash",
    #         f"#PBS -l select=1:ncpus={self.config.hpc_folding_job_num_cpus}:mem={self.config.hpc_folding_job_mem_gb}gb:ngpus=1",
    #         f"#PBS -l walltime={self.config.hpc_folding_job_time}",
    #         f"#PBS -N {fasta_file_name}",
    #         f"#PBS -e {logs_dir}/",
    #         f"#PBS -o {logs_dir}/",
    #         f"exec > >(tee -a {logs_dir}/" + "${PBS_JOBNAME}_${PBS_JOBID}_folding.out)",
    #         "cd $PBS_O_WORKDIR",
    #         'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
    #         f"conda activate {self.config.hpc_colabfold_conda_env}",
    #         fold_cmd_str,
    #         "wait",
    #     ]
    #     if scp_cmd_str:
    #         commands.append(scp_cmd_str)

    #     job_file_path = os.path.join(logs_dir, f"{fasta_file_name}_folding.pbs")
    #     with open(job_file_path, "w") as f:
    #         for command in commands:
    #             f.write(command + "\n")

    #     # submit_job
    #     try:
    #         subprocess.run(f"qsub {job_file_path}", shell=True)
    #     except subprocess.CalledProcessError as e:
    #         raise Exception(f"Error submitting job on HPC: {e}")
    #     return True

    # def submit_job_folding_hpc_hx1(
    #     self, job_details_list: list, copy_to_lilibet: bool = True
    # ):
    #     print(
    #         f"Running Folding for {len(job_details_list)} sequences on {self.device} parallely"
    #     )
    #     fold_cmd = []
    #     scp_cmd = [] if copy_to_lilibet else None
    #     for job_details in job_details_list:
    #         fasta_file_name = job_details["file_name"]
    #         logs_dir, predictions_dir, msa_output_dir, fasta_file_full_path = (
    #             self.job_housekeeping(job_details, device=DEVICE_HPC_HX1)
    #         )
    #         msa_file_path = os.path.join(msa_output_dir, f"{fasta_file_name}.a3m")
    #         template_file_path = os.path.join(
    #             msa_output_dir, f"{fasta_file_name}_pdb100_230517.m8"
    #         )
    #         drop_out_str = "--use-dropout" if self.config.colabfold_dropout else ""
    #         fold_cmd.append(
    #             f"colabfold_batch {msa_file_path} {predictions_dir} --num-recycle {self.config.colabfold_num_recycle} --num-models {self.config.colabfold_num_models} {drop_out_str}"
    #         )
    #         if copy_to_lilibet:
    #             scp_cmd.append(
    #                 f"scp -r -P 10002 {predictions_dir} {self.config.lilibet_host}:{self.config.lilibet_output_dir}/predictions/"
    #             )
    #     fold_cmd_str = " & ".join(fold_cmd)
    #     scp_cmd_str = " & ".join(scp_cmd) if copy_to_lilibet else None

    #     # Folding commands
    #     commands = [
    #         "#!/bin/bash",
    #         f"#PBS -l select=1:ncpus={self.config.hpc_folding_job_num_cpus}:mem={self.config.hpc_folding_job_mem_gb}gb:ngpus=1:gpu_type=A100",
    #         f"#PBS -l walltime={self.config.hpc_folding_job_time}",
    #         f"#PBS -N {fasta_file_name}",
    #         f"#PBS -e {logs_dir}/",
    #         f"#PBS -o {logs_dir}/",
    #         f"exec > >(tee -a {logs_dir}/" + "${PBS_JOBNAME}_${PBS_JOBID}_folding.out)",
    #         "cd $PBS_O_WORKDIR",
    #         'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
    #         f"conda activate {self.config.hpc_hx1_colabfold_conda_env}",
    #         fold_cmd_str,
    #         "wait",
    #     ]
    #     if scp_cmd_str:
    #         commands.append(scp_cmd_str)

    #     job_file_path = os.path.join(logs_dir, f"{fasta_file_name}_folding.pbs")
    #     with open(job_file_path, "w") as f:
    #         for command in commands:
    #             f.write(command + "\n")

    #     # submit_job
    #     try:
    #         subprocess.run(f"qsub -q hx {job_file_path}", shell=True)
    #     except subprocess.CalledProcessError as e:
    #         raise Exception(f"Error submitting job on HPC: {e}")
    #     return True

    @staticmethod
    def return_hpc_queue_status(queue_name="v1_gpu72"):
        """
        Returns the number of running and queued jobs in the specified queue.
        """
        running_jobs = subprocess.run(
            ["qstat", "-r", "-u", "ha1822", queue_name],
            capture_output=True,
            text=True,
            check=True,
        ).stdout

        running_job_count = sum(
            1
            for line in running_jobs.splitlines()
            if line.strip() and line[0].isdigit()
        )

        queued_jobs = subprocess.run(
            ["qstat", "-i", "-u", "ha1822", queue_name],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        queued_job_count = sum(
            1 for line in queued_jobs.splitlines() if line.strip() and line[0].isdigit()
        )
        return running_job_count, queued_job_count

    # def run_batch_jobs(self, job_details_list: list):
    #     """
    #     Run the folding jobs for the given list of jobs.
        
    #     Args:
    #         job_details_list (list): List of job details to run.
    #     """
    #     # Set the job status to queued
    #     for job_idx, job_details in enumerate(job_details_list):
    #         file_name = job_details["file_name"]
    #         job_id = self.find_db_id_for_file_name(file_name)
    #         self.logger.info(f"Setting {file_name} to {STATUS_QUEUED}")
    #         self.ref.child(job_id).child(f"folding_status").set(STATUS_QUEUED)
    #         self.ref.child(job_id).child(f"folding_device").set(self.device)

    #         # Add id to job details
    #         job_details_list[job_idx]["job_id"] = job_id

    #     # Submit the job
    #     # TODO: Add support for JEX here.
    #     if self.device == DEVICE_LILIBET:
    #         self.submit_job_folding_lilibet(job_details_list)
    #     elif self.device == DEVICE_HPC_CX3:
    #         self.submit_job_folding_hpc_cx3(job_details_list)
    #     elif self.device == DEVICE_HPC_HX1:
    #         self.submit_job_folding_hpc_hx1(job_details_list)
    #     else:
    #         raise ValueError(
    #             f"Folding on other devices not implemented yet: {self.device}"
    #         )

    # def submit_job_folding_lilibet_copy(self, job_details_list: list):
    #     """
    #     Run folding jobs on lilibet
        
    #     Args:
    #         job_details_list (list): List of job details to run.
    #     """
    #     print(
    #         f"Running Folding for {len(job_details_list)} sequences on lilibet parallely"
    #     )
    #     folding_commands = []
    #     for job_details in job_details_list:
    #         # check if fasta, msas exist. create predictions dir.
    #         self.job_housekeeping(job_details)
            
    #         # get all directories and paths
    #         fasta_file_name = job_details["file_name"]
    #         predictions_dir = os.path.join(self.output_dir, "predictions", fasta_file_name)
    #         msa_file_path = os.path.join(self.output_dir, "msas", f"{fasta_file_name}.a3m")
            
    #         # TODO: add support for templates 
    #         template_file_path = os.path.join(
    #             self.output_dir, "msas", f"{fasta_file_name}_pdb100_230517.m8"
    #         )
    #         folding_commands.append(
    #             f"python wrapper.py --task_id={job_details['job_id']} --msa_file_path={msa_file_path} --predictions_dir={predictions_dir} --num-recycle={self.config.colabfold_num_recycle} --num-models={self.config.colabfold_num_models} --use-dropout={self.config.colabfold_dropout}"
    #         )

    #     folding_commands_str = " & ".join(folding_commands)

    #     # define a common set of commands
    #     colabfold_exec_path = os.path.abspath(os.path.dirname(__file__))
    #     commands = [
    #         f"source ~/anaconda3/etc/profile.d/conda.sh",
    #         f"conda activate {self.config.lilibet_colabfold_conda_env}",
    #         f"cd {colabfold_exec_path}",
    #         folding_commands_str,
    #     ]

    #     # return commands
    #     subprocess.run(
    #         ["bash", "-c", "\n".join(commands)],
    #         text=True,
    #     )
    #     return True

    def generate_initial_folding_commands(self, job_idx: int):
        """
        Generate the folding commands for the given job details list.
        """
        logs_dir = os.path.join(self.output_dir, "logs")
        if self.device == DEVICE_LILIBET:
            conda_env = os.getenv("LILIBET_COLABFOLD_CONDA_ENV")
            colabfold_exec_path = os.path.abspath(os.path.dirname(__file__))
            commands = [
                f"source ~/anaconda3/etc/profile.d/conda.sh",
                f"conda activate {conda_env}",
                f"cd {colabfold_exec_path}"
            ]
        elif self.device in [DEVICE_HPC_CX3, DEVICE_HPC_HX1]:
            if self.device == DEVICE_HPC_CX3:
                conda_env = os.getenv("HPC_CX3_COLABFOLD_CONDA_ENV")
            else:
                conda_env = os.getenv("HPC_HX1_COLABFOLD_CONDA_ENV")
            commands = [
                "#!/bin/bash",
                f"#PBS -l select=1:ncpus={self.config.hpc_folding_job_num_cpus}:mem={self.config.hpc_folding_job_mem_gb}gb:ngpus=1",
                f"#PBS -l walltime={self.config.hpc_folding_job_time}",
                f"#PBS -N job_{self.device}_{job_idx}",
                f"#PBS -e {logs_dir}/",
                f"#PBS -o {logs_dir}/",
                f"exec > >(tee -a {logs_dir}/job_{job_idx}" + "_${PBS_JOBID}_folding.out)",
                "cd $PBS_O_WORKDIR",
                'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
                f"conda activate {conda_env}",
            ]
        else:
            raise ValueError(f"Invalid device: {self.device}")
        return commands

    def submit_job_folding_common(self, job_details_list: list, job_index: int):
        """
        Run folding jobs on lilibet
        
        Args:
            job_details_list (list): List of job details to run.
        """
        folding_commands = []
        scp_cmds = []

        db_jobs = self.ref.get()
        for job_details in job_details_list:
            # Set the job status to queued
            file_name = job_details["file_name"]
            job_id = self.find_db_id_for_file_name(file_name, db_jobs)
            self.logger.info(f"Setting {file_name} to {STATUS_QUEUED}")
            self.ref.child(job_id).child(f"folding_status").set(STATUS_QUEUED)
            self.ref.child(job_id).child(f"folding_device").set(self.device)

            # check if fasta, msas exist. create predictions dir.
            self.job_housekeeping(job_details)
            
            # get all directories and paths
            fasta_file_name = job_details["file_name"]
            predictions_dir = os.path.join(self.output_dir, "predictions", fasta_file_name)
            msa_file_path = os.path.join(self.output_dir, "msas", f"{fasta_file_name}.a3m")
            
            # TODO: add support for templates 
            template_file_path = os.path.join(
                self.output_dir, "msas", f"{fasta_file_name}_pdb100_230517.m8"
            )
            folding_commands.append(
                f"python wrapper.py --task_id={job_id} --msa_file_path={msa_file_path} --predictions_dir={predictions_dir} --num-recycle={self.args.colabfold_num_recycle} --num-models={self.args.colabfold_num_models} --use-dropout={self.args.colabfold_dropout}"
            )

            if self.args.copy_to_lilibet:
                lilibet_host = os.getenv("LILIBET_HOST")
                lilibet_port = os.getenv("LILIBET_PORT")
                scp_cmds.append(
                    f"scp -r -P {lilibet_port} {predictions_dir} {lilibet_host}:{self.args.lilibet_output_dir}/predictions/"
                )

        folding_commands_str = " & ".join(folding_commands)
        scp_cmds_str = " & ".join(scp_cmds) if len(scp_cmds) > 0 else None

        # define a common set of commands
        commands = self.generate_initial_folding_commands(job_index)
        commands.append(folding_commands_str)
        if scp_cmds_str:
            commands.append("wait")
            commands.append(scp_cmds_str)

        if self.device == DEVICE_LILIBET:
            result = subprocess.run(
                ["bash", "-c", "\n".join(commands)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            
        elif self.device in [DEVICE_HPC_CX3, DEVICE_HPC_HX1]:
            hx1_queue_str = "-q hx" if self.device == DEVICE_HPC_HX1 else ""
            job_file_path = os.path.join(self.output_dir, f"logs/job_{self.device}_{job_index}_folding.pbs")
            
            # create pbs file
            with open(job_file_path, "w") as f:
                for command in commands:
                    f.write(command + "\n")
            
            subprocess.run(f"qsub {hx1_queue_str} {job_file_path}", shell=True)
        
        else:
            raise ValueError(f"Invalid device: {self.device}")


def main(args):
    """
    Main function to run the scheduler
    """
    device_name = utils.obtain_device_name()
    # Set up logging
    logger = logging.getLogger("scheduler_main")
    logging.basicConfig(level=logging.INFO)
    logger.info(f"Running scheduler on: {device_name}. Num jobs per GPU: {args.num_jobs_per_gpu}")

    if args.generate_msas:
        if device_name != DEVICE_HPC_CX3:
            raise ValueError(
                f"Generating MSAs is only supported on HPC CX3, not {device_name}"
            )
        pass
    elif args.clear_db_and_upload_msa_jobs:
        pass
    elif args.set_all_queued_jobs_to_not_started:
        pass

    if args.copy_to_lilibet:
        assert (
            args.lilibet_output_dir is not None
        ), "Lilibet output directory is required."

    ts = TaskScheduler(args)

    current_job_idx = 1
    while current_job_idx <= args.max_jobs + 1:        
        # Only for HPC / JEX
        # TODO: Add support for JEX here.
        if device_name in [DEVICE_HPC_HX1, DEVICE_HPC_CX3]:
            queue_name = "v1_a100" if device_name == DEVICE_HPC_HX1 else "v1_gpu72"
            running_jobs, queued_jobs = ts.return_hpc_queue_status(queue_name)
            
            # Don't run more and sleep
            if queued_jobs >= args.max_queued_jobs or (running_jobs + queued_jobs) >= 20:
                logger.info(f"Already {queued_jobs} jobs in queue; sleeping for {args.sleep_itvl_mins} mins.")
                time.sleep(args.sleep_itvl_mins * 60)
                continue

        # Obtain new tasks to run
        db_tasks = ts.get_tasks_by_filter(
            {
                "folding_status": STATUS_NOT_STARTED,
                "msa_status": STATUS_COMPLETED,
            }
        )
        if len(db_tasks) == 0:
            if running_jobs == 0:
                # Update db if all jobs are completed
                logger.info("No more tasks to run! Finishing the scheduler.")
                break
        
        # Run the task
        batch_jobs = db_tasks[:args.num_jobs_per_gpu]
        ts.submit_job_folding_common(batch_jobs, current_job_idx)
        logger.info(f"Finished running all jobs for job idx: {current_job_idx}")
        current_job_idx += 1


if __name__ == "__main__":
    # run_scheduler()
    argparse_bool = lambda x: (str(x).lower() == "true")
    parser = ArgumentParser(description="Arguments for the scheduler")

    parser.add_argument(
        "--generate_msas",
        type=argparse_bool,
        default=False,
        help="If set, will generate MSAs for all sequences present in input/ folder. Only applicable for HPC CX3.",
    )
    parser.add_argument(
        "--clear_db_and_upload_msa_jobs",
        type=argparse_bool,
        default=False,
        help="If true, will clear the database and upload the MSAs generated in the msas/ folder as completed jobs in firebase.",
    )
    parser.add_argument(
        "--set_all_queued_jobs_to_not_started",
        type=argparse_bool,
        default=False,
        help="If true, will set all queued jobs in the database to 'not_started' status. \nAll completed jobs must be present in the predictions/ folder.",
    )
    parser.add_argument(
        "--copy_to_lilibet",
        type=argparse_bool,
        default=False,
        help="Copy the output files to lilibet. If true, then lilibet_output_dir is required.",
    )
    parser.add_argument(
        "--lilibet_output_dir",
        type=str,
        default=None,
        help="Path to the output directory on lilibet. Required if --copy_to_lilibet is true",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Path to the output directory on the current machine.",
    )
    parser.add_argument(
        "--hpc_msa_num_cpus",
        type=int,
        default=72,
        help="Number of CPUs to use for Colabfold search on HPC CX3",
    )
    parser.add_argument(
        "--hpc_msa_mem_gb",
        type=int,
        default=500,
        help="Memory in GB to use for Colabfold search on HPC CX3",
    )
    parser.add_argument(
        "--hpc_msa_time",
        type=str,
        default="04:00:00",
        help="Time limit for Colabfold search on HPC CX3",
    )
    parser.add_argument(
        "--hpc_folding_time",
        type=str,
        default="01:30:00",
        help="Time limit for Colabfold folding on HPC CX3/HX1",
    )
    parser.add_argument(
        "--hpc_folding_num_cpus",
        type=int,
        default=8,
        help="Number of CPUs to use for Colabfold folding on HPC CX3/HX1",
    )
    parser.add_argument(
        "--hpc_folding_mem_gb",
        type=int,
        default=64,
        help="Memory in GB to use for Colabfold folding on HPC CX3/HX1",
    )
    parser.add_argument(
        "--colabfold_num_recycle",
        type=int,
        default=3,
        help="Number of recycling cycles to use for Colabfold",
    )
    parser.add_argument(
        "--colabfold_num_models",
        type=int,
        default=5,
        help="Number of models to use for Colabfold",
    )
    parser.add_argument(
        "--colabfold_dropout",
        type=argparse_bool,
        default=True,
        help="If true, will use dropout for Colabfold",
    )
    parser.add_argument(
        "--max_jobs",
        type=int,
        default=100,
        help="Max number of queued jobs",
    )
    parser.add_argument(
        "--max_queued_jobs",
        type=int,
        default=5,
        help="Max number of queued jobs",
    )
    parser.add_argument(
        "--sleep_itvl_mins",
        type=int,
        default=2,
        help="Sleep interval in minutes",
    )
    parser.add_argument(
        "--num_jobs_per_gpu",
        type=int,
        default=2,
        help="Number of jobs to run per GPU",
    )

    args = parser.parse_args()
    main(args)
