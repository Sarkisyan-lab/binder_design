import os
from os import path
import subprocess
import firebase_admin

import json

config = {
    "HPC_BASE": "hpc_base",
    "HPC_CX3": "hpc_cx3",
    "LILIBET": "lilibet",
    "JEX": "jex",

    "HPC_PASS": "suwGyb-qasquq-7sephu",
    "HPC_BASE_HOST": "ha1822@login.hpc.imperial.ac.uk",
    "HPC_CX3_HOST": "ha1822@login.cx3.hpc.ic.ac.uk",
    "JEX_HOST": "afleiss@jex.cscdom.csc.mrc.ac.uk",
    "LILIBET_HOST": "harsh@fw5.sshreach.me",

    "LILIBET_OUTPUT_DIR": "/media/HDD3/harsh/v8/colabfold",
    "HPC_OUTPUT_DIR": "/rds/general/user/ha1822/ephemeral/outputs/colabfold/misc",
    "HPC_COLABFOLD_CONDA_ENV": "/rds/general/user/ha1822/home/code/sarkisyan_lab/localcolabfold/colabfold-conda",
    "HPC_DBS_LOC": "/rds/general/user/ha1822/ephemeral/colabsearch",
    "JEX_OUTPUT_DIR": "/home/afleiss/harsh/sarkisyan_lab/outputs/v1/colabfold",
    "JEX_COLABFOLD_CONDA_ENV": "/home/afleiss/harsh/sarkisyan_lab/localcolabfold/colabfold-conda",

    "HPC_MSA_JOB_NUM_CPUS_LOCAL": 128,
    "HPC_MSA_JOB_MEM_GB_LOCAL": 512,
    "HPC_MSA_JOB_TIME_LOCAL": "02:00:00",

    "HPC_GPU_JOB_NUM_CPUS": 8,
    "HPC_GPU_JOB_MEM_GB": 64,
    "HPC_GPU_JOB_TIME": "02:00:00",

    "COLABFOLD_NUM_MODELS": 2,
    "COLABFOLD_DROPOUT": True,
    "COLABFOLD_USE_GPU_RELAX": True,
    "COLABFOLD_NUM_RECYLE": 10
}

def scp(source: str, dest: str, from_server: str, to_server: str):
    if os.path.isfile(source):
        use_r_flag = ""
    else:
        use_r_flag = "-r"
    if from_server == config["LILIBET"] and to_server == config["HPC_BASE"]:
        try:
            subprocess.run(
                f"sshpass -p {config["HPC_PASS"]} scp {use_r_flag} {source} {config["HPC_BASE_HOST"]}:{dest}",
                shell=True,
            )
        except Exception as e:
            raise Exception(
                f"Error copying file from {from_server} to {to_server}: {e}"
            )
    elif from_server == config["HPC_BASE"] and to_server == config["LILIBET"]:
        try:
            subprocess.run(
                f"sshpass -p {config["HPC_PASS"]} scp {use_r_flag} {config["HPC_BASE_HOST"]}:{source} {dest}",
                shell=True,
            )
        except Exception as e:
            raise Exception(
                f"Error copying file from {from_server} to {to_server}: {e}"
            )
    elif from_server == config["LILIBET"] and to_server == config["JEX"]:
        try:
            subprocess.run(
                f"scp {use_r_flag} {source} {config["JEX_HOST"]}:{dest}",
                shell=True,
            )
        except Exception as e:
            raise Exception(
                f"Error copying file from {from_server} to {to_server}: {e}"
            )
    elif from_server == config["JEX"] and to_server == config["LILIBET"]:
        try:
            subprocess.run(
                f"scp {use_r_flag} {config["JEX_HOST"]}:{source} {dest}",
                shell=True,
            )
        except Exception as e:
            raise Exception(
                f"Error copying file from {from_server} to {to_server}: {e}"
            )


def mkdir(dir: str, server: str):
    if server == config["LILIBET"]:
        os.makedirs(dir, exist_ok=True)
    elif server == config["HPC_BASE"]:
        try:
            subprocess.run(
                f"sshpass -p {config["HPC_PASS"]} ssh {config["HPC_BASE_HOST"]} mkdir -p {dir}", shell=True
            )
        except Exception as e:
            raise Exception(f"Error creating directory on {server}: {e}")
    elif server == config["JEX"]:
        try:
            subprocess.run(f"ssh {config["JEX_HOST"]} mkdir -p {dir}", shell=True)
        except Exception as e:
            raise Exception(f"Error creating directory on {server}: {e}")


def rm(file: str, server: str):
    if server == config["HPC_BASE"]:
        try:
            subprocess.run(
                f"sshpass -p {config["HPC_PASS"]} ssh {config["HPC_BASE_HOST"]} rm -rf {file}", shell=True
            )
        except Exception as e:
            raise Exception(f"Error removing file on {server}: {e}")


def submit_job_cpu_hpc_cx3_local(fasta_file_path: str):
    fasta_file_name = path.basename(fasta_file_path)
    job_name = path.splitext(fasta_file_name)[0]

    # create the job directory and copy the fasta file
    mkdir(dir=f"{config["HPC_OUTPUT_DIR"]}/{job_name}/input/", server=config["HPC_BASE"])
    mkdir(dir=f"{config["HPC_OUTPUT_DIR"]}/{job_name}/output/", server=config["HPC_BASE"])
    mkdir(dir=f"{config["HPC_OUTPUT_DIR"]}/{job_name}/logs", server=config["HPC_BASE"])
    scp(
        source=fasta_file_path,
        dest=f"{config["HPC_OUTPUT_DIR"]}/{job_name}/input/",
        from_server=config["LILIBET"],
        to_server=config["HPC_BASE"],
    )

    commands = [
        "#!/bin/bash",
        f"#PBS -l select=1:ncpus={config["HPC_MSA_JOB_NUM_CPUS_LOCAL"]}:mem={config["HPC_MSA_JOB_MEM_GB_LOCAL"]}gb",
        f"#PBS -l walltime={config["HPC_MSA_JOB_TIME_LOCAL"]}",
        f"#PBS -N {job_name}",
        f"#PBS -e {config["HPC_OUTPUT_DIR"]}/{job_name}/logs/colabfold_batch.err",
        f"#PBS -o {config["HPC_OUTPUT_DIR"]}/{job_name}/logs/colabfold_batch.out",
        "exec > >(tee -a $PBS_O_WORKDIR/logs/$PBS_JOBNAME.out)",
        "cd $PBS_O_WORKDIR",
        'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
        f"conda activate {config["HPC_COLABFOLD_CONDA_ENV"]}",
        f"colabfold_search {config["HPC_OUTPUT_DIR"]}/{job_name}/input/{fasta_file_name} {config["HPC_DBS_LOC"]} {config["HPC_OUTPUT_DIR"]}/{job_name}/output/msas",
        f"scp -P 10002 -r {config["HPC_OUTPUT_DIR"]}/{job_name}/output {config["LILIBET_HOST"]}:{config["LILIBET_OUTPUT_DIR"]}/{job_name}/",
    ]

    with open(f"{config["LILIBET_OUTPUT_DIR"]}/temp/{job_name}.pbs", "w") as f:
        for command in commands:
            f.write(command + "\n")

    scp(
        source=f"{config["LILIBET_OUTPUT_DIR"]}/temp/{job_name}.pbs",
        dest=f"{config["HPC_OUTPUT_DIR"]}/{job_name}",
        from_server=config["LILIBET"],
        to_server=config["HPC_BASE"],
    )

    # submit_job
    try:
        full_cmd = " ".join(
            [
                f"sshpass -p {config["HPC_PASS"]} ssh -t {config["HPC_BASE_HOST"]}",
                f'"eval \\"\$(~/miniconda3/bin/conda shell.bash hook)\\"',
                f"&& conda activate base && sshpass -p {config["HPC_PASS"]} ssh {config["HPC_CX3_HOST"]}",
                f"'source /etc/profile; cd {config["HPC_OUTPUT_DIR"]}/{job_name} && qsub {job_name}.pbs'\"",
            ]
        )
        subprocess.run(full_cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error submitting job on HPC: {e}")


def submit_job_gpu_hpc_cx3(fasta_file_path: str):
    # Check if the output .pkl file exists
    fasta_file_name = path.basename(fasta_file_path)
    job_name = path.splitext(fasta_file_name)[0]
    if os.path.isfile(f"{config["LILIBET_OUTPUT_DIR"]}/{job_name}/{job_name}.pickle"):
        print("Pickle found locally.")
        # MSA .pkl file exists
        # Make output director in HPC
        mkdir(f"{config["HPC_OUTPUT_DIR"]}/{job_name}/output/", config["HPC_BASE"])
        scp(
            source=f"{config["LILIBET_OUTPUT_DIR"]}/{job_name}/{job_name}.pickle",
            dest=f"{config["HPC_OUTPUT_DIR"]}/{job_name}/output/",
            from_server=config["LILIBET"],
            to_server=config["HPC_BASE"],
        )

    commands = [
        "#!/bin/bash",
        f"#PBS -l select=1:ncpus={config["HPC_GPU_JOB_NUM_CPUS"]}:mem={config["HPC_GPU_JOB_MEM_GB"]}gb:ngpus=1",
        f"#PBS -l walltime={config["HPC_GPU_JOB_TIME"]}",
        "cd $PBS_O_WORKDIR",
        'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
        f"conda activate {config["HPC_COLABFOLD_CONDA_ENV"]}",
        f"colabfold_batch --amber --overwrite-existing-results --num-models 3 --num-recyle 5 {config["HPC_OUTPUT_DIR"]}/{job_name}/input/{fasta_file_name} {config["HPC_OUTPUT_DIR"]}/{job_name}/output/",
        f"scp -P 10002 -r {config["HPC_OUTPUT_DIR"]}/{job_name}/output harsh@fw5.sshreach.me:{config["LILIBET_OUTPUT_DIR"]}/{job_name}/",
    ]

    with open(f"{config["LILIBET_OUTPUT_DIR"]}/temp/{job_name}_gpu.pbs", "w") as f:
        for command in commands:
            f.write(command + "\n")

    scp(
        source=f"{config["LILIBET_OUTPUT_DIR"]}/temp/{job_name}_gpu.pbs",
        dest=f"{config["HPC_OUTPUT_DIR"]}/{job_name}",
        from_server=config["LILIBET"],
        to_server=config["HPC_BASE"],
    )

    # submit_job
    try:
        full_cmd = " ".join(
            [
                f"sshpass -p {config["HPC_PASS"]} ssh -t {config["HPC_BASE_HOST"]}",
                f'"eval \\"\$(~/miniconda3/bin/conda shell.bash hook)\\"',
                f"&& conda activate base && sshpass -p {config["HPC_PASS"]} ssh {config["HPC_CX3_HOST"]}",
                f"'source /etc/profile; cd {config["HPC_OUTPUT_DIR"]}/{job_name} && qsub {job_name}.pbs'\"",
            ]
        )
        subprocess.run(full_cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error submitting job on HPC: {e}")


def submit_job_gpu_jex(fasta_file_path: str):
    # Check if the output .pkl file exists
    fasta_file_name = path.basename(fasta_file_path)
    job_name = path.splitext(fasta_file_name)[0]
    if os.path.isfile(f"{config["LILIBET_OUTPUT_DIR"]}/{job_name}/{job_name}.pickle"):
        print("Pickle found locally.")
        # MSA .pkl file exists
        # Make output director in HPC
        mkdir(f"{config["JEX_OUTPUT_DIR"]}/{job_name}/output/", config["JEX"])
        scp(
            source=f"{config["LILIBET_OUTPUT_DIR"]}/{job_name}/{job_name}.pickle",
            dest=f"{config["JEX_OUTPUT_DIR"]}/{job_name}/output/",
            from_server=config["LILIBET"],
            to_server=config["JEX"],
        )

    # commands = [
    #     "#!/bin/bash",
    #     f"#PBS -l select=1:ncpus={HPC_GPU_JOB_NUM_CPUS}:mem={HPC_GPU_JOB_MEM_GB}gb:ngpus=1",
    #     f"#PBS -l walltime={HPC_GPU_JOB_TIME}",
    #     "cd $PBS_O_WORKDIR",
    #     'eval "$(~/miniconda3/bin/conda shell.bash hook)"',
    #     f"conda activate {HPC_COLABFOLD_CONDA_ENV}",
    #     f"colabfold_batch --amber --overwrite-existing-results --num-models 3 --num-recyle 5 {HPC_OUTPUT_DIR}/{job_name}/input/{fasta_file_name} {HPC_OUTPUT_DIR}/{job_name}/output/",
    #     f"scp -P 10002 -r {HPC_OUTPUT_DIR}/{job_name}/output harsh@fw5.sshreach.me:{LILIBET_OUTPUT_DIR}/{job_name}/",
    # ]

    # with open(f"{LILIBET_OUTPUT_DIR}/temp/{job_name}_gpu.pbs", "w") as f:
    #     for command in commands:
    #         f.write(command + "\n")

    # scp(
    #     source=f"{LILIBET_OUTPUT_DIR}/temp/{job_name}_gpu.pbs",
    #     dest=f"{HPC_OUTPUT_DIR}/{job_name}",
    #     from_server=LILIBET,
    #     to_server=HPC_BASE,
    # )

class TaskScheduler:
    def __init__(self):
        cred_obj = firebase_admin.credentials.Certificate('/home/harsh/code/sarkisyan_lab/binder_design/src/login_key.json')
        self.current_jobs = {}

    def inflate_jobs(self, fasta_dir: str):
        for fasta_file in os.listdir(fasta_dir):
            if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
    def add_job(self, job_id: str, job_type: str):
        self.current_jobs[job_id] = job_type

    def remove_job(self, job_id: str):
        del self.current_jobs[job_id]
        