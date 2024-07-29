import os, time, glob
not_done_2 = []

base_path = "/rds/general/user/ha1822/home/code/jobs_ligand/v5/"
fasta_path = base_path + "ligand_mpnn/processed/"
parafold_path = base_path + "parafold/"
msa_path = "/rds/general/user/ha1822/home/code/jobs_ligand/jobs/msas"

pbs_path = "/rds/general/user/ha1822/home/code/jobs_ligand/jobs/"

not_done = os.listdir(fasta_path)
for i in not_done:
    temp_string = f"""#!/bin/bash
#PBS -N jp_m2c590_{i[-5:-3]}
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=48:mem=512gb

cd ${{PBS_O_WORKDIR}}

module load anaconda3/personal

## Remember to change this value yo the name of your VEnv with pytorch
source activate parafold

echo "Starting Job...$(date)"

bash /rds/general/user/ha1822/home/code/sarkisyan_lab/ParallelFold/run_alphafold.sh -d /rds/general/user/ha1822/ephemeral/afdb -o {parafold_path} -m model_1_multimer_v3 -p multimer -i {fasta_path}{i} -t 2020-12-01 --Y common_chain -Z {msa_path} -f

echo "... Run finished $(date)"

    """
    with open(f"{pbs_path}/test_job.pbs", 'w') as f:
        f.write(temp_string)
    commands = [f"cd {pbs_path}", "qsub test_job.pbs"]
    result = os.popen('\n'.join(commands)).read()
    if not result.strip().endswith('.pbs'):
        not_done.append(i)
        for j in glob.glob(parafold_path + '*/*'):
            if j.endswith('features.pkl'):
                os.system(f"rsync -rv -e 'ssh -p 10002' --progress --ignore-existing {parafold_path}/{i.split('/')[-2]}/features.pkl harsh@fw5.sshreach.me:/media/HDD3/harsh/v5/parafold/{i.split('/')[-2]}/features.pkl")
        time.sleep(600)
        
    
    if len(not_done) == 0:
        break