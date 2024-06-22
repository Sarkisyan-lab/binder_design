import os
import glob
import shutil
import argparse
from argparse import Namespace
import pandas as pd
import subprocess
import wandb
from dotenv import load_dotenv
import numpy as np


# load .env
load_dotenv()

CHAIN_B_DESC = "\n>ChainB|common_chain|MRP-2"
CHAIN_B_SEQ = "\nMFPLVRELVCGDKTELLEPGWKNRSSIPHVTKCGQHTDFSTIPTLFLVVFSFIVIYQLLHTRTAQLRCFSPISFRIILGCLLVVDLIATVIYDLYLYVSQSPNFEVVHFYGDLVQFGGICLALILTVACKNKGIITSGVITLYWLLVVVCGIPEFRFYLSGFIYNEYALEGIRATLYIIAFTFSALELFLCCFADVPSDMYKSESSCPEYTASFINRLTFQWFTGLAYLGNKKSLENEDLWDLNEIDKAENLIPSFMQNLKPRIDEYHQNIKKDPSAALPKNHPSFVIPIFKTYKYTLLAGFFYKLCFDMLQFLAPQLLKQLIGFIEDKNQPVWIGCSIVGIMFFSSFLQSMFLHQYYHSMFRLGMHVRSVLTSAVYSKALNLSNEARKGKTIGAIVNLMSVDIQKIQDMAPTIMLFWSAPLQIFLSIYFLWKFLGVAALAGLVVLILALPVNGLIAIQMRKCQTEQMKLKDERIKMMSEILNGMKVLKLYSWERSMENMVLKIRERELHILKKLSYFMAAIVFSWICAPFLASVISFVVYVYLDPENNVLTPEITFVALSLFDILRMPLAMVAMVYGEAVQCSVSNTRLKEFFAAEEMSPQTSISHGETDSAIEVENGLFSWSSDEDPTLREISFKIQKGQLVAIVGKVGSGKSSLLHALLGEMNKLSGSVQINGNIAYVPQQAWIQNMSLRNNILFNKPYDLENYEDVVKNCALKEDLANLPAGDRTEIGEKGINLSGGQKQRVSLARAVYQNPDIILLDDPLSAVDSHVGKHIFENVISSSTGCLASKTRVLVTHGLTYLKHCDQLIVLKEGTISELGTYQELLNNSGAFAEFLEEFLIEESKTRGRVASIGDGSGEVDEILRDLGQVKPGILKRLESHLSQESDKEDTSARAIEYSRDSSRRSVLLHSPRSQHEENEALLGAISEDVPAQENTQLIEKETVETGKVKFEVYIAYFQAISIPITLLFFFLYVGSSGLGILSNFYLAKLSDHAKSGNRTSSDAKMELGIYAVLGMGQSFVVLIASIILTIGVLRASRILHAGLLGNIMRSPMAFFDVTPIGRILNRIGKDIEAIDRTLPDVIRHMSMTIFNVVATLVVIMWATPWAGIAFAILSVIYFIVLRFYISTSRQLKRLESASRSPIYSHFQESIQGASSIRAFGVVDNFIKQSQQRVDDHLIAYYPSIVANRWLAVRLEMVGNLIVLSAAGAAVYFRDSPGLSAGLVGLSVSYALNITQTLNWAVRMTSELETNIVSVERIKEYTVTPTEGNNSRRLAAKSWPEKGEISIKNFSVRYRPGLDLVLHGISAHIAPSEKVGIVGRTGAGKSSLTLALFRIIEADGGSIEIDGINIANLQLEQLRSCLTIVPQDPVLFSGTMKMNLDPFSAYSDSQVWEALENAHLKPFVKSLQDGLEHKISEGGENLSVGQRQLICLARALLRKTKVLVLDEAAAAVDVETDSLIQKTIREQFKECTVLTIAHRLNTVMDSDRLLVLDKGRVAEFDSPKNLLANPDGIFYSMAKDANVV"
TABLE_NAME = "parafold_jobs_table"
TABLE_COLUMNS = [
    "file_name",
    "msa_status",
    "msa_device",
    "folding_status",
    "folding_device",
]
ARTIFACT_NAME = "parafold_jobs"


class Scheduler:
    def __init__(self, args: Namespace) -> None:
        self.args = args
        self.wandb = wandb.init(project=os.getenv("WANDB_PROJECT_NAME"))
        self.output_file_dir = os.getenv("OUTPUT_DIR_PATH")

    def add_common_chain_to_fasta(self, fasta_file_path: str):
        with open(fasta_file_path, "r") as f:
            fasta_content = f.readlines()

        if len(fasta_content) == 2:
            print("\nAdding common chain to the fasta file...\n")
            with open(fasta_file_path, "a") as f:
                f.write(CHAIN_B_DESC)
                f.write(CHAIN_B_SEQ)
        elif len(fasta_content) == 4:
            print("\n Common Chain already exists\n")
        else:
            print("\nFasta file has more than 2 sequences\n")

    def get_fasta_file_paths(self):
        # check both .fasta and .fa files
        prot_mpnn_out_files = os.path.join(os.getenv("OUTPUT_DIR_PATH"), "protein_mpnn")
        fasta_files = glob.glob(
            f"{prot_mpnn_out_files}/**/*.fasta", recursive=True
        ) + glob.glob(f"{prot_mpnn_out_files}/**/*.fa", recursive=True)
        return fasta_files

    @staticmethod
    def obtain_output_file_path(fasta_file_dir: str):
        fasta_dirname = os.path.dirname(fasta_file_dir)
        # TODO: Remove outmulti once sto is done
        out_path = fasta_dirname.replace("protein_mpnn/", "parafold/")
        return out_path

    def submit_job_lilibet(self, fasta_file_dir, job_type="gpu"):
        fasta_file_name = os.path.basename(fasta_file_dir)
        final_output_dir = self.obtain_output_file_path(fasta_file_dir)

        if job_type == "gpu":
            # Check if msas exist locally
            
            af_multimer_command = f"bash run_alphafold.sh -d {os.getenv('AFDB_PATH')} -o {final_output_dir} -m model_1_multimer_v3 -p multimer -i {fasta_file_dir} -t 2020-12-01 -Y {args.common_chain_desc} -Z {args.common_chain_msa_path}"
        else:
            af_multimer_command = f"bash run_alphafold.sh -d {os.getenv('AFDB_PATH')} -o {final_output_dir} -m model_1_multimer_v3 -p multimer -i {fasta_file_dir} -t 2020-12-01"

        commands = [
            f"cd {os.getenv('PARAFOLD_REPO_PATH')}",
            f"conda activate {os.getenv('PARAFOLD_ENV')}",
            af_multimer_command,
        ]

        try:
            subprocess.run()
            # Execute each command
            for cmd in commands:
                # Execute the command in a new shell
                process = subprocess.run(cmd, shell=True, capture_output=True)
                process.wait()
            # print(f"Job submitted successfully in tmux session with id: {job_id}.")
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
            out_file_dir = self.obtain_output_file_path(fasta_file)
            if fasta_file.endswith(".fa"):
                file_name = os.path.basename(fasta_file).rstrip(".fa")
            elif fasta_file.endswith(".fasta"):
                file_name = os.path.basename(fasta_file).rstrip(".fasta")

            # check msa_gen status
            feat_pkl_file = os.path.join(out_file_dir, f"{file_name}/features.pkl")
            uniprot_sto_file = os.path.join(
                out_file_dir, f"{file_name}/msas/A/uniprot_hits.sto"
            )
            if (
                os.path.exists(feat_pkl_file)
                and os.path.isfile(feat_pkl_file)
                and os.path.isfile(uniprot_sto_file)
            ):
                msa_status = "completed"
                msa_device = self.args.device
            else:
                msa_status = "not_started"
                msa_device = "none"

            # check GPU status
            ranking_debug_file = os.path.join(
                out_file_dir, f"{file_name}/ranking_debug.json"
            )
            if os.path.exists(ranking_debug_file) and os.path.isfile(
                ranking_debug_file
            ):
                folding_status = "completed"
                folding_device = self.args.device
            else:
                folding_status = "not_started"
                folding_device = "none"

            # file_name, msa_status, msa_device, folding_status, folding_device
            table.append(
                {
                    "file_name": os.path.basename(fasta_file),
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

    def sync_table(self):
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
                        wandb_table[wandb_job_idx][f"{j_type}_status"] != "completed"
                        and local_table[local_job_idx][f"{j_type}_status"]
                        == "completed"
                    ):
                        print(f"Making {job} as complete in wandb table...")
                        wandb_table[wandb_job_idx][f"{j_type}_status"] = local_table[
                            local_job_idx
                        ][f"{j_type}_status"]
                        wandb_table[wandb_job_idx][
                            f"{j_type}_device"
                        ] = self.args.device

            # If job not found in wandb_table, add it
            else:
                print(f"Adding new job: {job} to wandb table...")
                wandb_table.append(local_table[local_job_idx])
        print("Finished Syncing...")
        self.upload_table_wandb(wandb_table)

    def update_table(self, table: np.array, file_name: str, status: str, job_type: str):
        # update row in the table
        for row_idx, row in enumerate(table):
            if row["file_name"] == file_name:
                table[row_idx][f"{job_type}_status"] = status
                table[row_idx][f"{job_type}_device"] = self.args.device
        self.upload_table_wandb(table)

    def check_if_msas_exist(self, out_file_dir: str, file_name: str):
        # check msa_gen status
        feat_pkl_file = os.path.join(out_file_dir, f"{file_name}/features.pkl")
        uniprot_sto_file = os.path.join(
            out_file_dir, f"{file_name}/msas/A/uniprot_hits.sto"
        )
        if (
            os.path.exists(feat_pkl_file)
            and os.path.isfile(feat_pkl_file)
            and os.path.isfile(uniprot_sto_file)
        ):
            return True
        else:
            return False

    def main(self, job_type: str, max_jobs=1):
        self.sync_table()
        curr_wandb_jobs = self.retrieve_wandb_jobs()

        OUTPUT_DIR_PATH = os.getenv("OUTPUT_DIR_PATH")
        fasta_files_local = self.get_fasta_file_paths()
        fasta_file_names = [file.split("/")[-1] for file in fasta_files_local]

        for job_idx in range(max_jobs):
            print(f"\nJob Idx: {job_idx}: Checking for jobs to submit...\n")
            for job in curr_wandb_jobs:
                # folding job
                if (
                    job_type == "folding"
                    and job["msa_status"] == "completed"  # completed msa generation
                    and job["folding_status"] == "not_started"  # not started folding
                ) or (
                    job_type == "msa"
                    and job["msa_status"] == "not_started"  # incomplete msa generation
                    and job["folding_status"] == "not_started"  # not started folding
                ):
                    self.update_table(
                        table=curr_wandb_jobs,
                        file_name=job["file_name"],
                        status="in_progress",
                        job_type=job_type,
                    )
                    full_inp_path = fasta_files_local[
                        fasta_file_names.index(job["file_name"])
                    ]
                    if not self.check_if_msas_exist(
                        self.obtain_output_file_path(full_inp_path),
                        job["file_name"],
                    ):
                        self.sync_msas(
                            file_name=job["file_name"],
                            from_device=job["msa_device"],
                            to_device=self.args.device,
                        )
                    # TODO: Submit Job
                    print(f"\nJob Idx: {job_idx}: Submitting {job[0]} \n")
                    self.add_common_chain_to_fasta(full_inp_path)
                    self.submit_job_lilibet(
                        fasta_file_dir=full_inp_path,
                        job_type="gpu",
                    )
                    print(f"\nFinished Job: {job["file_name"]} \n")
                    self.sync_table()

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
        "--common_chain_msa_path",
        type=str,
        help="Path to common chain MSA. Should end with the exact folder name containing MSA files",
        required=True,
    )
    parser.add_argument(
        "--common_chain_desc",
        type=str,
        help="Description of the common chain",
        default="common_chain",
    )

    args = parser.parse_args()
    scheduler = Scheduler(args)
