"""
This file contains utility functions for the Colabfold project.
"""

import os
import shutil
from dotenv import load_dotenv
import logging
from Bio.PDB import PDBParser, MMCIFParser
import numpy as np
import glob
logger = logging.getLogger("utils")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s :: %(levelname)s :: %(message)s",
)

# load environment variables
load_dotenv(os.path.abspath(os.path.join(__file__, "../../.env")))


def fasta_file_to_fasta_str(fasta_file_path: str) -> str:
    """
    Converts a fasta file to a fasta string

    Args:
        fasta_file_path (str): Path to the fasta file

    Returns:
        str: Fasta string
    """
    with open(fasta_file_path, "r") as f:
        fasta_lines = f.readlines()
    fasta_content = fasta_lines[0] + "".join(fasta_lines[1:]).replace("\n", "")
    return fasta_content


def save_fasta_from_fasta_str(fasta_str: str, file_name: str) -> str:
    """
    Saves a fasta string to a fasta file

    Args:
        fasta_str (str): Fasta string
        file_name (str): Path to the fasta file
    """
    if not os.path.exists(os.path.dirname(file_name)):
        os.makedirs(os.path.dirname(file_name), exist_ok=True)

    with open(file_name, "w") as f:
        for line in fasta_str.split():
            f.write(line + "\n")


def obtain_device_name() -> str:
    """
    Obtains the device name

    Returns:
        str: Device name
    """
    file_path = os.path.abspath(__file__)
    file_path_head = file_path.split("/")[1:3]

    jex_path_head = os.getenv("JEX_COLABFOLD_CONDA_ENV").split("/")[1:3]
    lilibet_path_head = os.getenv("LILIBET_COLABFOLD_CONDA_ENV").split("/")[1:3]
    hpc_cx3_path_head = os.getenv("HPC_CX3_COLABFOLD_CONDA_ENV").split("/")[1:3]
    hpc_hx1_path_head = os.getenv("HPC_HX1_COLABFOLD_CONDA_ENV").split("/")[1:3]

    if file_path_head == jex_path_head:
        return "jex"
    elif file_path_head == lilibet_path_head:
        return "lilibet"
    elif file_path_head == hpc_cx3_path_head:
        return "hpc_cx3"
    elif file_path_head == hpc_hx1_path_head:
        return "hpc_hx1"
    else:
        raise ValueError("Device not found")

def obtain_colabfold_conda_env_path(device: str = None) -> str:
    """
    Obtains the conda environment path for colabfold

    Returns:
        str: Path to the conda environment
    """
    if not device:
        device = obtain_device_name()
    
    if device == "jex":
        return os.getenv("JEX_COLABFOLD_CONDA_ENV")
    elif device == "lilibet":
        return os.getenv("LILIBET_COLABFOLD_CONDA_ENV")
    elif device == "hpc_cx3":
        return os.getenv("HPC_CX3_COLABFOLD_CONDA_ENV")
    elif device == "hpc_hx1":
        return os.getenv("HPC_HX1_COLABFOLD_CONDA_ENV")
    else:
        raise ValueError("Device not found")

def check_if_predictions_complete(predictions_folder: str) -> bool:
    """
    Checks if the predictions are complete

    Args:
        predictions_folder (str): Path to the predictions folder
        file_name (str): Name of the file to check

    Returns:
        bool: True if the predictions are complete, False otherwise
    """

    done_files = [f for f in os.listdir(predictions_folder) if f.endswith(".done.txt")]
    return len(done_files) > 0

def check_if_fasta_file_exists(fasta_file_path: str) -> bool:
    """
    Checks if the fasta file exists

    Args:
        fasta_file_path (str): Path to the fasta file

    Returns:
        bool: True if the fasta file exists, False otherwise
    """
    return os.path.exists(fasta_file_path)

def check_if_msas_exist(msa_folder_path: str, file_name: str) -> bool:
    """
    Checks if the msas exist

    Args:
        msas_folder (str): Path to the msas folder

    Returns:
        bool: True if the msas exist, False otherwise
    """
    if os.path.exists(msa_folder_path) and len(os.listdir(msa_folder_path)) > 0:
        msas_list = os.listdir(msa_folder_path)
        for msa_file in msas_list:
            if msa_file.endswith(".a3m") and file_name in msa_file:
                return True
    return False


def clean_predictions_folder(predictions_folder):
    """
    Cleans the predictions folder by removing all directories that do not have a .done.txt file

    Args:
        predictions_folder (str): Path to the predictions folder
    """
    # List all directories in the predictions folder
    directories = [
        d
        for d in os.listdir(predictions_folder)
        if os.path.isdir(os.path.join(predictions_folder, d))
    ]

    for directory in directories:
        dir_path = os.path.join(predictions_folder, directory)
        if not check_if_predictions_complete(dir_path):
            logger.info(f"Removing {dir_path} because the predictions are not complete")
            shutil.rmtree(dir_path)


def get_residue_by_id(chain, res_seq):
    for residue in chain:
        if residue.get_id()[1] == res_seq:
            return residue
    raise Exception(f"Residue {res_seq} not found in chain {chain.id}")


def check_if_chain_is_above_plane(
    constr_plane_chain_id: str,
    constr_plane_res_seq: list[int],
    check_chain_id: str,
    pdb_file_path: str,
) -> tuple[int, int, int]:
    """
    Checks if the check chain is above the constructed plane

    Args:
        constr_plane_chain_id (str): Chain id of the chain used to construct the plane
        constr_plane_res_seq (list[int]): The three residues used to construct the plane
        check_chain_id (str): Chain id of the chain to check
        pdb_file_path (str): Path to the pdb file
    
    Returns:
        tuple[int, int, int]: Number of residues above, below, and on the plane
    """
    # Parse the PDB file
    if pdb_file_path.endswith(".pdb"):
        parser = PDBParser()
        
    elif pdb_file_path.endswith(".cif"):
        parser = MMCIFParser()
    
    structure = parser.get_structure('protein', pdb_file_path)

    # Plane chain ID
    plane_chain_id = structure[0][constr_plane_chain_id]

    # Get the coordinates of the CA atoms of the specified residues
    coords = []

    for res_num in constr_plane_res_seq:
        residue = get_residue_by_id(plane_chain_id, res_num)
        if 'CA' not in residue:
            raise Exception(f"Residue {res_num} does not have a CA atom.")
        ca_atom = residue['CA']
        coords.append(ca_atom.get_coord())

    coords = np.array(coords)

    # Define the plane using the three points
    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[0]

    # Normal vector to the plane
    normal_vector = np.cross(v1, v2)
    normal_vector = normal_vector / np.linalg.norm(normal_vector)

    # Get chain to check
    chain_to_check = structure[0][check_chain_id]

    num_positive = 0
    num_negative = 0
    num_zero = 0

    for residue in chain_to_check:
        # Skip residues without a CA atom
        if 'CA' not in residue:
            continue
        ca_atom = residue['CA']
        coord = ca_atom.get_coord()
        # Compute the signed distance from the plane
        distance = np.dot(normal_vector, coord - coords[0])
        if distance > 0:
            num_positive += 1
        elif distance < 0:
            num_negative += 1
        else:
            num_zero += 1  # Residue lies exactly on the plane

    return num_positive, num_negative, num_zero

def check_chain_is_above_plane_for_dir(
    constr_plane_chain_id: str,
    constr_plane_res_seq: list[int],
    check_chain_id: str,
    pdb_dir_path: str,
    check_above_or_below: str = "above"
) -> list:
    """
    Checks if the check chain is above the constructed plane for a directory of PDB files
    """
    pdb_files = glob.glob(os.path.join(pdb_dir_path, "**", "*.pdb"), recursive=True) + glob.glob(os.path.join(pdb_dir_path, "**", "*.cif"), recursive=True)

    results = []
    for pdb_file in pdb_files:
        try:
            num_positive, num_negative, num_zero = check_if_chain_is_above_plane(
                constr_plane_chain_id,
                constr_plane_res_seq,
                check_chain_id,
                pdb_file
            )
            if check_above_or_below == "above" and num_positive > num_negative:
                logger.info(f"PDB file {pdb_file} has more residues above the plane")
                results.append(pdb_file)
            elif check_above_or_below == "below" and num_negative > num_positive:
                logger.info(f"PDB file {pdb_file} has more residues below the plane")
                results.append(pdb_file)
            else:
                logger.info(f"PDB file {pdb_file} has {num_positive} residues above, {num_negative} residues below.")
        except Exception as e:
            logger.error(f"Error checking {pdb_file}: {e}")
    
    return results