import os
import sys
from argparse import Namespace, ArgumentParser
import logging
import firebase_admin
from firebase_admin import db
import shutil
from dotenv import load_dotenv

sys.path.append(os.path.abspath(os.path.join(__file__, "../")))
import utils

# load dotenv
load_dotenv(os.path.abspath(os.path.join(__file__, "../../.env")))
import subprocess



# constants
STATUS_UNASSIGNED = "unassigned"
STATUS_NOT_STARTED = "not_started"
STATUS_RUNNING = "running"
STATUS_COMPLETED = "completed"


def main(args: Namespace):
    """
    Wrapper for the colabfold batch command.

    Args:
        args: Namespace
            - msa_file_path: str
            - predictions_dir: str
            - use-dropout: bool
            - num-recycle: int
            - num-models: int
            - lilibet_output_dir: str
    """
    logger = logging.getLogger("wrapper")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s :: %(levelname)s :: %(message)s",
    )

    # initialize firebase
    cred_obj = firebase_admin.credentials.Certificate(
        os.path.abspath(os.path.join(__file__, "../src/login_key.json"))
    )
    firebase_admin.initialize_app(
        cred_obj, {"databaseURL": os.environ["FIREBASE_DB_URL"]}
    )
    ref = db.reference("/")

    # Set the task status to running
    ref.child(args.task_id).child("folding_status").set(STATUS_RUNNING)

    if args.use_dropout:
        use_dropout = "--use-dropout"
    else:
        use_dropout = ""

    colabfold_cmd = f"colabfold_batch {args.msa_file_path} {args.predictions_dir} --num-recycle {args.num_recycle} --num-models {args.num_models} {use_dropout}"

    # Run the colabfold command
    subprocess.run(["bash", "-c", colabfold_cmd])

    # check if predictions are complete
    file_name = os.path.basename(args.msa_file_path)
    if not utils.check_if_predictions_complete(args.predictions_dir):
        logger.warning("PREDICTIONS FAILED!")
        logger.warning(f"Removing predictions directory: {args.predictions_dir}")
        shutil.rmtree(args.predictions_dir)
        logger.warning(
            f"Setting folding status for {file_name} to {STATUS_NOT_STARTED}"
        )
        ref.child(args.task_id).child("folding_status").set(STATUS_NOT_STARTED)

    else:
        logger.info(f"Setting folding status for {file_name} to {STATUS_COMPLETED}")
        ref.child(args.task_id).child("folding_status").set(STATUS_COMPLETED)

        if args.lilibet_output_dir != "":
            lilibet_host = os.environ["LILIBET_HOST"]
            lilibet_port = os.environ["LILIBET_PORT"]
            lilibet_output_dir = args.lilibet_output_dir

            scp_cmd = f"scp -r -P {lilibet_port} {args.predictions_dir} {lilibet_host}:{lilibet_output_dir}/predictions/"
            subprocess.run(["bash", "-c", scp_cmd])
            logging.debug(f"Copied predictions to lilibet: {file_name}")

    logging.info(f"Finished folding: {file_name}!")


if __name__ == "__main__":
    argparse_bool = lambda x: (str(x).lower() == "true")
    parser = ArgumentParser(description="Arguments for the colabfold wrapper")

    parser.add_argument(
        "--task_id",
        type=str,
        required=True,
        help="ID of the task.",
    )

    parser.add_argument(
        "--msa_file_path",
        type=str,
        required=True,
        help="Path to the MSA file.",
    )

    parser.add_argument(
        "--predictions_dir",
        type=str,
        required=True,
        help="Path to the directory where predictions will be saved.",
    )

    parser.add_argument(
        "--use-dropout",
        type=argparse_bool,
        default=False,
        help="If set, will use dropout in the colabfold model.",
    )

    parser.add_argument(
        "--num-recycle",
        type=int,
        default=3,
        help="Number of recycling steps.",
    )

    parser.add_argument(
        "--num-models",
        type=int,
        default=5,
        help="Number of models to use for prediction.",
    )

    parser.add_argument(
        "--lilibet_output_dir",
        type=str,
        default="",
        help="If set, will copy the predictions to lilibet.",
    )

    args = parser.parse_args()

    main(args)
