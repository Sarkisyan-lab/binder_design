import os
import sys
from argparse import Namespace, ArgumentParser
import logging
import firebase_admin
from firebase_admin import db
import shutil
from dotenv import load_dotenv
import subprocess

sys.path.append(os.path.abspath(os.path.join(__file__, "../")))
import utils
from retool_db import (
    RetoolDB,
    STATUS_UNASSIGNED,
    STATUS_NOT_STARTED,
    STATUS_RUNNING,
    STATUS_COMPLETED,
)

# load dotenv
load_dotenv(os.path.abspath(os.path.join(__file__, "../../.env")))


def main(args: Namespace):
    """
    Wrapper for the colabfold batch command.

    Args:
        args: Namespace
            - task_id: str
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

    # initialize RetoolDB
    db = RetoolDB()

    # Set the task status to running
    db.update_job_status(
        {
            "id": args.task_id,
            "folding_status": STATUS_RUNNING,
        }
    )

    if args.use_dropout:
        use_dropout = "--use-dropout"
    else:
        use_dropout = ""

    colabfold_cmd = f"colabfold_batch {args.msa_file_path} {args.predictions_dir} --num-recycle {args.num_recycle} --num-models {args.num_models} {use_dropout}"

    # Run the colabfold command
    subprocess.run(["bash", "-c", colabfold_cmd])

    # check if predictions are complete
    if not utils.check_if_predictions_complete(args.predictions_dir):
        logger.warning("PREDICTIONS FAILED!")
        logger.warning(f"Removing predictions directory: {args.predictions_dir}")
        shutil.rmtree(args.predictions_dir)
        logger.warning(
            f"Setting folding status for {args.task_id} to {STATUS_NOT_STARTED}"
        )
        db.update_job_status(
            {
                "id": args.task_id,
                "folding_status": STATUS_NOT_STARTED,
                "folding_device": STATUS_UNASSIGNED,
            }
        )

    else:
        logger.info(f"Setting folding status for {args.task_id} to {STATUS_COMPLETED}")
        db.update_job_status(
            {
                "id": args.task_id,
                "folding_status": STATUS_COMPLETED,
            }
        )

    logger.info(f"Finished folding: {args.task_id}!")


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
