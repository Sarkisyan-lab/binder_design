import os
import logging
import psycopg2

ORDERED_TABLE_COLUMNS = [
    "id",
    "msa_status",
    "msa_device",
    "folding_status",
    "folding_device",
    "seq_len",
]

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


class RetoolDB:
    def __init__(self, jobs_table_name="af2_jobs"):
        self.jobs_table_name = jobs_table_name
        self.logger = logging.getLogger("retool_db")
        self.logger.setLevel(logging.INFO)

        conn_params = {
            "dbname": os.getenv("RETOOL_DB_NAME"),
            "user": os.getenv("RETOOL_DB_USER"),
            "password": os.getenv("RETOOL_DB_PASSWORD"),
            "host": os.getenv("RETOOL_DB_HOST"),
            "sslmode": "require",
        }
        self.conn = psycopg2.connect(**conn_params)
        self.cur = self.conn.cursor()

    def create_jobs_table(self):
        self.logger.info(f"Creating jobs table {self.jobs_table_name}")
        self.cur.execute(
            f"""
            CREATE TABLE IF NOT EXISTS {self.jobs_table_name} (
                id VARCHAR(100) PRIMARY KEY,
                msa_status VARCHAR(20) NOT NULL,
                msa_device VARCHAR(20) NOT NULL,
                folding_status VARCHAR(20) NOT NULL,
                folding_device VARCHAR(20) NOT NULL,
                seq_len INT NOT NULL
            );
        """
        )

    def delete_jobs_table(self):
        self.logger.info(f"Deleting jobs table {self.jobs_table_name}")
        self.cur.execute(f"DROP TABLE IF EXISTS {self.jobs_table_name};")
        self.conn.commit()

    def insert_jobs(self, jobs: list[dict]):
        assert len(jobs) > 0, "No jobs to insert"

        # Order the jobs based on ORDERED_TABLE_COLUMNS
        ordered_jobs = []
        for job in jobs:
            assert len(job) == len(
                ORDERED_TABLE_COLUMNS
            ), f"{job} does not have the correct number of fields"
            ordered_job = tuple(job[col] for col in ORDERED_TABLE_COLUMNS)
            ordered_jobs.append(ordered_job)

        args_str = ",".join(
            self.cur.mogrify(
                f"({','.join(['%s'] * len(ORDERED_TABLE_COLUMNS))})", x
            ).decode("utf-8")
            for x in ordered_jobs
        )
        try:
            self.cur.execute(
                f"INSERT INTO af2_jobs ({','.join(ORDERED_TABLE_COLUMNS)}) VALUES "
                + args_str
            )
        except psycopg2.errors.UniqueViolation:
            self.cur.execute("ROLLBACK;")
            self.logger.info(
                "One or more of the provided jobs already exist. Skipping entirely."
            )
        self.conn.commit()

    def execute_query(self, query: str):
        try:
            self.cur.execute(query)
        except Exception as e:
            self.cur.execute("ROLLBACK;")
            self.logger.error(f"Error executing query: {e}")

    def fetch_jobs(
        self,
        filter: dict = None,
        limit_items: int = None,
        order_by_seq_len: str = "none",
    ) -> list[dict]:
        """
        Fetch jobs from the database
        Args:
            filter (dict): Filter criteria
            limit_items (int): Limit the number of items to fetch
            order_by_seq_len (str): Order by sequence length. Options are "asc"/"desc"/"none".
        """
        if limit_items == 0:
            self.logger.warning("Limit items is 0. Fetching all items.")
            return []
        limit_str = f"LIMIT {limit_items}" if limit_items else ""

        # assert order_by_seq_len is valid
        assert order_by_seq_len in [
            "none",
            "asc",
            "desc",
        ], "Invalid order_by_seq_len"
        if order_by_seq_len == "none":
            order_by_str = ""
        elif order_by_seq_len == "asc":
            order_by_str = "ORDER BY seq_len ASC"
        elif order_by_seq_len == "desc":
            order_by_str = "ORDER BY seq_len DESC"
        if filter and len(filter) > 0:
            filter_str = " AND ".join([f"{k}='{v}'" for k, v in filter.items()])
            self.execute_query(
                f"SELECT * FROM {self.jobs_table_name} WHERE {filter_str} {order_by_str} {limit_str};"
            )
        else:
            self.execute_query(
                f"SELECT * FROM {self.jobs_table_name} {order_by_str} {limit_str};"
            )
        items = self.cur.fetchall()
        # Convert fetched items into dictionaries with keys from ORDERED_TABLE_COLUMNS
        result = []
        for item in items:
            job_dict = {col: item[idx] for idx, col in enumerate(ORDERED_TABLE_COLUMNS)}
            result.append(job_dict)
        return result

    def update_job_status(self, updated_job_entry: dict):
        """
        Update the job status for a given job id
        """
        assert "id" in list(updated_job_entry.keys()), "Job id is required"
        assert len(updated_job_entry) > 1, "At least one field must be updated"

        # Obtain id by popping out that key from the dict
        job_id = updated_job_entry.pop("id")
        set_str = ",".join([f"{k}='{v}'" for k, v in updated_job_entry.items()])
        self.execute_query(
            f"UPDATE {self.jobs_table_name} SET {set_str} WHERE id='{job_id}';"
        )
        self.conn.commit()

    def revert_folding_queued_to_not_started(self):
        self.execute_query(
            "UPDATE af2_jobs SET folding_status = 'not_started', folding_device = 'unassigned' WHERE folding_status = 'queued';"
        )
        self.conn.commit()
