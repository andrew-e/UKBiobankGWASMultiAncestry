import pandas as pd
import logging
import sys

# logging
logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(filename)s:%(lineno)d %(message)s",
    datefmt="%Y-%m-%d:%H:%M:%S",
    level=logging.DEBUG,
)
logger = logging.getLogger(__name__)

job_cols = [
    "name",
    "application_id",
    "pheno_file",
    "pheno_col",
    "covar_file",
    "covar_col",
    "qcovar_col",
    "method",
    "ancestral_group",
    "trait_type",
]


def read_jobs(input_path, user):
    job_file = f"{input_path}/data/phenotypes/{user}/input/jobs.csv"
    logger.info(f"Reading {job_file}")

    try:
        job_df = pd.read_csv(job_file)
        logger.info(f"\n{job_df.head()}")
        # check columns
        if not list(job_df.columns) == job_cols:
            logger.error(f"File structure is not right: {list(job_df.columns)}")
            sys.exit()
        else:
            return job_df
    except:
        logger.error(f"Can't read {job_file}")
        sys.exit()


def check_data_size(row, ancestry, sample_size):
    logger.debug(row)
    logger.debug(row["trait_type"])
    logger.debug(row["name"])
    if row["trait_type"] == "continuous":
        logger.info("Continuous trait")
        logger.info(f"Ancestry: {ancestry}, Sample size: {sample_size}")
        if sample_size <= 300:
            logger.warning(f"{ancestry} sample size is too small, it is recommened by guildelines to have at least 300 samples")
    elif row["trait_type"] == "categorical":
        logger.info("Categorical trait")
        if 'get_case_size' in globals():
            case_size = get_case_size(row["name"]) # get case size by phenotype name
        else:
            case_size = -999
            logger.warning("get_case_size function not found, the case size is set to -999")
        logger.info(f"Ancestry: {ancestry}, Case size: {case_size}")
        if case_size <= 100:
            logger.warning(f"Case size for {ancestry} is too small, it is recommened by guildelines to have at least 100 cases")
