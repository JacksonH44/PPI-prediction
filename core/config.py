"""
Load the config file and create custom variables
that are available for ease of use.
"""


import os

import yaml

BASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")

# Check if config file exists
if os.path.exists('config/config.yml'):
    with open(BASE_DIR + "/config/config.yml", "r") as configFile:
        data = configFile.read()

    data = yaml.load(data, Loader=yaml.FullLoader)

    # Specify all sensitive data from the config file
    NCBI_API_KEY = data["ncbi"]["api_key"]
    BIOGRID_API_KEY = data["biogrid"]["api_key"]

    AWS_ID = data["aws"]["access_key_id"]
    AWS_SECRET_KEY = data["aws"]["secret_access_key"]
else:
    with open(BASE_DIR + "/config/config.template.yml", "r") as configFile:
        data = configFile.read()

    data = yaml.load(data, Loader=yaml.FullLoader)

    # Specify all sensitive data from the environment variables (e.g., GitHub Actions Secrets)
    NCBI_API_KEY = os.environ["NCBI_API_KEY"]
    BIOGRID_API_KEY = os.environ["BIOGRID_API_KEY"]

    AWS_ID = os.environ["AWS_ID"]
    AWS_SECRET_KEY = os.environ["AWS_SECRET_KEY"]

NCBI_BASE_URL = data["ncbi"]["base_url"]

BIOGRID_BASE_URL = data["biogrid"]["base_url"]
BIOGRID_STRICT_EVIDENCE = data["biogrid"]["strict_evidence"]
BIOGRID_RELAXED_EVIDENCE = data["biogrid"]["relaxed_evidence"]

ENSEMBL_BASE_URL = data["ensembl"]["base_url"]

GROUND_TRUTH_PATH = data["ground_truth"]["file_path"]
GROUND_TRUTH_SHEET = data["ground_truth"]["sheet_name"]
GROUND_TRUTH_COLUMN = data["ground_truth"]["column"]

SEED = data["random"]["seed"]

MANE_FILE = data["reference_files"]["mane"]

BIOMART_BASE_URL = data["biomart"]["base_url"]

AWS_REGION = data["aws"]["region"]
