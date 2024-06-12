"""
Load the config file and create custom variables
that are available for ease of use.
"""


import os

import yaml

BASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")


with open(BASE_DIR + "/config/config.yml", "r") as configFile:
    data = configFile.read()

data = yaml.load(data, Loader=yaml.FullLoader)

# Custom variables
NCBI_API_KEY = data["ncbi"]["api_key"]
NCBI_BASE_URL = data["ncbi"]["base_url"]

BIOGRID_API_KEY = data["biogrid"]["api_key"]
BIOGRID_BASE_URL = data["biogrid"]["base_url"]
BIOGRID_STRICT_EVIDENCE = data["biogrid"]["strict_evidence"]
BIOGRID_RELAXED_EVIDENCE = data["biogrid"]["relaxed_evidence"]

GROUND_TRUTH_PATH = data["ground_truth"]["file_path"]
GROUND_TRUTH_SHEET = data["ground_truth"]["sheet_name"]
GROUND_TRUTH_COLUMN = data["ground_truth"]["column"]
TRIPLET_FILE = data["ground_truth"]["triplet"]
