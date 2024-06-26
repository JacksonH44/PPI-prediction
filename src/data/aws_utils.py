"""
A collection of AWS s3 bucket utility functions.
"""

import os
import sys

import boto3

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from core import config as cfg


def download_msa_from_openfold(uniprot_id, save_dir, download=True) -> bool:
    """Download the MSA for the input uniprot_id."""
    s3 = boto3.client(
        "s3",
        region_name=cfg.AWS_REGION,
        aws_access_key_id=cfg.AWS_ID,
        aws_secret_access_key=cfg.AWS_SECRET_KEY,
    )
    results = s3.list_objects(
        Bucket="openfold", Prefix=f"uniclust30/{uniprot_id}/a3m/uniclust30.a3m"
    )
    if "Contents" in results:
        for obj in results["Contents"]:
            key = obj["Key"]
            if key == f"uniclust30/{uniprot_id}/a3m/uniclust30.a3m":
                if download:
                    os.makedirs(save_dir, exist_ok=True)
                    full_path = os.path.join(save_dir, f"{uniprot_id}.a3m")
                    s3.download_file("openfold", key, full_path)
                return True
    return False
