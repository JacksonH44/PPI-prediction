import pytest
import os

import aiohttp
import pandas as pd

from core import config as cfg
from src.data.aws_utils import download_msa_from_openfold
from src.data.bio_apis import (
    filter_for_uniref30,
    find_uniprot_ids,
    get_interactors,
    get_sequence_lengths,
)
from src.data.combine_msa import extract_header_sequence_pairs
from src.data.create_protein_triplets import find_triplets
from src.data.data_filtering import filter_out_long_sequences
from src.data.data_processing import (
    chunk_input_genes,
    find_unique_genes,
    map_symbols_to_transcripts,
    parse_input_genes,
    remove_ground_truth_data,
    UndersamplingError,
)
from src.data.fasta_one_to_many import _create_files, process_file
from src.data.fetch_msas import find_canonical_transcript
from src.data.generate_positive_dataset import ppi_in_MANE
from src.data.generate_negative_dataset import (
    count_gene_symbols,
    find_interacting_proteins,
    find_subcellular_proteins,
    find_too_long_proteins,
    find_unsuitable_partners,
    get_locations,
    randomly_select_partners,
)


@pytest.mark.asyncio
async def test_find_too_long_proteins_success():
    """Test that finding proteins that are too long
    (>2,000 AAs) to combine with a gene of interest is
    correct."""
    location_df = pd.read_csv(
        "tests/test_data/subcellular_location_test.csv", sep="\t", usecols=["Gene name"]
    )
    genes = ["ZFHX3", "PREX2", "MTCP1"]
    expected_output = {
        "ZFHX3": set(
            ["PREX2", "ZFHX3", "MTCP1", "FANCE", "CRLF2", "CIC", "CHEK2", "CHD4"]
        ),
        "PREX2": set(["CHD4", "CHEK2", "CIC", "FANCE", "PREX2", "ZFHX3"]),
        "MTCP1": set(["CHD4", "CIC", "ZFHX3"]),
    }
    actual_output = await find_too_long_proteins(genes, location_df)
    assert actual_output == expected_output


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "input, expected_output",
    [(["TP53*PTK2", "WAS*FNBP1", "BRCA1*RHEB"], ["TP53*PTK2", "WAS*FNBP1"])],
)
async def test_filter_out_long_sequences_success(input, expected_output):
    """Test the ability to filter out PPIs that have
    total amino acid length > 2,000."""
    async with aiohttp.ClientSession() as session:
        actual_output = await filter_out_long_sequences(session, input)
    assert actual_output == expected_output


def test_map_symbol_to_transcript_success():
    """Test the functionality of mapping a symbol to its
    canonical Ensembl transcript"""
    expected_result = {"CYB5RL": "ENST00000534324", "GABARAP": "ENST00000302386"}
    actual_result = map_symbols_to_transcripts(["CYB5RL", "GABARAP"])
    assert expected_result == actual_result


@pytest.mark.asyncio
async def test_get_sequence_length_success():
    """Test that querying the Ensembl API for transcript
    sequence length is correct."""
    expected_result = {"ENST00000534324": 315}
    async with aiohttp.ClientSession() as session:
        actual_result = await get_sequence_lengths(session, ["ENST00000534324"])
    assert expected_result == actual_result


def test_extract_header_sequence_pairs_success():
    """Test for proper return of header sequence pairs
    from an MSA file."""
    expected_output = {
        "ENST00000546591|RPL41": "MRAKWRKKRMRRLKRKRRKMRQRSK",
        "UniRef100_C8Z6C0": "MRAKWRKKKTRRLKRKRRKVRARSK",
    }
    actual_output = extract_header_sequence_pairs("tests/test_data/test_msa.a3m")
    assert actual_output == expected_output


@pytest.mark.parametrize(
    "uniprot_id, expected_result", [("A0A067CGT3", True), ("Q9HB96", False)]
)
def test_download_msa_from_openfold(uniprot_id, expected_result):
    """A test for proper recognition of download files
    (sets the download flag to false)."""
    actual_result = download_msa_from_openfold(uniprot_id, "dummy/", download=False)
    assert actual_result == expected_result


@pytest.mark.parametrize(
    "input, expected_result",
    [
        (
            ["A1BG", "NAT1", "CSTF3", "ITIH5", "SLC9A7", "POLR2A"],
            ["A1BG", "NAT1", "CSTF3", "SLC9A7"],
        ),
        (["A1BG*BRCA1", "ITIH5*BRCA1", "BRCA1*ITIH5"], ["A1BG*BRCA1"]),
    ],
)
def test_filter_for_uniref30(input, expected_result):
    """Filter a list of gene names or PPIs for only those
    that have proteins that are in the UniRef30 DB."""
    actual_result = filter_for_uniref30(input)
    assert set(actual_result) == set(expected_result)


@pytest.mark.parametrize(
    "input, expected_result",
    [
        (
            ["ENST00000263100", "ENST00000307719", "ENST00000323959"],
            {
                "ENST00000263100": "P04217",
                "ENST00000307719": "P18440",
                "ENST00000323959": "Q12996",
            },
        ),
        (["ENST00000643490", "ENST00000397146", "ENST00000616978"], {}),
    ],
)
def test_find_uniprot_ids_success(input, expected_result):
    """Find UniProt SwissProt IDs from the Ensembl transcript
    version."""
    actual_result = find_uniprot_ids(input)
    assert expected_result == actual_result


def test_find_canonical_transcripts_success():
    """Test the ability to find canonical transcripts from a set of genes."""
    expected_result = {
        "A1BG": "ENST00000263100.8",
        "NAT1": "ENST00000307719.9",
        "CSTF3": "ENST00000323959.9",
    }
    actual_result = find_canonical_transcript(["A1BG", "NAT1", "CSTF3"])
    assert actual_result == expected_result


@pytest.mark.parametrize(
    "ppi, expected_result", [("MARCH7*BRCA1", False), ("FANCE*BRCA1", True)]
)
def test_ppi_in_MANE_success(ppi, expected_result):
    """Test the ppi filtering function for proteins that are found
    in the MANE summary."""
    mane_df = pd.read_csv(cfg.MANE_FILE, sep="\t", usecols=["symbol"])
    actual_result = ppi_in_MANE(ppi, mane_df)
    assert actual_result == expected_result


def test_find_unique_genes_success():
    """Test the finding of unique genes from a dataset."""
    # NOTE: We are testing with two of the same file,
    # because that is sufficient for our testing purposes
    actual_result = find_unique_genes(
        [
            "tests/test_data/positive_ppis_test.csv",
            "tests/test_data/positive_ppis_test.csv",
        ]
    )
    expected_result = {
        "ABI1",
        "ABL1",
        "AR",
        "CACNA1A",
        "ENAH",
        "EPS8",
        "NCKAP1",
        "ABI2",
        "COPS2",
        "CREBBP",
        "CTNNB1",
    }
    assert expected_result == actual_result


def test_randomly_select_partners_error():
    """Test the handling of a set of available partners with a
    size less than the number of samples."""
    with pytest.raises(UndersamplingError):
        randomly_select_partners({"A", "B", "C"}, 4)


def test_count_gene_symbols_error():
    """Test that the function properly throws a
    FileNotFoundError when the positive dataset
    hasn't been generated yet."""
    with pytest.raises(FileNotFoundError):
        count_gene_symbols("dummy.txt")


@pytest.mark.parametrize(
    "positive_ppis, expected_result",
    [("tests/test_data/positive_ppis_test.csv", {"ABI1": 4, "ABL1": 2, "AR": 3})],
)
def test_count_gene_symbols_success(positive_ppis, expected_result):
    """Test the happy path for counting gene symbols
    from the positive dataset."""
    actual_result = count_gene_symbols(positive_ppis)
    assert actual_result == expected_result


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "genes, locations_data, all_ppis, expected_result",
    [
        (
            ["BRCA1", "FANCE", "HLF"],
            {
                "Gene name": ["BRCA1", "FANCE", "HLF", "SMURF1", "KIT"],
                "Reliability": [
                    "Supported",
                    "Supported",
                    "Enhanced",
                    "Approved",
                    "Supported",
                ],
                "Main location": [
                    "Nuclear bodies;Nucleoplasm",
                    "Nucleoplasm",
                    "Nucleoplasm",
                    "Vesicles",
                    "Plasma membrane",
                ],
            },
            ["BRCA1*SMURF1", "KIT*FANCE", "BRCA1*HLF", "HLF*HLF"],
            {
                "BRCA1": {"BRCA1", "FANCE", "HLF", "SMURF1", "KIT"},
                "FANCE": {"FANCE", "BRCA1", "HLF", "KIT"},
                "HLF": {"HLF", "BRCA1", "FANCE"},
            },
        )
    ],
)
async def test_find_unsuitable_partners_success(
    genes, locations_data, all_ppis, expected_result
):
    """Test the functionality of finding all unsuitable partners for a list of
    genes"""
    locations_df = pd.DataFrame(data=locations_data)
    actual_result = await find_unsuitable_partners(genes, locations_df, all_ppis)
    assert actual_result == expected_result


@pytest.mark.parametrize(
    "location_data, expected_result",
    [
        (
            {
                "Gene name": ["CHD4", "CHEK2", "CIC", "CRLF2"],
                "Reliability": ["Approved", "Enhanced", "Approved", "Approved"],
                "Main location": [
                    "Nucleoplasm",
                    "Nucleoplasm",
                    "Nucleoplasm",
                    "Plasma membrane",
                ],
            },
            {
                "CHD4": {"CHEK2", "CIC", "CHD4"},
                "CHEK2": {"CHD4", "CIC", "CHEK2"},
                "CIC": {"CHEK2", "CHD4", "CIC"},
                "CRLF2": {"CRLF2"},
            },
        ),
        (
            {
                "Gene name": ["MDN1", "SASH1", "CILK1", "HINT3"],
                "Reliability": ["Approved", "Approved", "Approved", "Enhanced"],
                "Main location": [
                    "Cytosol;Nucleoli",
                    "Cytosol;Nucleoplasm",
                    "Nucleoli fibrillar center",
                    "Mitochondria;Nucleoli",
                ],
            },
            {
                "MDN1": {"SASH1", "HINT3", "MDN1"},
                "SASH1": {"MDN1", "SASH1"},
                "CILK1": {"CILK1"},
                "HINT3": {"MDN1", "HINT3"},
            },
        ),
    ],
)
def test_find_subcellular_proteins(location_data, expected_result):
    """Test the ability to map a dataframe of proteins and their main locations to
    a hashmap of protein : proteins with same location pairs."""
    locations_df = pd.DataFrame(data=location_data)
    actual_result = find_subcellular_proteins(locations_df)
    assert actual_result == expected_result


@pytest.mark.parametrize(
    "ppis", [(["HLF*FANCE", "HLF*BRCA1", "BRCA1*FANCE", "HOOK3*HOOK3"])]
)
def test_find_interacting_proteins_success(ppis):
    """Test the function that maps a list of protein-protein interaction
    pairs to a hashmap of protein : partners pairs."""
    expected_result = {
        "HLF": {"FANCE", "BRCA1"},
        "FANCE": {"HLF", "BRCA1"},
        "BRCA1": {"HLF", "FANCE"},
        "HOOK3": {"HOOK3"},
    }
    actual_result = find_interacting_proteins(ppis)
    assert actual_result == expected_result


def test_get_locations_success():
    """Test the function that gets subcellular locations for all genes."""
    expected_data = {
        "Gene name": ["CHEK2", "CIC", "CRLF2"],
        "Reliability": ["Supported", "Supported", "Approved"],
        "Main location": [
            "Nucleoplasm",
            "Nucleoplasm",
            "Plasma membrane",
        ],
    }
    expected_result = pd.DataFrame(data=expected_data)
    actual_result = get_locations(
        "tests/test_data/MANE_summary_v3_test.csv",
        "tests/test_data/subcellular_location_test.csv",
    )
    assert expected_result.equals(actual_result)


def test_remove_ground_truth_data_success():
    """Test the pruning of data in a dataset that
    also appears in a ground truth dataset."""
    unpruned_genes = ["NAT2", "ABCA3", "FANCE", "SERPINA3", "BRCA1"]
    expected_output = ["FANCE", "BRCA1"]
    actual_output = remove_ground_truth_data(
        unpruned_genes,
        "tests/test_data/isoform_sequences_test.xlsx",
        "1A-Gene List",
        "Gene_Symbol",
    )
    assert actual_output == expected_output


def test_chunk_input_genes_success():
    """Test the chunking of a list of genes into chunks of
    genes (makes it easier for API calls)."""
    input_genes = ["FANCE", "BRCA1", "ARID1A", "BCL10", "ERBB4"]
    expected_chunked_genes = [["FANCE", "BRCA1"], ["ARID1A", "BCL10"], ["ERBB4"]]
    actual_chunked_genes = chunk_input_genes(input_genes, 2)
    assert actual_chunked_genes == expected_chunked_genes


@pytest.mark.parametrize(
    "infile", [("tests/test_data/cancer_driver_gene_list_test.csv")]
)
def test_parse_input_genes_success(infile):
    """Test the generation of reference gene list for PPI pairs."""
    expected_output = [
        "CHD4",
        "CHEK2",
        "CIC",
        "CIITA",
        "CLIP1",
        "CLTC",
        "CLTCL1",
        "CNBP",
        "CNOT3",
        "CNTRL",
        "COL1A1",
        "COL2A1",
        "CREB1",
        "CREB3L1",
        "CREB3L2",
        "CREBBP",
        "CRLF2",
        "CRTC1",
    ]
    actual_output = parse_input_genes(infile)
    assert expected_output == actual_output


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "gene_list, expected_interactions, threshold_level, relax_evidence",
    [
        (
            ["FANCE"],
            {
                "FANCE*FANCC",
                "FANCE*FANCA",
                "FANCE*FANCD2",
                "FANCE*FANCF",
                "FANCE*FANCM",
                "FANCE*FANCG",
                "FANCE*HES1",
                "FANCE*APITD1",
                "FANCE*STRA13",
            },
            2,
            False,
        ),
        (
            ["ASXL1", "SS18"],
            {
                "ASXL1*BAP1",
                "ASXL1*FOXK1",
                "ASXL1*FOXK2",
                "ASXL1*HCFC1",
                "ASXL1*HIST1H1C",
                "ASXL1*AKT1",
                "SS18*SMARCA2",
                "SS18*SMARCB1",
                "SS18*RBM14",
            },
            3,
            False,
        ),
        (["TMPRSS2"], set(), 2, False),
        (
            ["HLF"],
            {"HLF*CREBBP", "HLF*TET2", "HLF*TLK1", "HLF*DBP", "HLF*HNF4G", "HLF*MYB"},
            10,
            True,
        ),
    ],
)
async def test_get_interactors_success(
    gene_list, expected_interactions, threshold_level, relax_evidence
):
    """Test creation of PPI interactions."""
    async with aiohttp.ClientSession() as session:
        ppis = await get_interactors(
            session, gene_list, threshold_level, relax_evidence
        )
    actual_interactions = set(ppis)
    assert actual_interactions == expected_interactions


@pytest.mark.parametrize(
    "test_file_path, expected_data, positive",
    [
        (
            "tests/test_data/ppis_test1.xlsx",
            {
                "ref_ID": [
                    "ACTN4_1",
                    "ACTN4_1",
                    "AKT1_1",
                    "AKT1_1",
                    "AKT1_1",
                    "BAG1_1",
                ],
                "alt_ID": [
                    "ACTN4_4",
                    "ACTN4_4",
                    "AKT1_2",
                    "AKT1_2",
                    "AKT1_2",
                    "BAG1_2",
                ],
                "bait_ID": ["TRIM23", "MYOZ2", "TCL1A", "TMCC2", "MTUS2", "HSPA8"],
                "perturbation": [True, True, False, True, True, True],
            },
            False,
        ),
        (
            "tests/test_data/ppis_test2.xlsx",
            {
                "ref_ID": [
                    "BCL2L1_1",
                    "BCL2L1_1",
                    "BCL2L1_1",
                    "BCL2L1_1",
                    "BCL2L1_1",
                    "BCL2L1_1",
                    "BTC_1",
                ],
                "alt_ID": [
                    "BCL2L1_2",
                    "BCL2L1_2",
                    "BCL2L1_2",
                    "BCL2L1_3",
                    "BCL2L1_3",
                    "BCL2L1_3",
                    "BTC_2",
                ],
                "bait_ID": ["BAD", "BIK", "BMF", "BAD", "BIK", "BMF", "GMPPA"],
                "perturbation": [True, True, True, True, False, False, True],
            },
            True,
        ),
        (
            "tests/test_data/ppis_test3.xlsx",
            {
                "ref_ID": ["CLCN2_1", "CLCN2_1", "CLCN2_1", "CLCN2_1", "CLCN2_1"],
                "alt_ID": ["CLCN2_2", "CLCN2_3", "CLCN2_4", "CLCN2_4", "CLCN2_5"],
                "bait_ID": ["UBQLN1_ORF2", "FHL3", "FHL3", "UBQLN1_ORF2", "FHL3"],
                "perturbation": [True, True, False, False, False],
            },
            False,
        ),
    ],
)
def test_find_triplets_success(test_file_path, expected_data, positive):
    """Test that the observation dataframe creation function works, includes the
    correct labels for each triplet."""
    expected_df = pd.DataFrame(data=expected_data)
    actual_df = find_triplets(test_file_path, positive)
    assert expected_df.equals(actual_df)


def test_process_file_error():
    """Error handling test for a file that doesn't exist."""
    with pytest.raises(FileNotFoundError):
        actual_fname, actual_seq = process_file("../data/nonexistent.txt")


def test_process_file_success():
    """Test that FASTA files are correctly separated
    into the header and the protein sequence."""
    expected_header = ">ENSG00000002822|ENST00000455998"
    expected_sequence = """MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQREVDRNQELL
TRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRIS
ELQWSVMDQEMRVKRLESEKQELQ"""
    actual_headers, actual_sequences = process_file("tests/test_data/test_gene1.txt")
    assert (
        actual_headers[1] == expected_header
        and actual_sequences[1] == expected_sequence
    )


def test__create_files_success():
    """Test the creation of individual FASTA files
    from header sequence lists."""
    outfolder_name = os.path.join(os.getcwd(), "tests", "test_data", "test_gene1")
    fname = os.path.join(outfolder_name, "ENSG00000002822_ENST00000421113.fasta")
    expected_content = """>ENSG00000002822|ENST00000421113
MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE"""
    header = [">ENSG00000002822|ENST00000421113"]
    sequence = ["MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE"]
    _create_files(header, sequence, outfolder_name)
    with open(fname, "r") as file:
        actual_content = file.read()
    assert actual_content == expected_content
