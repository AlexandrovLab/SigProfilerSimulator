# import os
# import filecmp
# import pytest
# from SigProfilerSimulator import SigProfilerSimulator as sigsim

# # --- Dynamic Path Setup ---
# # This script's location: e.g., .../SigProfilerSimulator/test/
# this_dir = os.path.dirname(os.path.abspath(__file__))

# # Root of the test folder
# base_test_dir = this_dir  # Already in .../SigProfilerSimulator/test/

# # Path to seed file
# seed_file = os.path.join(base_test_dir, "seed_file", "Simulator_seeds.txt")


# # --- Folder Configs ---
# folder_configs = [
#     ("seed1234-chromT", True, seed_file),
#     ("seedNone-chromT", True, None),
#     ("seed1234-chromF", False, seed_file),
#     ("seedNone-chromF", False, None),
# ]


# # --- Helper Functions ---
# def run_simulator(project_path, chrom_based, seed_file=None):
#     print(
#         f">>> Running simulator at {project_path} | chrom_based={chrom_based} | seed={'yes' if seed_file else 'no'}"
#     )
#     sigsim.SigProfilerSimulator(
#         project="test_with_seed",
#         project_path=project_path,
#         genome="GRCh37",
#         contexts=["96"],
#         chrom_based=chrom_based,
#         seed_file=seed_file,
#         seqInfo=False,
#         exome=False,
#         cushion=0,
#     )


# def get_maf_paths(folder_name):
#     path1 = os.path.join(
#         base_test_dir,
#         folder_name,
#         "example",
#         "output",
#         "simulations",
#         "test_with_seed_simulations_GRCh37_96",
#         "1.maf",
#     )
#     path2 = os.path.join(
#         base_test_dir,
#         folder_name,
#         "example_copy",
#         "output",
#         "simulations",
#         "test_with_seed_simulations_GRCh37_96",
#         "1.maf",
#     )
#     return path1, path2


# def prepare_simulation_runs():
#     for folder_name, chrom_based, seed in folder_configs:
#         for subfolder in ["example", "example_copy"]:
#             out_path = os.path.join(base_test_dir, folder_name, subfolder)
#             os.makedirs(out_path, exist_ok=True)
#             run_simulator(
#                 project_path=out_path, chrom_based=chrom_based, seed_file=seed
#             )


# # --- Tests ---
# @pytest.mark.parametrize(
#     "folder_name, should_match",
#     [
#         ("seed1234-chromT", True),
#         ("seedNone-chromT", False),
#         ("seed1234-chromF", True),
#         ("seedNone-chromF", False),
#     ],
# )
# def test_maf_file_behavior(folder_name, should_match):
#     maf1, maf2 = get_maf_paths(folder_name)

#     assert os.path.exists(maf1), f"Missing file: {maf1}"
#     assert os.path.exists(maf2), f"Missing file: {maf2}"

#     if should_match:
#         assert filecmp.cmp(
#             maf1, maf2, shallow=False
#         ), f"Expected matching files in {folder_name}, but they differ!"
#     else:
#         assert not filecmp.cmp(
#             maf1, maf2, shallow=False
#         ), f"Expected different files in {folder_name}, but they match!"


# # --- Optional: Run simulations directly before testing ---
# if __name__ == "__main__":
#     prepare_simulation_runs()


import os
import filecmp
import pytest
from SigProfilerSimulator import SigProfilerSimulator as sigsim

# --- Dynamic Path Setup ---
this_dir = os.path.dirname(os.path.abspath(__file__))
base_test_dir = this_dir
seed_file = os.path.join(base_test_dir, "seed_file", "Simulator_seeds.txt")

# --- Folder Configs ---
folder_configs = [
    ("seed1234-chromT", True, seed_file),
    ("seedNone-chromT", True, None),
    ("seed1234-chromF", False, seed_file),
    ("seedNone-chromF", False, None),
]


# --- Helper Functions ---
def run_simulator(project_path, chrom_based, seed_file=None):
    print(
        f">>> Running simulator at {project_path} | chrom_based={chrom_based} | seed={'yes' if seed_file else 'no'}"
    )
    try:
        sigsim.SigProfilerSimulator(
            project="test_with_seed",
            project_path=project_path,
            genome="GRCh37",
            contexts=["96"],
            chrom_based=chrom_based,
            seed_file=seed_file,
            seqInfo=False,
            exome=False,
            cushion=0,
        )
    except Exception as e:
        print(f"Simulation failed at {project_path} with error:\n{e}")


def get_maf_paths(folder_name):
    path1 = os.path.join(
        base_test_dir,
        folder_name,
        "example",
        "output",
        "simulations",
        "test_with_seed_simulations_GRCh37_96",
        "1.maf",
    )
    path2 = os.path.join(
        base_test_dir,
        folder_name,
        "example_copy",
        "output",
        "simulations",
        "test_with_seed_simulations_GRCh37_96",
        "1.maf",
    )
    return path1, path2


def prepare_simulation_runs():
    for folder_name, chrom_based, seed in folder_configs:
        for subfolder in ["example", "example_copy"]:
            out_path = os.path.join(base_test_dir, folder_name, subfolder)
            os.makedirs(out_path, exist_ok=True)
            run_simulator(
                project_path=out_path, chrom_based=chrom_based, seed_file=seed
            )


# --- Pytest Fixture: Automatically run simulations before all tests ---
@pytest.fixture(scope="session", autouse=True)
def setup_simulations():
    print("\n>>> Preparing simulation runs before tests")
    prepare_simulation_runs()


# --- Tests ---
@pytest.mark.parametrize(
    "folder_name, should_match",
    [
        ("seed1234-chromT", True),
        ("seedNone-chromT", False),
        ("seed1234-chromF", True),
        ("seedNone-chromF", False),
    ],
)
def test_maf_file_behavior(folder_name, should_match):
    maf1, maf2 = get_maf_paths(folder_name)

    assert os.path.exists(maf1), f"Missing file: {maf1}"
    assert os.path.exists(maf2), f"Missing file: {maf2}"

    if should_match:
        assert filecmp.cmp(
            maf1, maf2, shallow=False
        ), f"Expected matching files in {folder_name}, but they differ!"
    else:
        assert not filecmp.cmp(
            maf1, maf2, shallow=False
        ), f"Expected different files in {folder_name}, but they match!"
