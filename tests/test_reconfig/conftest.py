import os 
import sys
import pytest
import shutil
import pathlib
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from tests.test_reconfig import benchmark_compare

output_dir = pathlib.Path.cwd() / "ExampleData" / "OrthoFinder" 

output_list = sorted([(entry.stat().st_mtime, entry.path) for entry in os.scandir(output_dir)])
latest_output = output_list[-1][1]

if os.path.exists(output_dir):
    for output in output_dir.iterdir():
        if output.name != os.path.basename(latest_output):
            shutil.rmtree(output)

# Please change this path accordingly before you run this test
# path_to_benchmark = pathlib.Path(r"/home/yiliu/benchmark/OrthoFinder/")
# bencmark_output_dir = path_to_benchmark / "ExampleData" / "OrthoFinder"

path_to_benchmark = pathlib.Path.cwd()
bencmark_output_dir = path_to_benchmark / "BENCHMARK_OrthoFinderV255" / "OrthoFinder"

benchmark_output_list = sorted([(entry.stat().st_mtime, entry.path) for entry in os.scandir(bencmark_output_dir)])

benchmark_latest_output = pathlib.Path(benchmark_output_list[-1][1])
curr_benchmark_name = benchmark_latest_output.name

if os.path.exists(bencmark_output_dir):
    for output in bencmark_output_dir.iterdir():
        if output.name != curr_benchmark_name:
            shutil.rmtree(output)

if "benchmark" not in curr_benchmark_name:
    new_benchmark_name = benchmark_latest_output.name + "_" + "benchmark"
    new_benchmark_path = benchmark_latest_output.parent / new_benchmark_name
    benchmark_latest_output.rename(new_benchmark_path)


"""
Before running this test, please make sure the test data is obtained with the `-efn`.
"""

@pytest.fixture(scope="module")
def orthogroup_stats():
    orthostats = benchmark_compare.orthogroups_analysis(pathlib.Path.cwd(), "ExampleData", path_to_benchmark)
    return orthostats

@pytest.fixture(scope="module")
def unassigned_gens_stats():
    orthostats = benchmark_compare.unassigned_gens_analysis(pathlib.Path.cwd(), "ExampleData", path_to_benchmark)
    return orthostats

@pytest.fixture(scope="module")
def orthologues_stats():
    orthostats = benchmark_compare.orthologues_analysis(pathlib.Path.cwd(), "ExampleData", path_to_benchmark)
    return orthostats