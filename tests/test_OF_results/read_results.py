import os
import re
from collections import defaultdict
# import glob 
from ete4 import Tree


node_pattern = re.compile(r"^n(\d+)$", re.IGNORECASE)

TEST_CATEGORIES = [
    "core",
    "assign"
]

TEST_ITEMS = [
    "statistics_overall",
    "duplications",
    "orthogroups",
    "resolved_gene_trees",
    "orthologues",
    "hogs"
]

EXPECTED_RESULTS = {
    test_cat: {
        item: {}
        for item in TEST_ITEMS
    }
    for test_cat in TEST_CATEGORIES
}

def get_stats_info(file_path):
    overall_stats_info = {}
    og_nspecies = {}
    in_og_table = False

    with open(file_path, "r", encoding="utf-8") as reader:
        for raw in reader:
            line = raw.strip()
            if not line or line.startswith("#"):
                if in_og_table:
                    break
                continue

            if not in_og_table:
                if line.startswith("Number of species in orthogroup"):
                    in_og_table = True
                    continue

                if "\t" in line:
                    name, val = line.rsplit("\t", 1)
                    name = name.strip()
                    val = val.strip()
                    try:
                        overall_stats_info[name] = float(val)
                    except ValueError:
                        overall_stats_info[name] = val
                continue

            # ---------- In OG table ----------
            parts = line.split()
            if len(parts) < 2:
                break
            try:
                og_idx = int(parts[0])
                count  = int(parts[-1])
            except ValueError:
                break

            og_nspecies[og_idx] = count

    return overall_stats_info, og_nspecies

def get_og_nspecies_info(file_path):
    _, og_nspecies = get_stats_info(file_path)

    return og_nspecies

def get_overall_stats_info(file_path):
    overall_stats_info, _ = get_stats_info(file_path)
    return overall_stats_info

def get_expected_results(expected_results_dir):
    # files = glob.glob(expected_results_dir, recursive = True)
    # for file in files:
    #     print(file)
    # expected_results_dir = r"./tests/data/proteome_expected_results"
    for cat_dir in os.listdir(expected_results_dir):
        test_category = cat_dir.lower()
        for item in os.listdir(os.path.join(expected_results_dir, cat_dir)):
            if cat_dir.lower() in ["core", "assign"]:
                test_item = os.path.join(expected_results_dir, cat_dir, item)
                if os.path.isfile(test_item):
                    filename, extension = item.rsplit(".", 1)
                    filename = filename.lower()
                    # ------------ Reconciliaton tree ----------------
                    if "resolved_gene_tree" in filename:
                        with open(test_item) as reader:
                            for line in reader:
                                if len(line.strip()) == 0 or "#" in line:
                                    continue
                                og, targets = line.strip().split(":", 1)
                                EXPECTED_RESULTS[cat_dir.lower()]["resolved_gene_trees"][og] = set(re.split(r"[,\s;]+", targets.strip()))
                    # --------------- Orthogroups -------------------
                    elif "orthogroup" in filename:
                        with open(test_item) as reader:
                            for line in reader:
                                if len(line.strip()) == 0 or "#" in line:
                                    continue
                                # if "txt" in extension:
                                #     og, targets = line.strip().split(":", 1)
                                #     EXPECTED_RESULTS[cat_dir.lower()]["orthogroups"][og] = set(re.split(r"[,\s;]+", targets.strip()))
                                # elif "tsv" in extension:
                                if "Orthogroup" in line or "#" in line:
                                    continue
                                line = re.split(r"[,\s;]+", line.strip())
                    
                                EXPECTED_RESULTS[cat_dir.lower()]["orthogroups"][line[0]] = \
                                tuple(
                                    frozenset(re.split(r"[,\s;]+", l.strip()))
                                    for l in line[1:]
                                )
                    # --------------- Duplications -------------------
                    elif "duplication" in filename:
                        with open(test_item) as reader:
                            for line in reader:
                                if "Orthogroup" in line or len(line) < 7 or "#" in line:
                                    continue
                                line = re.split(r"[\t]+", line.strip())
                                terminal_type = line[-3]
                                if line[0] not in EXPECTED_RESULTS[cat_dir.lower()]["duplications"]:
                                    EXPECTED_RESULTS[cat_dir.lower()]["duplications"][line[0]] = \
                                    [
                                        (
                                            terminal_type, 
                                            frozenset(re.split(r"[,\s;]+", line[-2].strip())), 
                                            frozenset(re.split(r"[,\s;]+", line[-1].strip()))
                                        )
                                    ]
                                elif line[0] in EXPECTED_RESULTS[cat_dir.lower()]["duplications"]:
                                    EXPECTED_RESULTS[cat_dir.lower()]["duplications"][line[0]].append(
                                        (
                                            terminal_type, 
                                            frozenset(re.split(r"[,\s;]+", line[-2].strip())), 
                                            frozenset(re.split(r"[,\s;]+", line[-1].strip()))
                                        )
                                    )
                    # --------------- Stats overall -------------------
                    elif "statistics_overall" in filename:
                        overall_stats_info, og_nspecies = get_stats_info(test_item)
                        
                        EXPECTED_RESULTS[cat_dir.lower()]["statistics_overall"]["stats_overall"] = overall_stats_info
                        EXPECTED_RESULTS[cat_dir.lower()]["statistics_overall"]["og_nspecies"] = og_nspecies

                elif os.path.isdir(test_item):
                    # --------------- Orthologues -------------------
                    if item.lower() == "orthologues":
                        for file in os.listdir(test_item):
                            file_path = os.path.join(test_item, file)
                            base_species = ""
                            with open(file_path) as reader:
                                for line in reader:
                                    line = re.split(r"[\t]+", line.strip())
                                    if "Orthogroup" in line or len(line) < 4 or "#" in line:
                                        base_species = line[2]
                                        EXPECTED_RESULTS[cat_dir.lower()]["orthologues"][base_species] = {}
                                    else:
                                        if line[0] not in EXPECTED_RESULTS[cat_dir.lower()]["orthologues"][base_species]:
                                             
                                            EXPECTED_RESULTS[cat_dir.lower()]["orthologues"][base_species][line[0]] = \
                                                [(line[1], frozenset(re.split(r"[,\s;]+", line[2])),  frozenset(re.split(r"[,\s;]+", line[-1])))]
                                        else:
                                            EXPECTED_RESULTS[cat_dir.lower()]["orthologues"][base_species][line[0]].append(
                                                (line[1], frozenset(re.split(r"[,\s;]+", line[2])),  frozenset(re.split(r"[,\s;]+", line[-1])))
                                            )

                    elif item.lower() == "hogs":
                        # --------------- Hogs -------------------
                        EXPECTED_RESULTS[cat_dir.lower()]["hogs"] = defaultdict(lambda: defaultdict(list))
                        for file in os.listdir(test_item):
                            file_path = os.path.join(test_item, file)
                            rootname = file.rsplit(".", 1)[0]
                            base_node = None

                            for part in rootname.split("_"):
                                m = node_pattern.match(part)
                                if m:
                                    base_node = f"N{int(m.group(1))}"
                                    break

                            if base_node is None:
                                base_node = rootname
        
                            with open(file_path) as reader:
                                for line in reader:
                                    if line.lower().startswith("hog") or len(line) < 4 or "#" in line:
                                        continue 
                                    line = re.split(r"[\t]+", line.strip())
                                    ogname = line[1]
                                    # if ogname not in EXPECTED_RESULTS[cat_dir.lower()]["hogs"] and \
                                    #     base_node not in EXPECTED_RESULTS[cat_dir.lower()]["hogs"][ogname]:
                                    #     EXPECTED_RESULTS[cat_dir.lower()]["hogs"][ogname][base_node] = [
                                    #         tuple(
                                    #             frozenset(re.split(r"[,\s;]+", l.strip()))
                                    #             for l in line[3:]
                                    #             if len(l) != 0
                                    #         )
                                    #     ]
                                    # elif ogname in EXPECTED_RESULTS[cat_dir.lower()]["hogs"] and \
                                    #     base_node in EXPECTED_RESULTS[cat_dir.lower()]["hogs"][ogname]:
                                    EXPECTED_RESULTS[cat_dir.lower()]["hogs"][ogname][base_node].append(
                                        tuple(
                                            frozenset(re.split(r"[,\s;]+", l.strip()))
                                            for l in line[3:]
                                            if len(l) != 0
                                        )
                                    )
    # print(EXPECTED_RESULTS["core"]["hogs"])
    # print()
    # print(EXPECTED_RESULTS["assign"]["hogs"])
    return EXPECTED_RESULTS
                                        
                        
if __name__ == "__main__":
    get_expected_results()
