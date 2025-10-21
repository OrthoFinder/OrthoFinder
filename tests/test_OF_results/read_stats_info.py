import os 


def get_stats_info(file_path):
    overall_stats_info = {}
    og_nspecies = {}
    in_og_table = False

    with open(file_path) as reader:
        for i, line in enumerate(reader):
            line = line.strip()
            if not line:
                continue
            if i <= 17:
                name, val = line.rsplit("\t", 1)
                overall_stats_info[name] = float(val)
                continue

            if line.startswith("Number of species in orthogroup"):
                in_og_table = True
                continue

            if in_og_table:
                parts = line.split()
                if len(parts) < 2:
                    break
                og, val = parts[0], parts[-1]
                og_nspecies[int(og)] = int(val)

    return overall_stats_info, og_nspecies

def get_og_nspecies_info(file_path):
    _, og_nspecies = get_stats_info(file_path)

    return og_nspecies

def get_overall_stats_info(file_path):
    overall_stats_info, _ = get_stats_info(file_path)
    return overall_stats_info


