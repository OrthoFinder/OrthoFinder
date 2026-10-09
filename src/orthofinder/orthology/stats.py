from __future__ import absolute_import

import os
import datetime
import numpy as np
from collections import Counter, defaultdict
try:
    from rich import print
except ImportError:
    ...

from rich.table import Table
from rich.console import Console

from ..utils import util, file_handler, file_io
from ..phylo import tree


def write_orthologs_stats(og_set, NOrthoSp):
    """
    NOrthoSp is a util.NOrthoSp object
    """
    species_to_use = og_set.species_to_use
    species_dict = og_set.species_dict()
    d = file_handler.handler.get_ogs_stats_results_directory()

    # print(NOrthoSp.n)
    # print(NOrthoSp.n_121)
    # print(NOrthoSp.n_12m)
    # print( NOrthoSp.n_m21)
    # print(NOrthoSp.n_m2m)
    
    write_orthologs_matrix(d + "OrthologuesStats_Totals.tsv", NOrthoSp.n, species_to_use, species_dict)
    write_orthologs_matrix(d + "OrthologuesStats_one-to-one.tsv", NOrthoSp.n_121, species_to_use, species_dict)
    write_orthologs_matrix(d + "OrthologuesStats_one-to-many.tsv", NOrthoSp.n_12m, species_to_use, species_dict)
    write_orthologs_matrix(d + "OrthologuesStats_many-to-one.tsv", NOrthoSp.n_m21, species_to_use, species_dict)
    write_orthologs_matrix(d + "OrthologuesStats_many-to-many.tsv", NOrthoSp.n_m2m, species_to_use, species_dict)
    
    # Duplications
    node_count = defaultdict(int)
    node_count50 = defaultdict(int)
    og_count = defaultdict(int)
    og_count50 = defaultdict(int)
    if not os.path.exists(file_handler.handler.get_duplications_fn()): return
    with open(file_handler.handler.get_duplications_fn(), util.csv_read_mode) as infile:
        reader = file_io.unquoted_reader(infile)   # written unquoted
        next(reader)
        # for line in reader:
        #     try:
        #         og, node, _, support, _, _, _ = line
        #     except:
        #         print(line)
        #         raise
        for og, node, _, support, _, _, _ in reader:
            support = float(support)
            node_count[node] += 1
            og_count[og] += 1
            if support >= 0.5:
                node_count50[node] += 1
                og_count50[og] += 1
    with open(d + "Duplications_per_Species_Tree_Node.tsv", util.csv_write_mode) as outfile:
        writer = file_io.writer(outfile)
        writer.writerow(["Species Tree Node", "Duplications (all)", "Duplications (50% support)"])
#        max_node = max([int(s[1:]) for s in nodeCount.keys()])    # Get largest node number
        for node in node_count:
            writer.writerow([node, node_count[node], node_count50[node]])
    # Write on species tree
    in_tree_fn = file_handler.handler.get_species_tree_results_node_labels_fn()
    out_tree_fn = os.path.split(file_handler.handler.get_duplications_fn())[0] + "/SpeciesTree_Gene_Duplications_0.5_Support.txt"
    t = tree.Tree(in_tree_fn, format=1)
    for n in t.traverse():
        n.name = n.name + "_" + str(node_count50[n.name])
    with open(out_tree_fn, 'w') as outfile:
        outfile.write(t.write(format=1)[:-1] + t.name + ";")
    with open(d + "Duplications_per_Orthogroup.tsv", util.csv_write_mode) as outfile:
        writer = file_io.writer(outfile)
        writer.writerow(["Orthogroup", "Duplications (all)", "Duplications (50% support)"])
        if len(og_count) > 0:
            max_og = max([int(s[2:]) for s in og_count.keys()]) 
            pat = file_handler.handler.base_og_format 
            for i in range(max_og + 1):
                og = pat % i
                writer.writerow([og, og_count[og], og_count50[og]])

def write_orthologs_matrix(fn, matrix, species_to_use, species_dict):
    with open(fn, util.csv_write_mode) as outfile:
        writer = file_io.writer(outfile)
        writer.writerow([""] + [species_dict[str(index)] for index in species_to_use])
        for ii, i_sp in enumerate(species_to_use):
            overlap = [matrix[ii, jj] for jj, j_sp in enumerate(species_to_use)]
            writer.writerow([species_dict[str(i_sp)]] + overlap)   
    
def orthogroups_matrix(i_species, proper_ogs):
    species_index_dict = {i_sp: i_col for i_col, i_sp in enumerate(i_species)}
    n_species = len(i_species)
    n_groups = len(proper_ogs)
    # (i, j)-th entry of ogMatrix gives the number of genes from i in orthogroup j
    og_matrix = np.zeros((n_groups, n_species))
    for i_og, og in enumerate(proper_ogs):
        for species, _ in og:
            og_matrix[i_og, species_index_dict[species]] += 1
    return og_matrix


# def Stats_SpeciesOverlaps(fn, speciesNamesDict, iSpecies, speciesPresence):
#     """ Number of orthogroups in which each species-pair is present. Called by Stats"""
#     with open(fn, util.csv_write_mode) as outfile:
#         writer = file_io.writer(outfile)
#         writer.writerow([""] + [speciesNamesDict[index] for index in iSpecies])
#         for iSp in iSpecies:
#             overlap = [len([1 for og in speciesPresence if (iSp in og and jSp in og)]) for jSp in iSpecies]
#             writer.writerow([speciesNamesDict[iSp]] + overlap)


def stats_size_table(writer_sum, writer_sp, proper_ogs, all_genes_counter, i_species, species_presence):
    """ Overall and per-species histogram tables of orthogroup sizes. Called by Stats"""
    bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 21, 51, 101, 151, 201, 501, 1001, 9e99]
    writer_sum.writerow([])
    writer_sp.writerow([])
    n_genes_per_og = [len(og) for og in proper_ogs]
    n_genes_per_og_per_species = [[len([1 for g in og if g[0] == i_sp]) for og in proper_ogs] for i_sp in i_species]
    counters_genes_per_og_per_species = [Counter(sp_count) for sp_count in n_genes_per_og_per_species]
    counters_genes_per_og = Counter(n_genes_per_og)
    n_sp = len(i_species)
    n_ogs = len(proper_ogs)
    n_genes_total = sum([len(og) for og in proper_ogs])
    percent_format = "%0.1f"
    writer_sum.writerow(
        ["Average number of genes per-species in orthogroup", "Number of orthogroups", "Percentage of orthogroups",
         "Number of genes", "Percentage of genes"])
    # Per-species tables
    table_no = [["Number of genes per-species in orthogroup"] + ["Number of orthogroups" for _ in
                                                                 i_species]]  # number of orthogroups
    table_po = [["Number of genes per-species in orthogroup"] + ["Percentage of orthogroups" for _ in
                                                                 i_species]]  # percentage of orthogroups
    table_ng = [
        ["Number of genes per-species in orthogroup"] + ["Number of genes" for _ in i_species]]  # Number of genes
    table_pg = [["Number of genes per-species in orthogroup"] + ["Percentage of genes" for _ in
                                                                 i_species]]  # percentage of genes
    for start, end in zip(bins, bins[1:]):
        bin_name = "<1" if start == 0 else (
                    "'%d" % start) if start + 1 == end else "'%d+" % start if end == 9e99 else "%d-%d" % (
        start, end - 1)
        n_orthogroups = sum([count for size, count in counters_genes_per_og.items() if start * n_sp <= size < end * n_sp])
        n_genes = sum([size * count for size, count in counters_genes_per_og.items() if start * n_sp <= size < end * n_sp])
        row_sum = [bin_name, n_orthogroups, percent_format % (100. * n_orthogroups / n_ogs), n_genes,
                   percent_format % (100. * n_genes / n_genes_total)]
        writer_sum.writerow(row_sum)
        # Per-species stats
        if bin_name == "<1": bin_name = "'0"
        n_orthogroups_ps = [sum([number for size, number in c.items() if start <= size < end]) for c in
                           counters_genes_per_og_per_species]
        n_genes_ps = [sum([size * number for size, number in c.items() if start <= size < end]) for c in
                     counters_genes_per_og_per_species]
        table_no.append([bin_name] + n_orthogroups_ps)
        table_po.append([bin_name] + [percent_format % (100. * n / n_ogs) for n in n_orthogroups_ps])
        table_ng.append([bin_name] + n_genes_ps)
        table_pg.append(
            [bin_name] + [percent_format % (100. * n / all_genes_counter[i_sp]) for i_sp, n in zip(i_species, n_genes_ps)])
    writer_sp.writerow([])
    for r in table_no: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_po: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_ng: writer_sp.writerow(r)
    writer_sp.writerow([])
    for r in table_pg: writer_sp.writerow(r)

    # Species presence
    n = list(map(len, species_presence))
    writer_sum.writerow([])
    writer_sum.writerow(["Number of species in orthogroup", "Number of orthogroups"])
    for i in range(1, n_sp + 1):
        writer_sum.writerow([i, n.count(i)])


def stats(ogs, species_names_dict, i_species, i_results_version, fasta_writer, ids_dict, q_fast_add=False):
    """ Top-level method for calculation of stats for the orthogroups"""
    all_ogs = [[list(map(int, g.split("_"))) for g in og] for og in ogs]
    proper_ogs = [og for og in all_ogs if len(og) > 1]
    iogs_proper_ogs = [iog for iog, og in enumerate(all_ogs) if len(og) > 1]
    all_genes = [g for og in all_ogs for g in og]
    og_stats_results_dir = file_handler.handler.get_ogs_stats_results_directory()
    filename_sp = og_stats_results_dir + "Statistics_PerSpecies" + (
        "" if i_results_version == 0 else "_%d" % i_results_version) + ".tsv"
    filename_sum = og_stats_results_dir + "Statistics_Overall" + (
        "" if i_results_version == 0 else "_%d" % i_results_version) + ".tsv"
    filename_overlap = og_stats_results_dir + "Orthogroups_SpeciesOverlaps" + (
        "" if i_results_version == 0 else "_%d" % i_results_version) + ".tsv"
    filename_single_copy = file_handler.handler.get_orthogroup_results_fn_base() + "_SingleCopyOrthologues.txt"
    percent_format = "%0.1f"
    with open(filename_sp, util.csv_write_mode) as outfile_species, \
        open(filename_sum, util.csv_write_mode) as outfile_sum:
        writer_sp = file_io.writer(outfile_species)
        writer_sum = file_io.writer(outfile_sum)
        # header
        writer_sp.writerow([""] + [species_names_dict[index] for index in i_species])

        # Number of genes
        all_genes_counter = Counter([g[0] for g in all_genes])
        n_genes = sum(all_genes_counter.values())
        writer_sum.writerow(["Number of species", len(i_species)])
        writer_sp.writerow(["Number of genes"] + [all_genes_counter[i_sp] for i_sp in i_species])
        writer_sum.writerow(["Number of genes", n_genes])

        # Number of assigned/unassigned genes
        assigned_genes_counter = Counter([g[0] for og in proper_ogs for g in og])
        n_assigned = sum(assigned_genes_counter.values())
        writer_sp.writerow(["Number of genes in orthogroups"] + [assigned_genes_counter[i_sp] for i_sp in i_species])
        writer_sum.writerow(["Number of genes in orthogroups"] + [n_assigned])
        writer_sp.writerow(
            ["Number of unassigned genes"] + [all_genes_counter[i_sp] - assigned_genes_counter[i_sp] for i_sp in i_species])
        writer_sum.writerow(["Number of unassigned genes"] + [n_genes - n_assigned])
        # Percentages
        p_assigned = 100. * n_assigned / n_genes
        writer_sp.writerow(["Percentage of genes in orthogroups"] + [
            percent_format % (100. * assigned_genes_counter[i_sp] / all_genes_counter[i_sp]) for i_sp in i_species])
        writer_sum.writerow(["Percentage of genes in orthogroups", percent_format % p_assigned])
        writer_sp.writerow(["Percentage of unassigned genes"] + [
            percent_format % (100 * (1. - (float(assigned_genes_counter[i_sp]) / all_genes_counter[i_sp]))) for i_sp in
            i_species])
        writer_sum.writerow(
            ["Percentage of unassigned genes", percent_format % (100 * (1. - (float(n_assigned) / n_genes)))])

        # Number of Orthogroups
        species_presence = [set([g[0] for g in og]) for og in proper_ogs]
        n_ogs = len(proper_ogs)
        # one pass over the orthogroups rather than one per species
        n_ogs_with_species = Counter(i_sp for og_sp in species_presence for i_sp in og_sp)
        writer_sum.writerow(["Number of orthogroups", n_ogs])
        writer_sp.writerow(
            ["Number of orthogroups containing species"] + [n_ogs_with_species[i_sp] for i_sp in i_species])
        writer_sp.writerow(["Percentage of orthogroups containing species"] + [percent_format % (
            (100. * n_ogs_with_species[i_sp] / len(proper_ogs)) if len(proper_ogs) > 0 else 0.)
                                                                               for i_sp in i_species])

        # Species specific orthogroups - orthogroups-based
        species_specific_ogs_counter = Counter([next(iter(og_sp)) for og_sp in species_presence if len(og_sp) == 1])
        writer_sp.writerow(
            ["Number of species-specific orthogroups"] + [species_specific_ogs_counter[i_sp] for i_sp in i_species])
        writer_sum.writerow(["Number of species-specific orthogroups", sum(species_specific_ogs_counter.values())])

        # Species specific orthogroups - gene-based
        i_species_specific_ogs = [i for i, og_sp in enumerate(species_presence) if len(og_sp) == 1]
        n_genes_in_species_specific_ogs = Counter()
        for iog in i_species_specific_ogs:
            n_genes_in_species_specific_ogs[proper_ogs[iog][0][0]] += len(proper_ogs[iog])
        i_sp_specific_ogs_gene_counts = [n_genes_in_species_specific_ogs[i_sp] for i_sp in i_species]
        writer_sp.writerow(["Number of genes in species-specific orthogroups"] + i_sp_specific_ogs_gene_counts)
        writer_sum.writerow(["Number of genes in species-specific orthogroups", sum(i_sp_specific_ogs_gene_counts)])
        writer_sp.writerow(["Percentage of genes in species-specific orthogroups"] + [
            percent_format % (100. * n_ss / all_genes_counter[i_sp]) for n_ss, i_sp in
            zip(i_sp_specific_ogs_gene_counts, i_species)])
        writer_sum.writerow(["Percentage of genes in species-specific orthogroups",
                             percent_format % (100. * sum(i_sp_specific_ogs_gene_counts) / n_genes)])

        # 'averages'
        l = list(sorted(list(map(len, proper_ogs))))
        writer_sum.writerow(["Mean orthogroup size", "%0.1f" % np.mean(l)])
        writer_sum.writerow(["Median orthogroup size", np.median(l)])
        L = np.cumsum(l)
        j, _ = next((i, x) for i, x in enumerate(L) if x > n_assigned / 2)
        writer_sum.writerow(["G50 (assigned genes)", l[j]])
        l2 = list(reversed(list(map(len, ogs))))
        L2 = np.cumsum(l2)
        j2, _ = next((i, x) for i, x in enumerate(L2) if x > n_genes / 2)
        G50 = l2[j2]
        writer_sum.writerow(["G50 (all genes)", G50])
        writer_sum.writerow(["O50 (assigned genes)", len(l) - j])
        O50 = len(l2) - j2
        writer_sum.writerow(["O50 (all genes)", O50])

        # Single-copy orthogroups
        og_matrix = orthogroups_matrix(i_species, proper_ogs)  # use iogs_properOGs for indexing
        n_species = len(i_species)
        n_present = (og_matrix > np.zeros((1, n_species))).sum(1)
        n_complete_ogs = list(n_present).count(n_species)
        single_copy_ogs = (og_matrix == np.ones((1, n_species))).all(1).nonzero()[0]
        n_single_copy = len(single_copy_ogs)
        writer_sum.writerow(["Number of orthogroups with all species present", n_complete_ogs])
        writer_sum.writerow(["Number of single-copy orthogroups", n_single_copy])
        with open(filename_single_copy, 'w') as outfile_singlecopy:
            # outfile_singlecopy.write("\n".join(["N0.HOG%07d" % iogs_properOGs[i_] for i_ in singleCopyOGs]))
            outfile_singlecopy.write("\n".join(["OG%07d" % iogs_proper_ogs[i_] for i_ in single_copy_ogs]))
        # Link single-copy orthologs
        # g_fmt = files.FileHandler.get_results_seqs_dir_single_copy() +  "N0.HOG%07d.fa"
        g_fmt = file_handler.handler.get_results_seqs_dir_single_copy() +  "OG%07d.fa"

        for i in single_copy_ogs:
            out_fn = g_fmt % iogs_proper_ogs[i]
            og = ["%d_%d" % tuple(g) for g in proper_ogs[i]]
            fasta_writer.write_seqs_to_fasta_with_new_accessions(og, out_fn, ids_dict)

        # Results filenames
        writer_sum.writerow(["Date", str(datetime.datetime.now()).split()[0]])
        writer_sum.writerow(
            ["Orthogroups file", "Orthogroups" + ("" if i_results_version == 0 else "_%d" % i_results_version) + ".tsv"])
        writer_sum.writerow(["Unassigned genes file", "Orthogroups" + (
            "" if i_results_version == 0 else "_%d" % i_results_version) + "_UnassignedGenes.tsv"])
        writer_sum.writerow(["Per-species statistics", os.path.split(filename_sp)[1]])
        writer_sum.writerow(["Overall statistics", os.path.split(filename_sum)[1]])
        writer_sum.writerow(["Orthogroups shared between species", os.path.split(filename_overlap)[1]])

        # Sizes
        stats_size_table(writer_sum, writer_sp, proper_ogs, all_genes_counter, i_species, species_presence)
        # Stats_SpeciesOverlaps(filename_overlap, speciesNamesDict, iSpecies, speciesPresence)
    # 
    util.print_time("Done writing files")
    print("\n[bold]Summary: [/bold]")
    summary_text1 = f"[dark_goldenrod]OrthoFinder[/dark_goldenrod] assigned [deep_sky_blue2]{n_assigned}[/deep_sky_blue2] genes ([deep_sky_blue2]{p_assigned:0.1f}%[/deep_sky_blue2] of total) to [deep_sky_blue2]{n_ogs}[/deep_sky_blue2] orthogroups. "
    summary_text2 = f"[deep_sky_blue2]50%[/deep_sky_blue2] of genes were in orthogroups with [deep_sky_blue2]{G50}[/deep_sky_blue2] or more genes (G50 was [deep_sky_blue2]{G50}[/deep_sky_blue2]) and were contained in the largest [deep_sky_blue2]{O50}[/deep_sky_blue2] orthogroups (O50 was [deep_sky_blue2]{O50}[/deep_sky_blue2]). "
    summary_text3 = f"There were [deep_sky_blue2]{n_complete_ogs}[/deep_sky_blue2] orthogroups with all species present and [deep_sky_blue2]{n_single_copy}[/deep_sky_blue2] of these consisted entirely of single-copy genes."

    console = Console()
    table = Table(show_header=False, box=None, expand=False)
    # table_options.add_column("Option", justify="left", width=width, no_wrap=False, overflow="fold")
    table.add_column("Description", justify="left", overflow="fold")
    table.add_row(summary_text1)
    table.add_row(summary_text2)
    table.add_row(summary_text3)
    # console.print(table)

#     summaryText = """  [dark_goldenrod]OrthoFinder[/dark_goldenrod] assigned %d genes (%0.1f%% of total) to %d orthogroups. 
#   50%% of genes were in orthogroups with %d or more genes (G50 was %d) and were contained in the largest %d orthogroups (O50 was %d). 
#   There were %d orthogroups with all species present and %d of these consisted entirely of single-copy genes.""" % (nAssigned, pAssigned, nOgs, G50, G50, O50, O50, nCompleteOGs, nSingleCopy)
    
    if q_fast_add:
        # print("  The majority of genes have been assigned to existing orthogroups, however, the remaining clade-specific genes not seen in the core species were also analysed with the following results:")
        table.add_row(
            "  The majority of genes have been assigned to existing orthogroups, however, the remaining clade-specific genes not seen in the core species were also analysed with the following results:"
        )
    console.print(table)
    # print(summaryText)


def add_unassigned_genes(ogs, all_seq_ids):
    """Extend OGs with unassigned genes as singletons"""
    all_assigned = set([g for og in ogs for g in og])
    # in sequence-ID order, not set order (which differs between runs)
    unassigned = [g for g in dict.fromkeys(all_seq_ids) if g not in all_assigned]
    ogs.extend([{g,} for g in unassigned])
    return ogs

