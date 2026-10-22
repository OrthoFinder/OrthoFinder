"""The end of a run: its results files and the message that it finished."""

import os  # Y
import os.path  # Y
import shutil

from rich.markup import escape

from collections import Counter

from ..mcl_orthogroups import infer_orthogroups, mcl
from ..orthology import ortholog_writers
from ..utils import file_handler, file_io, gene_ids, orthoxml, paths, printing
from . import previous_run, views

try:
    from rich import print
except ImportError:
    ...


def write_mcl_orthogroups(options):
    """
    The MCL orthogroups (before the HOGs) in the results, from the clusters in
    the working directory: Orthogroups/MCL_Orthogroups.txt, .tsv, .GeneCount.tsv,
    _UnassignedGenes.tsv, and with --orthoxml MCL_Orthogroups.orthoxml. For a
    stop option (-og, -ogt, -ost), -mog, or a run that failed after the
    clustering. With --assign: its orthogroups before the HOGs (the core run's
    with the genes assigned to them, and the new clade-specific ones). Returns
    the .tsv file, or None (no clusters, or already written).
    """
    handler = file_handler.handler
    if options is None or handler.rd1 is None or handler.clusters_filename is None:
        return None
    clusters_fn = handler.clusters_filename + "_id_pairs.txt"
    if not os.path.isfile(clusters_fn):
        return None
    base = handler.get_mcl_orthogroup_results_fn_base()
    if os.path.isfile(base + ".tsv"):
        return None
    species_ids_fn, sequence_ids_fn = handler.get_species_ids_fn(), handler.get_sequence_ids_fn()
    species_names_dict = previous_run.species_name_dict(species_ids_fn)
    species_to_use = sorted(species_names_dict)
    with file_io.open_text_any(sequence_ids_fn) as infile:
        n_seqs = Counter(int(line.split("_", 1)[0]) for line in infile if line.strip())
    ogs = mcl.get_predicted_ogs(clusters_fn)
    infer_orthogroups.add_unclustered_genes(ogs, species_to_use, n_seqs)
    try:
        ids_dict = mcl.id_full_dict([sequence_ids_fn], func=gene_ids.FirstWordExtractor)
    except RuntimeError:   # the first words are not unique: the full accessions
        ids_dict = mcl.id_full_dict([sequence_ids_fn], func=gene_ids.FullAccession)
    # (no species given: <species>_<gene> if two species have genes of the same name)
    prefixes = gene_ids.prefixes_if_shared(ids_dict, species_names_dict, species_to_use, options.q_add_species_to_ids)
    mcl.MCL.create_ogs(ogs, base + ".txt", None, ids_dict, prefixes)
    mcl.MCL.create_orthogroup_table(ogs, ids_dict, species_names_dict, species_to_use, base, write_ids=False)
    if options.orthoxml:
        species = orthoxml.species_for_analysis(species_ids_fn, species_to_use)
        gene_ids_xml = orthoxml.read_genes([base + ".tsv", base + "_UnassignedGenes.tsv"])
        orthoxml.write_flat(base + ".tsv", gene_ids_xml, base + ".orthoxml", species)
    return base + ".tsv"


def remove_empty_results_dirs():
    """
    Remove the results directories with no files in them (not WorkingDirectory),
    of a run that stopped early or failed: a step may make one before it knows
    whether it will write to it, e.g. the orthologs before -ogt or a failure.
    """
    results_dir = file_handler.handler.rd1
    if not results_dir or not os.path.isdir(results_dir):
        return
    for entry in os.scandir(results_dir):
        if entry.is_dir(follow_symlinks=False) and entry.name != "WorkingDirectory" \
                and not any(files for _, _, files in os.walk(entry.path)):
            shutil.rmtree(entry.path)


def clean_up_at_exit(options, q_failed):
    """
    At the end of a run that stopped early (a stop option) or failed: the MCL
    orthogroups, and no empty results directories, without hiding the run's own
    error if this fails.
    """
    stopped = options is not None and (options.q_stop_after_mcl_groups
                                       or options.q_stop_after_trees or options.q_stop_after_species_trees)
    try:
        if q_failed or stopped:
            write_mcl_orthogroups(options)
    except Exception as e:
        print("WARNING: The MCL orthogroups could not be written to the results: " + escape(repr(e)))
    try:
        remove_empty_results_dirs()
    except OSError:
        pass


def process_output_files(ctx, species_info_obj, species_names_dict):
    """The last step of every run: the results files, -x, --tree-view and the run end."""
    options = ctx.options
    with ctx.step("Process output files"):
        # The results directories are not made here: only those a step wrote (e.g. none
        # of the orthologs and trees after -og)
        d_orthologs = os.path.join(file_handler.handler.rd1, "Orthologs")
        if not options.save_space and os.path.isdir(d_orthologs):
            # split up the orthologs into one file per species-pair
            ortholog_writers.split_ortholog_files(d_orthologs, n_processes=options.n_process_alg or 1)

        ### ------------- Compress the Gene_Trees --------------
        gene_tree_dir = file_handler.handler.get_ogs_tree_dir(q_results=True)
        # usr_gene_tree_fn = files.FileHandler.GetUserTreeFN()
        # paths.compress_files(gene_tree_dir, usr_gene_tree_fn)

        ### ------------- Comprees the Resolved_Gene_Trees ------------
        resolved_gene_tree_dir = file_handler.handler.get_ogs_recon_tree_dir(q_results=True)
        if any(entry.is_file() for entry in os.scandir(resolved_gene_tree_dir)):
            usr_resolved_gene_tree_fn = file_handler.handler.get_user_resolved_tree_fn()
            paths.compress_files(resolved_gene_tree_dir, usr_resolved_gene_tree_fn)

        ### ---------- Clean up WorkingDirectory ---------------
        if not options.keep_gene_trees:
            paths.cleanup_path(gene_tree_dir)

        if not options.keep_resolved_gene_trees:
            paths.cleanup_path(resolved_gene_tree_dir)

        d_results = (
            os.path.normpath(file_handler.handler.get_results_directory1()) + os.path.sep
        )

        if options.mcl_orthogroups or options.q_stop_after_mcl_groups:
            # -mog, -og: the MCL orthogroups as MCL_Orthogroups.*
            write_mcl_orthogroups(options)

        if options.orthoxml:
            # --orthoxml (-x): Orthogroups.orthoxml, from the final results files
            try:
                species = orthoxml.species_for_analysis(
                    file_handler.handler.get_species_ids_fn(), species_info_obj.species_to_use)
                orthoxml.write_results(d_results, species)
            except (OSError, KeyError, ValueError) as e:
                print("WARNING: The orthogroups could not be written in OrthoXML format: " + escape(repr(e)))

        if options.q_stop_after_mcl_groups:
            # -og stops here (the other stop options and failures: clean_up_at_exit); a full
            # run's empty directories are results, e.g. no phylogenetically misplaced genes
            remove_empty_results_dirs()

        # the pages to view the results (--view: all three), linked to each other, and with
        # two or more, OrthoFinder_Results.html to open them, a tab for each
        pages, hub, every = views.write_views(d_results, options.tree_view, options.stats_view, options.groups_view)
        printing.print_pages(every if hub else pages, hub)

        # printer.print("\nResults:\n    %s" % d_results, style="path")
        if options.q_stop_after_mcl_groups:      # stopped after the orthogroups (-og)
            printing.print_run_end(d_results, "-fg", "its orthogroups")
        else:
            printing.print_run_end(d_results)
        file_handler.handler.write_to_log("OrthoFinder run completed\n", True)
