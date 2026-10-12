"""The end of a run: its results files and the message that it finished."""

import os  # Y
import os.path  # Y

from rich.markup import escape

from ..orthology import ortholog_writers
from ..utils import file_handler, orthoxml, paths, printing
from . import views

try:
    from rich import print
except ImportError:
    ...


def process_output_files(ctx, species_info_obj, species_names_dict):
    """The last step of every run: the results files, -x, --tree-view and the run end."""
    options = ctx.options
    with ctx.step("Process output files"):
        if not options.save_space:
            # split up the orthologs into one file per species-pair
            ortholog_writers.split_ortholog_files(file_handler.handler.get_orthologs_directory())

        ### ------------- Compress the Gene_Trees --------------
        gene_tree_dir = file_handler.handler.get_ogs_tree_dir(q_results=True)
        # usr_gene_tree_fn = files.FileHandler.GetUserTreeFN()
        # paths.compress_files(gene_tree_dir, usr_gene_tree_fn)

        ### ------------- Comprees the Resolved_Gene_Trees ------------
        resolved_gene_tree_dir = file_handler.handler.get_ogs_recon_tree_dir(q_results=True)
        usr_resolved_gene_tree_fn = file_handler.handler.get_user_resolved_tree_fn()
        paths.compress_files(resolved_gene_tree_dir, usr_resolved_gene_tree_fn)

        ### ---------- Clean up WorkingDirectory ---------------
        if options.rm_gene_trees:
            paths.cleanup_path(gene_tree_dir)

        if options.rm_resolved_gene_trees:
            paths.cleanup_path(resolved_gene_tree_dir)

        d_results = (
            os.path.normpath(file_handler.handler.get_results_directory1()) + os.path.sep
        )

        if options.fewer_open_files and options.save_space:
            for i in range(len(species_info_obj.species_to_use)):
                sp0 = species_info_obj.species_to_use[i]
                sp0_name = species_names_dict[sp0]
                sp_path = os.path.join(
                    file_handler.handler.get_orthologs_directory(), f"{sp0_name}.tsv"
                )
                # the uncompressed file left from writing the .tsv.gz; only a
                # regular file is removed, never e.g. a folder a user extracted to
                if os.path.isfile(sp_path):
                    os.remove(sp_path)

        if options.orthoxml:
            # --orthoxml (-x): the HOGs as OrthoXML, from the final results files (the
            # directory is not made if a run removed it)
            orthogroups_orthoxml = file_handler.handler.get_orthogroup_results_fn_base() + ".orthoxml"
            if os.path.exists(orthogroups_orthoxml):
                try:
                    orthoxml.write_hog_files(
                        orthogroups_orthoxml, os.path.join(d_results, "Phylogenetic_Hierarchical_Orthogroups"))
                except (OSError, KeyError, ValueError) as e:
                    print("WARNING: The HOGs could not be written in OrthoXML format: " + escape(repr(e)))

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
