"""The end of a run: its results files and the message that it finished."""

import os  # Y
import os.path  # Y

from rich.markup import escape

from ..orthology import ortholog_writers
from ..phylo import tree_pages
from ..utils import file_handler, orthoxml, paths, printing

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
            # -x: the HOGs and orthologs as OrthoXML and web pages, from the final
            # results files (the directories are not made if a run removed them)
            orthogroups_orthoxml = file_handler.handler.get_orthogroup_results_fn_base() + ".orthoxml"
            orthologs_dir = os.path.join(d_results, "Orthologs")
            if os.path.exists(orthogroups_orthoxml):
                try:
                    orthoxml.write_hogs_and_orthologs(
                        orthogroups_orthoxml,
                        os.path.join(d_results, "Phylogenetic_Hierarchical_Orthogroups"),
                        # the orthologs and xenologs only with -pof: with many species they are very large
                        None if options.save_space else orthologs_dir,
                        None if options.save_space else os.path.join(d_results, "Putative_Xenologs"),
                    )
                except (OSError, KeyError, ValueError) as e:
                    print("WARNING: The HOGs and orthologs could not be written in OrthoXML format: "
                          + escape(repr(e)))
                if options.save_space and os.path.isdir(orthologs_dir) and os.listdir(orthologs_dir):
                    print("\nAdd -pof to also write the orthologs and putative xenologs in OrthoXML format "
                          "and as web pages (Orthologs/Orthologs.*, Putative_Xenologs/Putative_Xenologs.*).")

        if options.tree_view:
            # --tree-view: a page to look at the species tree and the resolved gene trees
            n_trees = tree_pages.write_tree_pages(d_results)
            if n_trees:
                page = os.path.join(d_results, "Resolved_Gene_Trees", "Resolved_Gene_Trees.html")
                print(f"\nGene tree page ({n_trees} trees): {escape(page)}")

        # printer.print("\nResults:\n    %s" % d_results, style="path")
        if options.q_stop_after_mcl_groups:      # stopped after the orthogroups (-og)
            printing.print_run_end(d_results, "-fg", "its orthogroups")
        else:
            printing.print_run_end(d_results)
        file_handler.handler.write_to_log("OrthoFinder run completed\n", True)
