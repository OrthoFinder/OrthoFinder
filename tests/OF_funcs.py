import os
import contextlib
import re
from collections import defaultdict

import helper
from orthofinder.utils import file_handler, fasta_processor
from orthofinder.run import process_args
from orthofinder.phylo import tree
from orthofinder.species_tree import stride
from orthofinder.run.process_args import Options
from orthofinder.hog_orthogroups.convert import HogOrthogroupSet, update_ogs
from orthofinder.mcl_orthogroups.mcl import id_full_dict
from orthofinder.hog_orthogroups.update_trees import read_tree_file
from orthofinder.hog_orthogroups.finalise import index_files
from orthofinder.run.previous_run import species_name_dict, process_previous_files
from orthofinder.orthology.reconciliation import get_species_neighbours, reconciliation_and_orthologs, get_all_orthologs

from orthofinder.orthology.hog_writer import HogWriter
from orthofinder.orthology.tree_analyser import get_lines_for_ortholog_files, gene_to_species_dash, get_orthologs_from_tree
from orthofinder.orthology.rooting import check_and_root_tree
from orthofinder.utils import species_info
from orthofinder.utils import gene_ids
from orthofinder.phylo import taxa
from orthofinder.run import previous_run


class OrthoFinderTestFuncs:
    def __init__(
        self,
        projects,  # FASTA dir for core
        of_args_dict,
        baseline_options,
        assign=None,      
        species_tree=None,
        species_tree_assign=None,
        
    ):
        self.projects = projects
        self.assign = assign
        self.options = Options()
        self.options.q_start_from_fasta = True
        if "msa" in of_args_dict:
            if of_args_dict["msa"]:
                self.options.msa_program = of_args_dict["msa"]
        if "gene_tree_method" in of_args_dict:
            self.options.tree_program = of_args_dict["gene_tree_method"]
            
        if "species_tree_method" in of_args_dict:
            self.options.species_tree_program = of_args_dict["species_tree_method"]
            
        if "-t" in of_args_dict:
            self.options.n_blast = of_args_dict["sequence_search_threads"]
        if "-a" in of_args_dict:
            self.options.n_process_alg = of_args_dict["analysis_threads"]

        core_results_dir = os.path.join(self.projects, "OrthoFinder")

        core_name = baseline_options[0]
        assign_name = baseline_options[1] if len(baseline_options) > 1 else None

        self.core_results = os.environ.get(
            f"ORTHOFINDER_TEST_RESULTS_{core_name}",
            ""
        )

        if not self.core_results:
            self.core_results = helper._find_output_dir(
                core_results_dir,
                test_filename="Results_" + core_name,
                find_core=True
            )

        if not self.core_results:
            raise AssertionError(
                f"Could not find valid core OrthoFinder results for '{core_name}' "
                f"in {core_results_dir}. Run the core test first or check output naming."
            )

        self.core_working_dir = os.path.join(self.core_results, "WorkingDirectory")

        if self.assign:
            self.assign_results = ""

            if assign_name:
                self.assign_results = os.environ.get(
                    f"ORTHOFINDER_TEST_RESULTS_{assign_name}",
                    ""
                )

            if not self.assign_results and assign_name:
                self.assign_results = helper._find_output_dir(
                    core_results_dir,
                    test_filename="Results_" + assign_name,
                    find_core=False
                )

            if not self.assign_results:
                raise AssertionError(
                    f"Could not find valid assign OrthoFinder results for '{assign_name}' "
                    f"in {core_results_dir}. Run the assign test first or check output naming."
                )

            self.assign_working_dir = os.path.join(self.assign_results, "WorkingDirectory")

        if self.assign:

            self.current_results_dir = self.assign_results
            self.current_working_dir = self.assign_working_dir

            args = ["--core", self.core_results, "--assign", self.assign, "--no-print-info"]
            if species_tree_assign:
                args += ["-s", species_tree_assign]
        else:
            self.current_results_dir = self.core_results
            self.current_working_dir = self.core_working_dir
            args = ["-f", self.projects, "--no-print-info"]
            if species_tree:
                args += ["-s", species_tree]

        self.species_id_fn  = os.path.join(self.current_working_dir, "SpeciesIDs.txt")
        self.sequence_id_fn = os.path.join(self.current_working_dir, "SequenceIDs.txt")
        self.ogs_all_fn = os.path.join(self.current_working_dir, "OGsAll.tsv")

        (
            self.prog_caller,
            _opts,
            fasta_dir,
            continuation_dir,
            results_dir_non_default,
            pickle_dir_non_default,
            _user_m
        ) = process_args.process_args(args)

        # ---------- Initialise FileHandler to the **current** WD ----------
        previous_run.initialise_file_handler(
            self.options,
            fasta_dir=self.projects,
            continuation_dir=continuation_dir,
            results_dir_non_default=results_dir_non_default,
            pickle_dir_non_default=pickle_dir_non_default,
            working_dir=self.current_working_dir,
        )
        self._fh_init = self._snapshot_fh()

        # ---------- Caching / run state ----------
        self._species_cache = None          # caches speciesInfoObj
        self.active_wd_list = None          # WDs actually used this run
        # self.wd_scope_list = None
        self.added_any_new  = False         # whether assign added new species

    # ---------------- helpers ----------------
    def _snapshot_fh(self):
        return dict(
            wd_base = list(getattr(file_handler.handler, "wd_base", [])),
            wd_current= getattr(file_handler.handler, "wd_current", None),
            rd1 = getattr(file_handler.handler, "rd1", None),
            wd_trees = getattr(file_handler.handler, "wd_trees", None),
        )

    def _restore_fh(self, snap):
        file_handler.handler.wd_base = snap["wd_base"]
        file_handler.handler.wd_current = snap["wd_current"]
        file_handler.handler.rd1 = snap["rd1"]
        file_handler.handler.wd_trees = snap["wd_trees"]

    @contextlib.contextmanager
    def _use_fh_wds(self, wds):
        wds = [wd if wd.endswith(os.sep) else wd + os.sep for wd in wds]
        snap = self._snapshot_fh()
        file_handler.handler.wd_base = wds
        file_handler.handler.wd_current = wds[-1]
        file_handler.handler.rd1 = os.path.dirname(file_handler.handler.wd_current.rstrip(os.sep))
        file_handler.handler.wd_trees = file_handler.handler.wd_current
        try:
            yield
        finally:
            self._restore_fh(snap)

    def get_species_info_obj(self):
        if self._species_cache is not None:
            return self._species_cache

        if not self.assign:
            self.active_wd_list = [self.core_working_dir]
            # self.wd_scope_list = [self.core_working_dir]
            with self._use_fh_wds(self.active_wd_list):
                species_info_obj, _ = process_previous_files(
                    self.active_wd_list,
                    self.options.q_double_blast,
                    check_blast=False,
                )
            self._species_cache = species_info_obj
            return species_info_obj

        # self.wd_scope_list = [self.assign_working_dir, self.core_working_dir]
        self.active_wd_list = [self.core_working_dir, self.assign_working_dir]

        with self._use_fh_wds(self.active_wd_list[::-1]):
            species_info_obj, species_to_use_names = process_previous_files(
                self.active_wd_list[::-1],
                self.options.q_double_blast,
                check_blast=False,
            )

            # newSpeciesIDs = speciesInfoObj.speciesToUse
            # iSpecies = len(newSpeciesIDs)

            # with open(self.species_id_fn, 'a') as speciesFile:
            #     for fastaFilename in os.listdir(self.assign):
            #         fastaFilename = os.path.join(self.assign, fastaFilename)
            #         newSpeciesIDs.append(iSpecies)
            #         fastaFilename = fastaFilename.rstrip()
            #         speciesFile.write("%d: %s\n" % (iSpecies, os.path.basename(fastaFilename)))
            #         iSpecies += 1

            # # for i in range(len(speciesInfoObj.speciesToUse), len(speciesInfoObj.speciesToUse) + num_new_species):
            # #     speciesInfoObj.speciesToUse.append(i)
            # # speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1 
            # speciesInfoObj.speciesToUse = newSpeciesIDs
            # speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
    

        self._species_cache = species_info_obj
        return species_info_obj

    def get_sequence_info_obj(self):
        species_info_obj = self.get_species_info_obj()
        return species_info.get_seqs_info(self.active_wd_list[::-1], species_info_obj.species_to_use, species_info_obj.n_sp_all)

    def get_og_obj(self):
        species_info_obj = self.get_species_info_obj()
        og_set = HogOrthogroupSet(
            self.options.min_seq,
            self.active_wd_list[::-1],
            species_info_obj.species_to_use,
            species_info_obj.n_sp_all,
            self.options.q_add_species_to_ids,
            self.options.tree_program,
            id_extractor=gene_ids.FirstWordExtractor,
            species_id_fn=self.species_id_fn,
            sequence_id_fn=self.sequence_id_fn,
            ogs_all_fn=self.ogs_all_fn,        # core get_ogs_all.tsv in assign mode
            # results_dir=self.current_results_dir,
        )
        return og_set

    def get_species_name_dict(self):
        return species_name_dict(self.species_id_fn)

    def get_othologues_obj(self):
        og_set = self.get_og_obj()
        return reconciliation_and_orthologs(
            self.options.recon_method,
            og_set,
            self.options.n_blast,
            self.options.n_process_alg,
            i_species_tree=None,
            stride_dups=None,
            q_split_para_clades=False,
            save_space=False,
            print_info=True,
            exist_msa=True,
            write_hog_tree=True,
            write_to_rd=True,
        )

    def get_unique_hogs(self):
        """The orthogroups (the HOGs of N0, Orthogroups/Orthogroups.tsv) of 2 or more genes."""
        orthogroups_fn = os.path.join(self.current_results_dir, "Orthogroups", "Orthogroups.tsv")
        unique_ogs = set()
        with open(orthogroups_fn) as infile:
            next(infile)
            for line in infile:
                og, *cells = line.rstrip("\n").split("\t")
                if sum(len(c.split(", ")) for c in cells if c) >= 2:
                    unique_ogs.add(og)
        return unique_ogs
    
    def get_file_indexes(self, input_dir, extension=".fa"):
        file_index = index_files(input_dir, extension=extension)
        return file_index

    def get_recon_gene_trees(self):
        recon_tree_path = os.path.join(self.current_results_dir, "Resolved_Gene_Trees", "Resolved_Gene_Trees.txt")
        # recon_tree_wd_path = os.path.join(self.current_working_dir, "Resolved_Gene_Trees")
        tree_file_index = self.get_file_indexes(recon_tree_path, ".txt")
        unique_ogs = self.get_unique_hogs()
        gene_tree_dict = {}
        for og in unique_ogs:
            if og in tree_file_index:   # orthogroups of 2 or 3 genes have no tree
                gene_tree_dict[og] = read_tree_file(og, tree_file_index)

        return gene_tree_dict

    def get_orthologs(self):
        og_set = self.get_og_obj()
        species_tree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
        labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
        taxa.rename_tree_taxa(
            species_tree_ids_fn, labeled_tree_fn,
            og_set.species_dict(), q_support=False, q_fix_negatives=True, label='N'
        )
        species_tree_rooted_labelled = tree.Tree(species_tree_ids_fn)
        species_tree_rooted_labelled.name = "N0"
        i_node = 1
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = f"N{i_node}"
                i_node += 1
        q_no_recon  = (self.options.recon_method == "only_overlap")
        neighbours = get_species_neighbours(species_tree_rooted_labelled)

        recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
        iogs4 = og_set.get_i_ogs4()
        orthologs_alltrees = []
        for iog in iogs4:
            og_id = f"OG{iog:07d}"
            recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id+".txt")
            rooted_tree_ids, _ = check_and_root_tree(
                recon_tree_file,
                species_tree_rooted_labelled,
                gene_to_species_dash,
            )
            orthologs, _recon_tree, suspect_genes, _dups = get_orthologs_from_tree(
                iog,
                rooted_tree_ids,
                species_tree_rooted_labelled,
                gene_to_species_dash,
                neighbours,
                q_get_dups=True,
                q_no_recon=q_no_recon,
            )
            orthologs_alltrees.append((iog, orthologs))

        two_three_all, nspecies, sp_to_index = get_all_orthologs(og_set)  # ensure n==1 uses "continue"
        orthologs_alltrees.extend(two_three_all)

        ortholog_lines_tot     = ["" for _ in range(nspecies)]
        ortholog_sus_lines_tot = ["" for _ in range(nspecies)]
        get_lines_for_ortholog_files(
            orthologs_alltrees,
            og_set.species_dict(),
            og_set.species_to_use,
            og_set.sequence_dict(),
            q_contains_suspect_orthologs=False,
            ortholog_lines=ortholog_lines_tot,
            ortholog_sus_lines=ortholog_sus_lines_tot,
        )

        species_names = [og_set.species_dict()[str(sp)] for sp in og_set.species_to_use]
        orthologs_dict = defaultdict(lambda: defaultdict(list))

        for i in range(nspecies):
            species_name = species_names[i] 
            block = ortholog_lines_tot[i]   
            if not block:
                continue
            tokens = [t for t in re.split(r"[\t\n]+", block.strip()) if t]
            if not tokens:
                continue
            for k in range(0, len(tokens), 4):
                if k + 3 >= len(tokens):
                    break
                ogname, partner_name, a, b = tokens[k:k+4]
                # if "," in a:
                #     a = ", ".join([g.rsplit("_", 1)[0] for g in a.split(", ")])
                # else:
                #     a = a.rsplit("_", 1)[0] 
                # if "," in b:
                #     b = ", ".join([g.rsplit("_", 1)[0] for g in b.split(", ")])
                # else:
                #     b = b.rsplit("_", 1)[0] 
                aset = frozenset(re.split(r"[,\s;]+", a)) - {""}
                bset = frozenset(re.split(r"[,\s;]+", b)) - {""}
                
                orthologs_dict[species_name][ogname].append((partner_name, aset, bset))
        return orthologs_dict

    def get_duplications(self):
        dup_fn = os.path.join(self.current_results_dir, "Gene_Duplication_Events", "Duplications.tsv")
        
        og_set = self.get_og_obj()
        species_tree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
        labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
        taxa.rename_tree_taxa(
            species_tree_ids_fn, labeled_tree_fn,
            og_set.species_dict(), q_support=False, q_fix_negatives=True, label='N'
        )
        species_tree_rooted_labelled = tree.Tree(species_tree_ids_fn)
        species_tree_rooted_labelled.name = "N0"
        i_node = 1
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = f"N{i_node}"
                i_node += 1
        q_no_recon  = (self.options.recon_method == "only_overlap")
        neighbours = get_species_neighbours(species_tree_rooted_labelled)

        recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
        sp_tree_fn_ids = os.path.join(self.current_working_dir, "SpeciesTree_unrooted_ids.txt") 
        n_high_parallel = self.options.n_blast 
        trees_id_dir = os.path.join(self.current_working_dir, "Trees_ids")
        roots, clusters_counter, rooted_species_tree_fn, n_support, _, _, stride_dups = \
            stride.get_root(
                sp_tree_fn_ids, 
                trees_id_dir, 
                stride.gene_to_species_dash, 
                n_high_parallel, 
                q_write_rooted_tree=True
            )
        spec_seq_dict = og_set.spec_seq_dict()
        nspecies = len(og_set.species_to_use)
        spec_seq_dict[">%d genes" % (5*nspecies)] = ">%d genes" % (5*nspecies)   # used in Duplications
        iogs4 = og_set.get_i_ogs4()
        duplication_dict = {}
        for iog in iogs4:
            og_id = f"OG{iog:07d}"
            recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id+".txt")
            rooted_tree_ids, _ = check_and_root_tree(
                recon_tree_file,
                species_tree_rooted_labelled,
                gene_to_species_dash,
            )
            orthologs, _recon_tree, suspect_genes, duplications = get_orthologs_from_tree(
                iog,
                rooted_tree_ids,
                species_tree_rooted_labelled,
                gene_to_species_dash,
                neighbours,
                q_get_dups=True,
                q_no_recon=q_no_recon,
            )
            duplication_dict[og_id] = []
            for sp_node_id, gene_node_name, frac, genes0, genes1 in duplications:
                q_terminal = not sp_node_id.startswith("N")
                if stride_dups is None:
                    is_stride = "Terminal" if q_terminal else "Non-Terminal"
                else:
                    is_stride = "Terminal" if q_terminal else "Non-Terminal: STRIDE" if frozenset(genes0 + genes1) in stride_dups else "Non-Terminal"
                # gene_list0 = ", ".join([seqIDs[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                # gene_list1 = ", ".join([seqIDs[g] for g in genes1])

                # gene_set0 = frozenset([spec_seq_dict[g].rsplit("_", 1)[0] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                # gene_set1 = frozenset([spec_seq_dict[g].rsplit("_", 1)[0] for g in genes1])

                gene_set0 = frozenset([spec_seq_dict[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                gene_set1 = frozenset([spec_seq_dict[g] for g in genes1])


                duplication_dict[og_id].append(
                        (is_stride, gene_set0, gene_set1)
                )

        return duplication_dict

    def get_orthogroups(self):
        # the HOGs of N0 of the MCL orthogroups, from which the run updated the
        # orthogroups (kept with --keep-legacy)
        old_hog_n0_file = os.path.join(self.current_working_dir, "Legacy", "HOGs", "N0_mcl.ids.tsv")
        og_set = self.get_og_obj()
        seqs_info = self.get_sequence_info_obj()
        species_info_obj = self.get_species_info_obj()
        species_names_dict = self.get_species_name_dict()
        all_seq_ids = og_set.all_used_sequence_ids()
        new_ogs, name_dictionary =  update_ogs(old_hog_n0_file)
        all_assigned = set([g for og in new_ogs for g in og])
        unassigned = set(all_seq_ids).difference(all_assigned)
        single_ogs_list = [{g,} for g in unassigned]
        new_ogs.extend(single_ogs_list)
        ids_filenames = [self.sequence_id_fn]
        try:
            id_to_name_dict = id_full_dict(ids_filenames, func=gene_ids.FirstWordExtractor)
        except:
            id_to_name_dict = id_full_dict(ids_filenames, func=gene_ids.FullAccession)

        n_species = len(species_names_dict)
        
        ogs_names = [[id_to_name_dict[seq] for seq in og] for og in new_ogs]
        ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in new_ogs]

        orthogroups_dict = {}
        for i_og, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
            og_dict = defaultdict(list)
            og_id = "OG%07d" % i_og
            
            for (i_species, i_sequence), name in zip(og, og_names):
                og_dict[species_info_obj.species_to_use.index(i_species)].append(name)
            
            # orthogroups_list = [
            #     (species_id, frozenset([g.rsplit("_", 1)[0] for g in genes]))
            #     for species_id, genes in ogDict.items()
            # ]
            orthogroups_list = [
                (species_id, frozenset(genes))
                for species_id, genes in og_dict.items()
            ]

            sorted_orthogroups_list = sorted(orthogroups_list)
            _, orthogroups_set = zip(*sorted_orthogroups_list)
            orthogroups_dict[og_id] = tuple(orthogroups_set)
        return orthogroups_dict

        
    def get_hogs(self):

        og_set = self.get_og_obj()
        species_tree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
        labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
        taxa.rename_tree_taxa(
            species_tree_ids_fn, labeled_tree_fn,
            og_set.species_dict(), q_support=False, q_fix_negatives=True, label='N'
        )
        species_tree_rooted_labelled = tree.Tree(species_tree_ids_fn)
        species_tree_rooted_labelled.name = "N0"
        i_node = 1
        node_names = [species_tree_rooted_labelled.name]
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = "N%d" % i_node
                node_names.append(n.name)
                i_node += 1

        q_no_recon  = (self.options.recon_method == "only_overlap")
        neighbours = get_species_neighbours(species_tree_rooted_labelled)

        recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
        
        species_dict = og_set.species_dict()
        sequence_dict = og_set.sequence_dict()
        hog_writer = HogWriter(
            species_tree_rooted_labelled, 
            node_names, 
            sequence_dict, 
            species_dict, 
            og_set.species_to_use,
            write_to_rd=True,
            write_output=False
        )
        q_split_paralogous_clades = self.options.q_split_para_clades
        
        iogs4 = og_set.get_i_ogs4()
        hogs_dict = defaultdict(lambda: defaultdict(list))
        for iog in iogs4:
            og_id = f"OG{iog:07d}"
            recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id+".txt")
            rooted_tree_ids, _ = check_and_root_tree(
                recon_tree_file,
                species_tree_rooted_labelled,
                gene_to_species_dash,
            )
            orthologs, _recon_tree, suspect_genes, _dups = get_orthologs_from_tree(
                iog,
                rooted_tree_ids,
                species_tree_rooted_labelled,
                gene_to_species_dash,
                neighbours,
                q_get_dups=True,
                q_no_recon=q_no_recon,
            )

            recon_tree = hog_writer.mark_dups_below(_recon_tree)
            cached_hogs = []
            for n in recon_tree.traverse("preorder"):
                cached_hogs.extend(hog_writer.write_clade_v2(n, og_id, q_split_paralogous_clades))
            
            
            for h, row in cached_hogs:
                if h not in ["N0.ids", "N0"]:
                    hogs = tuple(
                        # frozenset(re.split(r"[,\s;]+", item.rsplit("_", 1)[0].strip()))
                        frozenset(re.split(r"[,\s;]+", item.strip()))
                        for item in row[2:]
                        if len(item) != 0
                    )
                    hogs_dict[og_id][h].append(hogs)

        return hogs_dict
