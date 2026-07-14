import os
import contextlib
import re
from collections import defaultdict

import helper
from orthofinder.utils import files, util, fasta_processor
from orthofinder.run import process_args
from orthofinder.tools import tree, stride
from orthofinder.run.process_args import Options
from orthofinder.file_updates.ogs import OrthoGroupsSet, update_ogs, IDFullDict
from orthofinder.file_updates.trees import read_tree_file
from orthofinder.file_updates.file_updates import id_converter, hogs_converter, read_hog_file, hog_file_over4genes, index_files
from orthofinder.run.species_info import SpeciesNameDict, ProcessPreviousFiles
from orthofinder.gene_tree_inference.trees2ologs_of import (
    GetSpeciesNeighbours, 
    ReconciliationAndOrthologues, 
    AllOrthologues
)

from orthofinder.gene_tree_inference.hog_processor import HogWriter
from orthofinder.gene_tree_inference.tree_processor import (
    GetLinesForOlogFiles, 
    CheckAndRootTree, 
    GeneToSpecies_dash, 
    GetOrthologues_from_tree, 
)


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
        self.options.qStartFromFasta = True
        if "msa" in of_args_dict:
            if of_args_dict["msa"]:
                self.options.msa_program = of_args_dict["msa"]
        if "gene_tree_method" in of_args_dict:
            self.options.tree_program = of_args_dict["gene_tree_method"]
            
        if "species_tree_method" in of_args_dict:
            self.options.species_tree_program = of_args_dict["species_tree_method"]
            
        if "-t" in of_args_dict:
            self.options.nBlast = of_args_dict["sequence_search_threads"]
        if "-a" in of_args_dict:
            self.options.nProcessAlg = of_args_dict["analysis_threads"]

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
            fastaDir,
            continuationDir,
            resultsDir_nonDefault,
            pickleDir_nonDefault,
            _userM
        ) = process_args.ProcessArgs(args)

        # ---------- Initialise FileHandler to the **current** WD ----------
        files.InitialiseFileHandler(
            self.options,
            fastaDir=self.projects,
            continuationDir=continuationDir,
            resultsDir_nonDefault=resultsDir_nonDefault,
            pickleDir_nonDefault=pickleDir_nonDefault,
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
            wd_base = list(getattr(files.FileHandler, "wd_base", [])),
            wd_current= getattr(files.FileHandler, "wd_current", None),
            rd1 = getattr(files.FileHandler, "rd1", None),
            wd_trees = getattr(files.FileHandler, "wd_trees", None),
        )

    def _restore_fh(self, snap):
        files.FileHandler.wd_base = snap["wd_base"]
        files.FileHandler.wd_current = snap["wd_current"]
        files.FileHandler.rd1 = snap["rd1"]
        files.FileHandler.wd_trees = snap["wd_trees"]

    @contextlib.contextmanager
    def _use_fh_wds(self, wds):
        wds = [wd if wd.endswith(os.sep) else wd + os.sep for wd in wds]
        snap = self._snapshot_fh()
        files.FileHandler.wd_base = wds
        files.FileHandler.wd_current = wds[-1]
        files.FileHandler.rd1 = os.path.dirname(files.FileHandler.wd_current.rstrip(os.sep))
        files.FileHandler.wd_trees = files.FileHandler.wd_current
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
                speciesInfoObj, _ = ProcessPreviousFiles(
                    self.active_wd_list,
                    self.options.qDoubleBlast,
                    check_blast=False,
                )
            self._species_cache = speciesInfoObj
            return speciesInfoObj

        # self.wd_scope_list = [self.assign_working_dir, self.core_working_dir]
        self.active_wd_list = [self.core_working_dir, self.assign_working_dir]

        with self._use_fh_wds(self.active_wd_list[::-1]):
            speciesInfoObj, speciesToUse_names = ProcessPreviousFiles(
                self.active_wd_list[::-1],
                self.options.qDoubleBlast,
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
    

        self._species_cache = speciesInfoObj
        return speciesInfoObj

    def get_sequence_info_obj(self):
        speciesInfoObj = self.get_species_info_obj()
        return util.GetSeqsInfo(self.active_wd_list[::-1], speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)

    def get_og_obj(self):
        speciesInfoObj = self.get_species_info_obj()
        ogSet = OrthoGroupsSet(
            self.options.min_seq,
            self.active_wd_list[::-1],
            speciesInfoObj.speciesToUse,
            speciesInfoObj.nSpAll,
            self.options.qAddSpeciesToIDs,
            self.options.tree_program,
            idExtractor=util.FirstWordExtractor,
            species_id_fn=self.species_id_fn,
            sequence_id_fn=self.sequence_id_fn,
            ogs_all_fn=self.ogs_all_fn,        # core OGsAll.tsv in assign mode
            # results_dir=self.current_results_dir,
        )
        return ogSet

    def get_species_name_dict(self):
        return SpeciesNameDict(self.species_id_fn)

    def get_othologues_obj(self):
        ogSet = self.get_og_obj()
        return ReconciliationAndOrthologues(
            self.options.recon_method,
            ogSet,
            self.options.nBlast,
            self.options.nProcessAlg,
            iSpeciesTree=None,
            stride_dups=None,
            q_split_para_clades=False,
            fewer_open_files=False,
            save_space=False,
            old_version=False,
            print_info=True,
            exist_msa=True,
            write_hog_tree=True,
            fix_files=True,
            write_to_rd=True,
        )

    def get_unique_hogs(self):
        hog_n0_file = os.path.join(self.current_results_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
        hog_n0 = read_hog_file(hog_n0_file)
        hog_n0_over4genes = hog_file_over4genes(hog_n0, 2)

        unique_ogs = set(d['OG'] for d in hog_n0_over4genes)
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
            gene_tree_dict[og] = read_tree_file(og, tree_file_index)

        return gene_tree_dict

    def get_orthologues(self):
        ogSet = self.get_og_obj()
        speciesTree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
        labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
        util.RenameTreeTaxa(
            speciesTree_ids_fn, labeled_tree_fn,
            ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N'
        )
        species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
        species_tree_rooted_labelled.name = "N0"
        iNode = 1
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = f"N{iNode}"
                iNode += 1
        qNoRecon  = (self.options.recon_method == "only_overlap")
        neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)

        recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
        iogs4 = ogSet.Get_iOGs4()
        orthologues_alltrees = []
        for iog in iogs4:
            og_id = f"OG{iog:07d}"
            recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id+".txt")
            rooted_tree_ids, _ = CheckAndRootTree(
                recon_tree_file,
                species_tree_rooted_labelled,
                GeneToSpecies_dash,
            )
            ologs, _recon_tree, suspect_genes, _dups = GetOrthologues_from_tree(
                iog,
                rooted_tree_ids,
                species_tree_rooted_labelled,
                GeneToSpecies_dash,
                neighbours,
                q_get_dups=True,
                qNoRecon=qNoRecon,
            )
            orthologues_alltrees.append((iog, ologs))

        two_three_all, nspecies, sp_to_index, _, _ = AllOrthologues(ogSet)  # ensure n==1 uses "continue"
        orthologues_alltrees.extend(two_three_all)

        fewer_open_files = True
        dim2 = 1 if fewer_open_files else nspecies
        olog_lines_tot     = [["" for _ in range(dim2)] for _ in range(nspecies)]
        olog_sus_lines_tot = ["" for _ in range(nspecies)]
        GetLinesForOlogFiles(
            orthologues_alltrees,
            ogSet.SpeciesDict(),
            ogSet.speciesToUse,
            ogSet.SequenceDict(),
            qContainsSuspectOlogs=False,
            olog_lines=olog_lines_tot,
            olog_sus_lines=olog_sus_lines_tot,
            fewer_open_files=fewer_open_files,
        )

        species_names = [ogSet.SpeciesDict()[str(sp)] for sp in ogSet.speciesToUse]
        orthologues_dict = defaultdict(lambda: defaultdict(list))

        for i in range(nspecies):
            species_name = species_names[i] 
            block = olog_lines_tot[i][0]   
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
                
                orthologues_dict[species_name][ogname].append((partner_name, aset, bset))
        return orthologues_dict

    def get_duplications(self):
        dup_fn = os.path.join(self.current_results_dir, "Gene_Duplication_Events", "Duplications.tsv")
        
        ogSet = self.get_og_obj()
        speciesTree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
        labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
        util.RenameTreeTaxa(
            speciesTree_ids_fn, labeled_tree_fn,
            ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N'
        )
        species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
        species_tree_rooted_labelled.name = "N0"
        iNode = 1
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = f"N{iNode}"
                iNode += 1
        qNoRecon  = (self.options.recon_method == "only_overlap")
        neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)

        recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
        spTreeFN_ids = os.path.join(self.current_working_dir, "SpeciesTree_unrooted_ids.txt") 
        nHighParallel = self.options.nBlast 
        trees_id_dir = os.path.join(self.current_working_dir, "Trees_ids")
        roots, clusters_counter, rootedSpeciesTreeFN, nSupport, _, _, stride_dups = \
            stride.GetRoot(
                spTreeFN_ids, 
                trees_id_dir, 
                stride.GeneToSpecies_dash, 
                nHighParallel, 
                qWriteRootedTree=True
            )
        spec_seq_dict = ogSet.Spec_SeqDict()
        nspecies = len(ogSet.speciesToUse)
        spec_seq_dict[">%d genes" % (5*nspecies)] = ">%d genes" % (5*nspecies)   # used in Duplications
        iogs4 = ogSet.Get_iOGs4()
        duplication_dict = {}
        for iog in iogs4:
            og_id = f"OG{iog:07d}"
            recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id+".txt")
            rooted_tree_ids, _ = CheckAndRootTree(
                recon_tree_file,
                species_tree_rooted_labelled,
                GeneToSpecies_dash,
            )
            ologs, _recon_tree, suspect_genes, duplications = GetOrthologues_from_tree(
                iog,
                rooted_tree_ids,
                species_tree_rooted_labelled,
                GeneToSpecies_dash,
                neighbours,
                q_get_dups=True,
                qNoRecon=qNoRecon,
            )
            duplication_dict[og_id] = []
            for sp_node_id, gene_node_name, frac, genes0, genes1 in duplications:
                q_terminal = not sp_node_id.startswith("N")
                if stride_dups is None:
                    isSTRIDE = "Terminal" if q_terminal else "Non-Terminal"
                else:
                    isSTRIDE = "Terminal" if q_terminal else "Non-Terminal: STRIDE" if frozenset(genes0 + genes1) in stride_dups else "Non-Terminal"
                # gene_list0 = ", ".join([seqIDs[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                # gene_list1 = ", ".join([seqIDs[g] for g in genes1])

                # gene_set0 = frozenset([spec_seq_dict[g].rsplit("_", 1)[0] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                # gene_set1 = frozenset([spec_seq_dict[g].rsplit("_", 1)[0] for g in genes1])

                gene_set0 = frozenset([spec_seq_dict[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                gene_set1 = frozenset([spec_seq_dict[g] for g in genes1])


                duplication_dict[og_id].append(
                        (isSTRIDE, gene_set0, gene_set1)
                )

        return duplication_dict

    def get_orthogroups(self):
        old_hog_n0_file = os.path.join(self.current_working_dir, "Legacy", "HOGs", "N0.tsv")
        ogSet = self.get_og_obj()
        species_id_dict, sequence_id_dict = id_converter(ogSet.SpeciesDict(),  ogSet.SequenceDict())
        species_to_use = ogSet.speciesToUse
        sp_ids = ogSet.SpeciesDict()
        iSps = list(map(str, sorted(species_to_use)))   # list of strings
        species_names = [sp_ids[i] for i in iSps]

        seqsInfo = self.get_sequence_info_obj()
        speciesInfoObj = self.get_species_info_obj()
        speciesNamesDict = self.get_species_name_dict()
        all_seq_ids = ogSet.AllUsedSequenceIDs()
        hogs_converter(old_hog_n0_file, sequence_id_dict, species_id_dict, species_names)
        new_ogs, name_dictionary =  update_ogs(old_hog_n0_file)
        all_assigned = set([g for og in new_ogs for g in og])
        unassigned = set(all_seq_ids).difference(all_assigned)
        single_ogs_list = [{g,} for g in unassigned]
        new_ogs.extend(single_ogs_list)
        idsFilenames = [self.sequence_id_fn]
        try:
            idToNameDict = IDFullDict(idsFilenames, func=util.FirstWordExtractor)
        except:
            idToNameDict = IDFullDict(idsFilenames, func=util.FullAccession)

        nSpecies = len(speciesNamesDict)
        
        ogs_names = [[idToNameDict[seq] for seq in og] for og in new_ogs]
        ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in new_ogs]

        orthogroups_dict = {}
        for iOg, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
            ogDict = defaultdict(list)
            og_id = "OG%07d" % iOg
            
            for (iSpecies, iSequence), name in zip(og, og_names):
                ogDict[speciesInfoObj.speciesToUse.index(iSpecies)].append(name)
            
            # orthogroups_list = [
            #     (species_id, frozenset([g.rsplit("_", 1)[0] for g in genes]))
            #     for species_id, genes in ogDict.items()
            # ]
            orthogroups_list = [
                (species_id, frozenset(genes))
                for species_id, genes in ogDict.items()
            ]

            sorted_orthogroups_list = sorted(orthogroups_list)
            _, orthogroups_set = zip(*sorted_orthogroups_list)
            orthogroups_dict[og_id] = tuple(orthogroups_set)
        return orthogroups_dict

        
    def get_hogs(self):

        ogSet = self.get_og_obj()
        speciesTree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
        labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
        util.RenameTreeTaxa(
            speciesTree_ids_fn, labeled_tree_fn,
            ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N'
        )
        species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
        species_tree_rooted_labelled.name = "N0"
        iNode = 1
        node_names = [species_tree_rooted_labelled.name]
        for n in species_tree_rooted_labelled.traverse():
            if (not n.is_leaf()) and (not n.is_root()):
                n.name = "N%d" % iNode
                node_names.append(n.name)
                iNode += 1

        qNoRecon  = (self.options.recon_method == "only_overlap")
        neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)

        recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
        
        speciesDict = ogSet.SpeciesDict()
        SequenceDict = ogSet.SequenceDict()
        hog_writer = HogWriter(
            species_tree_rooted_labelled, 
            node_names, 
            SequenceDict, 
            speciesDict, 
            ogSet.speciesToUse,
            write_to_rd=True,
            write_output=False
        )
        q_split_paralogous_clades = self.options.qSplitParaClades
        
        iogs4 = ogSet.Get_iOGs4()
        hogs_dict = defaultdict(lambda: defaultdict(list))
        for iog in iogs4:
            og_id = f"OG{iog:07d}"
            recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id+".txt")
            rooted_tree_ids, _ = CheckAndRootTree(
                recon_tree_file,
                species_tree_rooted_labelled,
                GeneToSpecies_dash,
            )
            ologs, _recon_tree, suspect_genes, _dups = GetOrthologues_from_tree(
                iog,
                rooted_tree_ids,
                species_tree_rooted_labelled,
                GeneToSpecies_dash,
                neighbours,
                q_get_dups=True,
                qNoRecon=qNoRecon,
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
