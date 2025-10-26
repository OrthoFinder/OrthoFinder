import os, shutil, tempfile
import contextlib
import re
from collections import defaultdict

import helper
from orthofinder.utils import fasta_processor
from orthofinder.utils import files, util
from orthofinder.run import process_args
from orthofinder.tools import tree, stride
from orthofinder.run.process_args import Options
from orthofinder.file_updates.ogs import OrthoGroupsSet
from orthofinder.file_updates.trees import read_tree_file
from orthofinder.file_updates.file_updates import read_hog_file, hog_file_over4genes, index_files
from orthofinder.run.species_info import SpeciesNameDict, ProcessPreviousFiles
from orthofinder.comparative_genomics.orthologues import AllOrthologues, ReconciliationAndOrthologues
from orthofinder.gene_tree_inference.trees2ologs_of import GetLinesForOlogFiles, OrthologsFiles, GetSpeciesNeighbours, GeneToSpecies_dash, TreeAnalyser, CheckAndRootTree, GetOrthologues_from_tree, GetLinesForOlogFiles


def _read_species_names(species_id_fn: str):

    names = []
    if species_id_fn and os.path.exists(species_id_fn):
        with open(species_id_fn, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if ": " in line:
                    _, name = line.split(": ", 1)
                    names.append(name)
    return names


class OrthoFinderTestFuncs:
    def __init__(
        self,
        projects,  # FASTA dir for core
        msa,
        gene_tree_method,
        recon_method,
        sequence_search_threads,
        analysis_threads,
        assign=None,      
        species_tree=None,
        species_tree_assign=None,
    ):
        self.projects = projects
        self.assign = assign
        self.options = Options()
        self.options.qStartFromFasta = True
        self.options.msa_program = msa
        self.options.tree_program = gene_tree_method
        self.options.recon_method = recon_method
        self.options.nBlast = sequence_search_threads
        self.options.nProcessAlg = analysis_threads

        core_results_dir = os.path.join(self.projects, "OrthoFinder")

        self.core_results = helper._latest_output_dir(core_results_dir, fileno=-2)
        self.core_working_dir = os.path.join(self.core_results, "WorkingDirectory")

        if self.assign:
            self.assign_results = helper._latest_output_dir(core_results_dir, fileno=-1)
            self.assign_working_dir = os.path.join(self.assign_results, "WorkingDirectory")

            self.current_results_dir = self.assign_results
            self.current_working_dir = self.assign_working_dir

            args = ["--core", self.core_results, "--assign", self.assign]
            if species_tree_assign:
                args += ["-s", species_tree_assign]
        else:
            self.current_results_dir = self.core_results
            self.current_working_dir = self.core_working_dir
            args = ["-f", self.projects]
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
        files.FileHandler.wd_current = wds[0]
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
            with self._use_fh_wds(self.active_wd_list):
                speciesInfoObj, _ = ProcessPreviousFiles(
                    self.active_wd_list,
                    self.options.qDoubleBlast,
                    check_blast=False,
                )
            self._species_cache = speciesInfoObj
            return speciesInfoObj

        wd_scope_for_species = [self.assign_working_dir, self.core_working_dir]
        self.active_wd_list = [self.core_working_dir, self.assign_working_dir]  # used elsewhere (order as you prefer)

        with self._use_fh_wds(wd_scope_for_species):
            speciesInfoObj, _ = ProcessPreviousFiles(
                wd_scope_for_species,
                self.options.qDoubleBlast,
                check_blast=False,
            )

        self._species_cache = speciesInfoObj
        return speciesInfoObj

    def get_sequence_info_obj(self):
        speciesInfoObj = self.get_species_info_obj()
        return util.GetSeqsInfo(self.active_wd_list, speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)

    def get_og_obj(self):
        speciesInfoObj = self.get_species_info_obj()
        ogSet = OrthoGroupsSet(
            self.options.min_seq,
            self.active_wd_list,
            speciesInfoObj.speciesToUse,
            speciesInfoObj.nSpAll,
            self.options.qAddSpeciesToIDs,
            self.options.tree_program,
            idExtractor=util.FirstWordExtractor,
            species_id_fn=self.species_id_fn,
            sequence_id_fn=self.sequence_id_fn,
            ogs_all_fn=self.ogs_all_fn,        # core OGsAll.tsv in assign mode
            results_dir=self.current_results_dir,
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

                gene_set0 = frozenset([spec_seq_dict[g] for g in genes0])   # line can read ">1234 genes" for example, but this has been added to dict
                gene_set1 = frozenset([spec_seq_dict[g] for g in genes1])
                duplication_dict[og_id].append(
                        (isSTRIDE, gene_set0, gene_set1)
                )

        return duplication_dict


    # def get_orthologues(self):
    #     px_dir_path = os.path.join(self.current_results_dir, "Putative_Xenologs")
    #     dResultsOrthologues = os.path.join(self.current_results_dir, "Orthologues")
    #     reconTreesRenamedDir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees")
    #     ogSet = self.get_og_obj()
    #     iogs4 = ogSet.Get_iOGs4()

    #     speciesTree_ids_fn = os.path.join(self.current_working_dir, "SpeciesTree_rooted_ids.txt")
    #     labeled_tree_fn = os.path.join(self.current_results_dir, "Species_Tree", "SpeciesTree_rooted_node_labels.txt")
    #     util.RenameTreeTaxa(speciesTree_ids_fn, labeled_tree_fn, ogSet.SpeciesDict(), qSupport=False, qFixNegatives=True, label='N')
    #     species_tree_rooted_labelled = tree.Tree(speciesTree_ids_fn)
    #     species_tree_rooted_labelled.name = "N0"   
    #     iNode = 1
    #     node_names = [species_tree_rooted_labelled.name]
    #     for n in species_tree_rooted_labelled.traverse():
    #         if (not n.is_leaf()) and (not n.is_root()):
    #             n.name = "N%d" % iNode
    #             node_names.append(n.name)
    #             iNode += 1
    #     qNoRecon = ("only_overlap" == self.options.recon_method)
    #     neighbours = GetSpeciesNeighbours(species_tree_rooted_labelled)
    #     recon_gene_tree_id_dir = os.path.join(self.current_working_dir, "Resolved_Gene_Trees_ids")
    #     nspecies = len(ogSet.speciesToUse) 
    #     orthologues_dict = defaultdict(lambda: defaultdict(list))
    #     sp_to_index = {str(sp):i for i, sp in enumerate(ogSet.speciesToUse)}
    #     fewer_open_files = True
    #     save_space = self.options.save_space
    #     with OrthologsFiles(
    #                 dResultsOrthologues, 
    #                 ogSet.SpeciesDict(), 
    #                 ogSet.speciesToUse, 
    #                 nspecies, 
    #                 sp_to_index,
    #                 save_space, 
    #                 fewer_open_files, 
    #                 read_only=True,
    #                 putative_xenolog_dir=px_dir_path
    #         ) as (ologs_files_handles, putative_xenolog_file_handles):

    #         for iog in iogs4:
    #             og_id = "OG%07d.txt" % iog
    #             recon_tree_file = os.path.join(recon_gene_tree_id_dir, og_id)
    #             # if not os.path.exists(recon_tree_file):
    #             #     return None

    #             rooted_tree_ids, qHaveSupport = \
    #                 CheckAndRootTree(
    #                     recon_tree_file, 
    #                     species_tree_rooted_labelled, 
    #                     GeneToSpecies_dash
    #                 ) # this can be parallelised easily
                

    #             ologs, recon_tree, suspect_genes, dups = \
    #                 GetOrthologues_from_tree(
    #                     iog, 
    #                     rooted_tree_ids, 
    #                     species_tree_rooted_labelled,
    #                     GeneToSpecies_dash, 
    #                     neighbours, 
    #                     q_get_dups=True, 
    #                     qNoRecon=qNoRecon
    #             )
    #             dim2 = 1 if fewer_open_files else nspecies
    #             olog_lines = [["" for j in range(dim2)] for i in range(nspecies)]
    #             olog_sus_lines = ["" for i in range(nspecies)]

    #             nOrthologues_SpPair = GetLinesForOlogFiles(
    #                 [(iog, ologs)], 
    #                 ogSet.SpeciesDict(), 
    #                 ogSet.speciesToUse,
    #                 ogSet.SequenceDict(), 
    #                 len(suspect_genes) > 0, 
    #                 olog_lines,
    #                 olog_sus_lines, 
    #                 fewer_open_files=fewer_open_files
    #             )

    #             for i in range(nspecies):
    #                 if len(olog_lines[i][0]) > 0:
    #                     species_name = os.path.basename(ologs_files_handles[i][0].name).rsplit(".", 1)[0]
    #                     line = re.split(r"[\t\n]+", olog_lines[i][0].strip())
    #                     ogname = set([line[i] for i in range(0, len(line), 4)]).pop()
    #                     ologs = [line[i+1:i+4] for i in range(0, len(line), 4)]
    #                     ologs_set_list = [
    #                         [olog[0], set(re.split(r"[,\s;]+", olog[1])), set(re.split(r"[,\s;]+", olog[2]))]
    #                         for olog in ologs
    #                     ]
    #                     # print(species_name, ogname, ologs_set_list)
    #                     if species_name not in orthologues_dict:
    #                         orthologues_dict[species_name][ogname] = ologs_set_list
    #                     elif species_name in orthologues_dict:
    #                         orthologues_dict[species_name][ogname].extend(ologs_set_list)


    #     _, nspecies, sp_to_index, olog_lines_tot, olog_sus_lines_tot = AllOrthologues(ogSet)
    #     for i in range(nspecies):
    #         for j in range(nspecies):
    #             if len(olog_lines_tot[i][j]) > 0:
    #                 species_name = os.path.basename(olog_files_handles[i][j].name).rsplit(".", 1)[0]
    #                 line = re.split(r"[\t\n]+", olog_lines_tot[i][j].strip())
    #                 ogname = set([line[i] for i in range(0, len(line), 4)]).pop()
    #                 ologs = [line[i+1:i+4] for i in range(0, len(line), 4)]
    #                 ologs_set_list = [
    #                     [olog[0], set(re.split(r"[,\s;]+", olog[1])), set(re.split(r"[,\s;]+", olog[2]))]
    #                     for olog in ologs
    #                 ]
    #                 orthologues_dict[species_name][ogname].extend(ologs_set_list)

    #     return orthologues_dict



        # ta = TreeAnalyser(
        #         len(iogs4), 
        #         dResultsOrthologues, 
        #         reconTreesRenamedDir, 
        #         species_tree_rooted_labelled,
        #         ogSet.speciesToUse, 
        #         GeneToSpecies_dash, 
        #         ogSet.SequenceDict(), 
        #         ogSet.SpeciesDict(), 
        #         ogSet.Spec_SeqDict(), 
        #         neighbours, 
        #         qNoRecon, 
        #         outfile_dups, 
        #         stride_dups, 
        #         ologs_file_handles, 
        #         putative_xenolog_file_handles, 
        #         hog_writer, 
        #         q_split_paralogous_clades, 
        #         fewer_open_files=fewer_open_files,
        #         exist_msa=exist_msa,
        #         write_hog_tree=write_hog_tree,
        #         fix_files=fix_files
        #     )

