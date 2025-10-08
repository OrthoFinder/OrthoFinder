import os 

from orthofinder.utils import fasta_processor
from orthofinder.utils import files, util, matrices, split_ortholog_files
from orthofinder.run import run_commands, run_info, process_args, helpinfo

from orthofinder.tools import dendroblast
from orthofinder.run.process_args import Options
from orthofinder.file_updates.ogs import OrthoGroupsSet
from orthofinder.file_updates.file_updates import read_hog_file
from orthofinder.run.species_info import SpeciesNameDict
from orthofinder.comparative_genomics.orthologues import ReconciliationAndOrthologues
from orthofinder.orthogroups import gathering
from orthofinder.utils.util import CreateNewWorkingDirectory

class OrthoFinderTestFuncs:

    def __init__(
            self, 
            projects,
            projects_results,
            msa, 
            gene_tree, 
            recon_method, 
            sequence_search_threads,
            analysis_threads,
        ):

        self.options = Options()
        self.options.qStartFromFasta = True
        self.options.msa_program = msa
        self.options.tree_program = gene_tree
        self.options.recon_method = recon_method
        self.options.nBlast = sequence_search_threads
        self.options.nProcessAlg = analysis_threads

        # baseDirectoryName = os.path.join(projects, "OrthoFinder")
        # if baseDirectoryName[-1] == os.sep:
        #     baseDirectoryName += "Results_"
        # else:
        #     baseDirectoryName = baseDirectoryName + os.sep + "Results_"

        # projects_results = CreateNewWorkingDirectory(
        #     baseDirectoryName,
        #     qDate=True,
        #     search_program=None,
        #     msa_program=None,
        #     tree_program=None,
        #     scorematrix=None,
        #     gapopen=None,
        #     gapextend=None,
        #     extended_filename=False,
        #     makedir=False
        # )

        self.working_dir = os.path.join(projects_results, "WorkingDirectory")
        self.species_id_fn = os.path.join(self.working_dir, "SpeciesIDs.txt")
        self.sequence_id_fn = os.path.join(self.working_dir, "SequenceIDs.txt")
        self.ogs_all_fn = os.path.join(self.working_dir, "OGsAll.tsv")
        self.projects = projects
        self.projects_results = projects_results
        args = [
            "-f",
            os.path.basename(projects),
        ]
        self.args = args.copy()

        (   
            prog_caller, 
            options,
            fastaDir,
            continuationDir,
            resultsDir_nonDefault,
            pickleDir_nonDefault,
            user_specified_M,
        ) = process_args.ProcessArgs(args)

        # self.prog_caller = prog_caller
        # helpinfo.PrintHelp(prog_caller)

        files.InitialiseFileHandler(
                self.options, 
                fastaDir=self.projects, 
                continuationDir=None, 
                resultsDir_nonDefault=None, 
                pickleDir_nonDefault=None,
                working_dir=self.working_dir,
            )
        
    
        
    def get_species_info_obj(self):

        speciesInfoObj = fasta_processor.ProcessesNewFasta(
            self.projects, 
            self.options.dna, 
            species_id_fn=self.species_id_fn,
            sequence_id_fn=self.sequence_id_fn,
            working_dir=self.working_dir
            )
        return speciesInfoObj
    
    def get_sequence_info_obj(self):
        speciesInfoObj = self.get_species_info_obj()
        sequenceInfoObj = util.GetSeqsInfo([self.working_dir], speciesInfoObj.speciesToUse, speciesInfoObj.nSpAll)
        return sequenceInfoObj

    def get_og_obj(self):
        
        speciesInfoObj = self.get_species_info_obj()

        ogSet = OrthoGroupsSet(
            self.options.min_seq,
            [self.working_dir], 
            speciesInfoObj.speciesToUse,
            speciesInfoObj.nSpAll,
            self.options.qAddSpeciesToIDs,
            self.options.tree_program,
            idExtractor=util.FirstWordExtractor,
            species_id_fn=self.species_id_fn,
            sequence_id_fn=self.sequence_id_fn,
            ogs_all_fn=self.ogs_all_fn,
            results_dir=self.projects_results,
        )

        return ogSet 
    
    def get_species_name_dict(self):
        return SpeciesNameDict(self.species_id_fn)
    
    def get_othologues_obj(self):
        ogSet = self.get_og_obj()
        nHighParallel = self.options.nBlast
        nLowParallel = self.options.nProcessAlg
        orthologgues_obj = ReconciliationAndOrthologues(
            self.options.recon_method,
            ogSet,
            nHighParallel,
            nLowParallel,
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
            write_to_rd=True
        )
    
        return orthologgues_obj

