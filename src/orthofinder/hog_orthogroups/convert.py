import sys
import os
from operator import itemgetter

from ..gene_tree import msa_trees
from ..mcl_orthogroups.mcl import id_full_dict, orthogroup_table_cells
from ..mcl_orthogroups.orthogroups_set import Seq
from ..utils import util, file_handler, file_io

from ..utils import species_info
from ..utils import gene_ids


def post_hogs_processing(
    all_seq_ids,
    species_info_obj,
    seqs_info,
    species_names_dict,
    options,
    q_incremental=False,
):
    """
    Write OGs & statistics to results files, write Fasta files.
    Args:
        q_incremental - These are not the final orthogroups, don't write results
    """
    new_ogs, name_dictionary = update_ogs(
        file_handler.handler.n0_mcl_ids_fn()
    )
    results_base_filename = file_handler.handler.get_orthogroup_results_fn_base()
    # printing.print_underline("Writing orthogroups to file")
    all_assigned = {g for og in new_ogs for g in og}   # (not a list of every gene first)
    # Unassigned genes in sequence-ID order: iterating a set of strings gives a
    # different order (and so different OG numbers) in every run.
    unassigned = [g for g in dict.fromkeys(all_seq_ids) if g not in all_assigned]
    single_ogs_list = [{g} for g in unassigned]
    new_ogs.extend(single_ogs_list)

    with open(file_handler.handler.ogs_all_id_fn(), "w", encoding="utf-8", newline="\n") as outfile:
        for og in new_ogs:
            outfile.write(", ".join(sorted(og, key=species_info.seq_id_key)) + "\n")

    # (Orthogroups.txt gives no species: <species>_<gene>, gene_ids.species_prefixes)
    prefixes = gene_ids.species_prefixes(species_names_dict) if options.q_add_species_to_ids else None
    ids_dict = HogOrthogroupFiles.write_orthogroup_files(
        new_ogs, [file_handler.handler.get_sequence_ids_fn()], results_base_filename, prefixes
    )

    if not q_incremental:
        HogOrthogroupFiles.create_orthogroup_table(
            new_ogs,
            ids_dict,
            species_names_dict,
            species_info_obj.species_to_use,
            results_base_filename,
            write_ids=options.keep_legacy,
        )

    # Write Orthogroup FASTA files
    og_set = HogOrthogroupSet(
        options.min_seq,
        file_handler.handler.get_working_directory1_read(),
        species_info_obj.species_to_use,
        species_info_obj.n_sp_all,
        options.q_add_species_to_ids,
        options.tree_program,
        id_extractor=gene_ids.FirstWordExtractor,
    )

    ## ------------------ Fix Orthogroup_Sequences and Sequences_ids --------------------
    tree_gen = msa_trees.TreesForOrthogroups(None, None, None)

    d_seqs = file_handler.handler.get_results_seqs_dir()
    if not os.path.exists(d_seqs):
        os.mkdir(d_seqs)

    # Update Orthogroup_Sequneces (as the alignments and trees: <species>_<gene>)
    id_dict = og_set.spec_seq_dict()
    tree_gen.write_fasta_files(
        file_handler.handler.get_species_seqs_dir(), species_info_obj.species_to_use,
        og_set.get_ogs_all(), id_dict, q_id=False, q_results=True,
        n_processes=options.n_process_alg or 1,
    )

    id_dict.update(
        og_set.species_dict()
    )  # same code will then also convert concatenated alignment for species tree

    # Update Orthogroup_Sequneces and Sequences_ids
    # treeGen.write_fasta_files(fastaWriter, ogSet.OGssAll(), idDict, True)
    # treeGen.write_fasta_files(fastaWriter, ogSet.get_ogs_all(), idDict, False) # Set to False, only update the Orthogroup_Sequneces

    # if not q_incremental:
    ## stats.Stats(ogs, speciesNamesDict, speciesInfoObj.speciesToUse, files.FileHandler.iResultsVersion)
    # util.print_time("Done orthogroups")
    # files.FileHandler.log_ogs()

    return og_set, id_dict, name_dictionary, new_ogs


def update_ogs(input_path):
    """
    The N0 HOGs as orthogroups, most genes first: (list of gene sets, {OG:
    [[new OG name, HOG, parent clade], ...]}).
    """
    name_dictionary = {}
    new_og_list = []
    # In sorted HOG order, name each HOG's orthogroup by its position
    for pos, (count, hog, og, parent, genes) in enumerate(read_hogs_to_matrix(input_path)):
        name_dictionary.setdefault(og, []).append(["OG%07d" % pos, hog, parent])
        new_og_list.append(genes)
    return new_og_list, name_dictionary


def read_hogs_to_matrix(input_path):
    """
    The HOGs of an N0.tsv as (gene count, HOG, OG, parent clade, set of
    genes), most genes first (equal counts in file order).

    Only the genes are kept, not the row of cells: with many species almost
    every cell is empty, and a list of them per HOG would take most of the
    memory (e.g. 1,000 species x 100,000 HOGs).
    """
    rows = []
    with open(input_path, encoding="utf-8") as input_file:
        next(input_file, None)      # header
        for line in input_file:
            line_split = line.strip().split("\t")
            genes = [g for cell in line_split[3:] if cell for g in cell.split(", ") if g]
            rows.append((len(genes), line_split[0], line_split[1], line_split[2], set(genes)))
    rows.sort(key=itemgetter(0), reverse=True)
    return rows


# ==============================================================================================================================


class HogOrthogroupSet(object):
    def __init__(
        self,
        min_seq,
        orthofinder_working_dir_list,
        species_to_use,
        n_sp_all,
        q_add_species_to_ids,
        tree_program="fasttree",
        id_extractor=gene_ids.FirstWordExtractor,
        species_id_fn="",
        sequence_id_fn="",
        ogs_all_fn="",
    ):

        if not species_id_fn:
            self.species_id_fn = file_handler.handler.get_species_ids_fn()
        else:
            self.species_id_fn = species_id_fn

        if not sequence_id_fn:
            self.sequence_id_fn = file_handler.handler.get_sequence_ids_fn()
        else:
            self.sequence_id_fn = sequence_id_fn

        if not ogs_all_fn:
            self.ogs_all_fn = file_handler.handler.ogs_all_id_fn()
        else:
            self.ogs_all_fn = ogs_all_fn

        self.species_names = gene_ids.species_names(self.species_id_fn)
        self._spec_seq_ids = None
        self._extractor = id_extractor
        self.seq_ids_ex = None
        self.ogs_all = None
        self.i_ogs4 = None
        self.species_to_use = species_to_use  # list of ints
        self.seqs_info = species_info.get_seqs_info(
            orthofinder_working_dir_list, self.species_to_use, n_sp_all
        )
        self.id_to_og = None
        self.q_add_species_to_ids = q_add_species_to_ids
        self.cached_seq_ids_dict = None
        self.min_seq = min_seq
        self.tree_program = tree_program

    def sequence_dict(self):
        """returns Dict[str, str]"""
        if self.cached_seq_ids_dict is not None:
            return self.cached_seq_ids_dict
        if self.seq_ids_ex == None:
            try:
                self.seq_ids_ex = self._extractor(self.sequence_id_fn)
            except RuntimeError as error:
                print(str(error))
                if str(error).startswith("ERROR"):
                    file_handler.handler.log_fail_and_exit()
                else:
                    print(
                        "Tried to use only the first part of the accession in order to list the sequences in each orthogroup"
                    )
                    print(
                        "more concisely but these were not unique. The full accession line will be used instead.\n"
                    )
                    self.seq_ids_ex = gene_ids.FullAccession(self.sequence_id_fn)
        self.cached_seq_ids_dict = self.seq_ids_ex.get_id_to_name_dict()
        return self.cached_seq_ids_dict

    def species_dict(self):
        """returns Dict[str, str]"""
        return dict(self.species_names)

    def spec_seq_dict(self):
        """returns Dict[str, str]"""
        if self._spec_seq_ids != None:
            return self._spec_seq_ids
        # (a view of the names of the species used, with their species' prefix unless -X:
        # not a second copy of every gene's name)
        prefixes = gene_ids.species_prefixes(self.species_dict()) if self.q_add_species_to_ids else None
        self._spec_seq_ids = gene_ids.SpeciesNames(self.sequence_dict(), self.species_to_use, prefixes)
        return self._spec_seq_ids

    def get_i_ogs4(self):
        if self.i_ogs4 is None:
            ogs = self.get_ogs_all()
            self.i_ogs4 = [i for i, og in enumerate(ogs) if len(og) >= self.min_seq]
        return self.i_ogs4

    def get_ogs_all(self):
        if self.ogs_all is None:
            # (all of them: see orthogroups_set.OrthogroupSet.get_ogs_all), each line's
            # genes as Seqs straight away: not every gene held twice (as text, then Seqs)
            with open(self.ogs_all_fn, encoding="utf-8") as infile:
                self.ogs_all = [[Seq(g) for g in line.strip().split(", ")] for line in infile]

            # self.ogs_all = sorted(self.ogs_all, key=len, reverse=True)
        return self.ogs_all

    def all_ogs(self):
        with open(self.ogs_all_fn, encoding="utf-8") as infile:
            ogs = [og.strip().split(", ") for og in infile]
        return [[g for g in og] for og in ogs]

    def all_used_sequence_ids(self):
        ids_dict = self.sequence_dict()
        species_to_use_strings = set(map(str, self.species_to_use))   # (a set: a test per gene)
        all_ids = [
            s for s in ids_dict.keys() if s.split("_")[0] in species_to_use_strings
        ]
        return all_ids


class HogOrthogroupFiles:
    @staticmethod
    def create_ogs(predicted_ogs, output_filename, id_dict, prefixes=None):
        """prefixes: of the gene names (gene_ids.species_prefixes), None: the names alone (-X)"""
        with open(output_filename, "w", encoding="utf-8", newline="\n") as output_file:
            for i_og, og in enumerate(predicted_ogs):
                output_file.write("OG%07d: " % i_og)
                accessions = sorted(
                    [(prefixes[seq.split("_", 1)[0]] if prefixes else "") + id_dict[seq]
                     for seq in og if id_dict.get(seq) is not None]
                )
                output_file.write(" ".join(accessions))
                output_file.write("\n")

    # @staticmethod

    @staticmethod
    def write_orthogroup_files(ogs, ids_filenames, results_base_filename, prefixes=None):
        output_fn = results_base_filename + ".txt"
        try:
            # fullDict = dict()
            # for idsFilename in idsFilenames:
            #     idExtract = gene_ids.FirstWordExtractor(idsFilename)
            #     idDict = idExtract.get_id_to_name_dict()
            #     fullDict.update(idDict)

            full_dict = id_full_dict(ids_filenames, func=gene_ids.FirstWordExtractor)
            HogOrthogroupFiles.create_ogs(ogs, output_fn, full_dict, prefixes)
        except KeyError as e:
            sys.stderr.write("ERROR: Sequence ID not found in %s\n" % ids_filenames)
            sys.stderr.write(str(e) + "\n")
            file_handler.handler.log_fail_and_exit(
                ("ERROR: Sequence ID not found in %s\n" % ids_filenames) + str(e) + "\n"
            )
        except RuntimeError as error:
            print(str(error))
            if str(error).startswith("ERROR"):
                err_text = "ERROR: %s contains a duplicate ID. " % (ids_filenames)
                file_handler.handler.log_fail_and_exit(err_text)
            else:
                print(
                    "Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n"
                )
                try:
                    # fullDict = dict()
                    # for idsFilename in idsFilenames:
                    #     idExtract = gene_ids.FullAccession(idsFilename)
                    #     idDict = idExtract.get_id_to_name_dict()
                    #     fullDict.update(idDict)
                    full_dict = id_full_dict(ids_filenames, func=gene_ids.FullAccession)
                    HogOrthogroupFiles.create_ogs(ogs, output_fn, full_dict, prefixes)
                except RuntimeError as e:
                    # (a duplicate ID: gene_ids raises a RuntimeError; any other error is its own)
                    err_text = "ERROR: %s contains a duplicate ID (%s). " % (ids_filenames, e)
                    file_handler.handler.log_fail_and_exit(err_text)
        return full_dict

    @staticmethod
    def create_orthogroup_table(
        ogs, id_to_name_dict, species_names_dict, species_to_use, results_base_filename, write_ids=True
    ):
        """write_ids: also the orthogroups as sequence IDs (a legacy file, --keep-legacy: Legacy/Fixed_Orthogroups_id.tsv)"""


        # write out
        orthogroups_id_fn = file_handler.handler.ogs_id_fn(fixed=True) if write_ids else os.devnull
        output_filename = results_base_filename + ".tsv"
        output_filename_counts = results_base_filename + ".GeneCount.tsv"
        single_gene_filename = results_base_filename + "_UnassignedGenes.tsv"
        with (
            open(output_filename, util.csv_write_mode, encoding="utf-8", newline="\n") as output_file,
            open(orthogroups_id_fn, util.csv_write_mode, encoding="utf-8", newline="\n") as ogidfile,
            open(single_gene_filename, util.csv_write_mode, encoding="utf-8", newline="\n") as single_gene_file,
            open(output_filename_counts, util.csv_write_mode, encoding="utf-8", newline="\n") as out_file_counts,
        ):
            ogid_filewriter = file_io.writer(ogidfile)
            file_writer = file_io.writer(output_file)

            file_writer_counts = file_io.writer(out_file_counts)
            single_gene_writer = file_io.writer(single_gene_file)
            for writer in [ogid_filewriter, file_writer, single_gene_writer]:
                row = ["Orthogroup"] + [
                    species_names_dict[index] for index in species_to_use
                ]
                writer.writerow(row)
            file_writer_counts.writerow(row + ["Total"])

            column = {i_species: i for i, i_species in enumerate(species_to_use)}
            n_columns = len(species_to_use)
            for i_og, og in enumerate(ogs):
                og_name = "OG%07d" % i_og
                name_cells, id_cells, counts = orthogroup_table_cells(og, id_to_name_dict, column, n_columns)
                if len(og) == 1:
                    single_gene_writer.writerow([og_name] + name_cells)
                else:
                    file_writer.writerow([og_name] + name_cells)
                    ogid_filewriter.writerow([og_name] + id_cells)
                    file_writer_counts.writerow([og_name] + counts + [sum(counts)])

    # @staticmethod
    # def SingleGeneWriter(
    #     ogs,
    #     idToNameDict,
    #     speciesNamesDict,
    #     speciesToUse,
    #     resultsBaseFilename
    # ):

    #     nSpecies = len(speciesNamesDict)

    #     ogs_names = [[idToNameDict[seq] for seq in og] for og in ogs]
    #     ogs_ints = [[list(map(int, sequence.split("_"))) for sequence in og] for og in ogs]

    #     # write out
    #     outputFilename = resultsBaseFilename + ".tsv"
    #     outputFilename_counts = resultsBaseFilename + ".GeneCount.tsv"
    #     singleGeneFilename = resultsBaseFilename + "_UnassignedGenes.tsv"
    #     with open(outputFilename, util.csv_write_mode) as outputFile, \
    #         open(singleGeneFilename, util.csv_write_mode) as singleGeneFile, \
    #             open(outputFilename_counts, util.csv_write_mode) as outFile_counts:
    #         fileWriter = file_io.writer(outputFile)
    #         fileWriter_counts = file_io.writer(outFile_counts)
    #         singleGeneWriter = file_io.writer(singleGeneFile)
    #         for writer in [fileWriter, singleGeneWriter]:
    #             row = ["Orthogroup"] + [speciesNamesDict[index] for index in speciesToUse]
    #             writer.writerow(row)
    #         fileWriter_counts.writerow(row + ['Total'])
    #         for iOg, (og, og_names) in enumerate(zip(ogs_ints, ogs_names)):
    #             ogDict = defaultdict(list)
    #             row = ["OG%07d" % iOg]
    #             thisOutputWriter = fileWriter
    #             # separate it into sequences from each species
    #             if len(og) == 1:
    #                 row.extend(['' for x in range(nSpecies)])
    #                 row[speciesToUse.index(og[0][0]) + 1] = og_names[0]
    #                 thisOutputWriter = singleGeneWriter
    #             else:
    #                 for (iSpecies, iSequence), name in zip(og, og_names):
    #                     ogDict[speciesToUse.index(iSpecies)].append(name)
    #                 for iSpecies in range(nSpecies):
    #                     row.append(", ".join(sorted(ogDict[iSpecies])))
    #                 counts = Counter([iSpecies for iSpecies, _ in og])
    #                 counts_row = [counts[iSpecies] for iSpecies in speciesToUse]
    #                 fileWriter_counts.writerow(row[:1] + counts_row + [sum(counts_row)])
    #             thisOutputWriter.writerow(row)
