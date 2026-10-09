"""Reading gene identifiers and names from sequence accessions."""


def get_id_pair_from_string(line):
    return list(map(int, line.split("_")))


class IDExtractor(object):
    """IDExtractor deals with the fact that for different datasets a user will
    want to extract a unique sequence ID from the fasta file accessions uin different
    ways."""

    def get_id_to_name_dict(self):
        raise NotImplementedError("Should not be implemented")


def clean_accession(accession):
    """
    A FASTA header as a gene name: characters that would break the output
    files are replaced by "_" (":" "," "(" ")" in trees; a tab, possible in a
    full header line, would split the name across columns in the unquoted
    HOG files).
    """
    return (
        accession.replace(":", "_")
        .replace(",", "_")
        .replace("(", "_")
        .replace(")", "_")
        .replace("\t", "_")
    )  # .replace(".", "_")


class GeneNames(object):
    """
    The gene names of one species, kept unique, also as they are written in
    the gene trees (where the characters newick cannot hold, ; [ ] = ..., are
    replaced by "_"). Two headers that differ only in such characters
    ("gene:1" and "gene,1", or "x[1]" and "x_1_") would give the same name,
    and the genes could not be told apart in the trees and orthogroup
    files: the later one gets a suffix ("gene_1_2"), the next free "_<n>" in
    the order of the input file.
    """

    _TREE_UNSAFE = None

    def __init__(self):
        self.used = set()            # the names as written in trees

    @classmethod
    def _tree_form(cls, name):
        if cls._TREE_UNSAFE is None:
            import re

            from ..phylo import newick
            cls._TREE_UNSAFE = re.compile("[" + newick._ILEGAL_NEWICK_CHARS + "]")
        return cls._TREE_UNSAFE.sub("_", name)

    def add(self, accession):
        """The unique name of the next gene with this header."""
        name = clean_accession(accession)
        tree_form = self._tree_form(name)
        if tree_form in self.used:
            n = 2
            while "%s_%d" % (tree_form, n) in self.used:
                n += 1
            name, tree_form = "%s_%d" % (name, n), "%s_%d" % (tree_form, n)
        self.used.add(tree_form)
        return name


class FullAccession(IDExtractor):
    def __init__(self, ids_filename):
        # only want the first part and nothing else (easy!)
        self.id_to_name_dict = dict()
        names_in_species = {}
        with open(ids_filename, "r") as ids_file:
            for line in ids_file:
                line = line.rstrip()
                if not line:
                    continue
                #                if line.startswith("#"): continue
                id, accession = line.split(": ", 1)
                id = id.replace("#", "")
                id = id.strip()
                if id in self.id_to_name_dict:
                    raise RuntimeError(
                        "ERROR: A duplicate id was found in the fasta files: % s" % id
                    )
                # Names are unique within a species (see GeneNames)
                i_sp = id.split("_")[0]
                if i_sp not in names_in_species:
                    names_in_species[i_sp] = GeneNames()
                self.id_to_name_dict[id] = names_in_species[i_sp].add(accession)

    def get_id_to_name_dict(self):
        return self.id_to_name_dict


class FirstWordExtractor(IDExtractor):
    def __init__(self, ids_filename):
        # only want the first part and nothing else (easy!)
        self.id_to_name_dict = dict()
        accs_in_species = []
        with open(ids_filename, "r") as ids_file:
            for line in ids_file:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                i_sp = int(id.split("_")[0])
                while len(accs_in_species) < i_sp + 1:
                    accs_in_species.append(set())
                # Replace problematic characters
                accession = (
                    accession.replace(":", "_")
                    .replace(",", "_")
                    .replace("(", "_")
                    .replace(")", "_")
                )  # .replace(".", "_")
                # Only ensure the accessions are unique within the species, there's no need to ensure global uniqueness
                if accession in accs_in_species[i_sp]:
                    raise RuntimeError(
                        "A duplicate accession was found using just first part: % s"
                        % accession
                    )
                accs_in_species[i_sp].add(accession)
                if id in self.id_to_name_dict:
                    raise RuntimeError(
                        "ERROR: A duplicate id was found in the fasta files: % s" % id
                    )
                self.id_to_name_dict[id] = accession

    def get_id_to_name_dict(self):
        return self.id_to_name_dict
