import itertools
from collections import defaultdict


POLYTOMY_CLEAN = 0
POLYTOMY_MINOR_COLLAPSIBLE = 1
POLYTOMY_MODERATE_COLLAPSIBLE = 2
POLYTOMY_HEAVY_COLLAPSIBLE = 3
POLYTOMY_NON_COLLAPSIBLE = 4


class PolytomyHandler(object):
    """
    Handles orthologue calling from non-binary gene-tree nodes.

    Policy:

        Level 0: clean unresolved speciation
            No species is repeated across child branches.
            Safe to assign orthology normally.

            Example:
                (A1, B1, C1, D1)

        Level 1: minor collapsible terminal-copy ambiguity
            Exactly one species is repeated once, and the repeated copies are
            species-pure child branches.

            This is the simplest case of:
                "is this a terminal duplication in A, or an unresolved deeper event?"

            We collapse repeated species-specific children and assign orthology.

            Example:
                (A1, A2, B1)

            Collapsed groups:
                {A1, A2}, {B1}

        Level 2: moderate collapsible terminal-copy ambiguity
            Repeated species are still confined to species-pure child branches,
            so the polytomy is structurally collapsible.

            However, more than one species may be repeated.

            Acceptance is controlled by scalable thresholds:
                max_children
                max_repeated_species
                max_total_extra_species_copies
                max_copies_per_species
                max_overlap_pair_ratio

            Example accepted in a 4-species run:
                (A1, A2, B1, B2, C1, C2, D1, D2)

            Collapsed groups:
                {A1, A2}, {B1, B2}, {C1, C2}, {D1, D2}

        Level 3: heavy collapsible terminal-copy ambiguity
            The repeated species are still species-pure, so the polytomy is
            technically collapsible, but it exceeds the moderate-collapsible thresholds.

            This is rejected by default because it is too large or too repetitive.

            Example:
                (A1, A2, A3, A4, B1, C1, D1)

            Here one species has too many copies.

        Level 4: non-collapsible lineage-assignment conflict
            At least one repeated species occurs inside a mixed-species child clade.

            In this case, simply collapsing species-specific children is not enough.
            Correct orthology would require resolving which lineage the duplicated
            copy belongs to.

            Reject.

            Example:
                ((A1, B1), A2, C1)
    """

    def __init__(
            self,
            GeneToSpecies,
            SpeciesAndGene,
            max_keep_level=POLYTOMY_MODERATE_COLLAPSIBLE,
            n_species_run=None,
            max_children=None,
            max_repeated_species=None,
            max_total_extra_species_copies=None,
            max_copies_per_species=None,
            max_overlap_pair_ratio=None,
        ):
        self.GeneToSpecies = GeneToSpecies
        self.SpeciesAndGene = SpeciesAndGene
        self.max_keep_level = max_keep_level

        if n_species_run is None:
            n_species_run = 4

        self.n_species_run = n_species_run

        defaults = self._default_moderate_thresholds(n_species_run)

        if max_children is not None:
            defaults["max_children"] = max_children

        if max_repeated_species is not None:
            defaults["max_repeated_species"] = max_repeated_species

        if max_total_extra_species_copies is not None:
            defaults["max_total_extra_species_copies"] = max_total_extra_species_copies

        if max_copies_per_species is not None:
            defaults["max_copies_per_species"] = max_copies_per_species

        if max_overlap_pair_ratio is not None:
            defaults["max_overlap_pair_ratio"] = max_overlap_pair_ratio

        self.max_children = defaults["max_children"]
        self.max_repeated_species = defaults["max_repeated_species"]
        self.max_total_extra_species_copies = defaults["max_total_extra_species_copies"]
        self.max_copies_per_species = defaults["max_copies_per_species"]
        self.max_overlap_pair_ratio = defaults["max_overlap_pair_ratio"]


    @staticmethod
    def _default_moderate_thresholds(n_species_run):
        """
        Scalable default thresholds for moderate collapsible polytomies.

        The defaults scale with the number of species in the full run.

        Thresholds:

            max_children
                Maximum number of direct child branches.

                Example:
                    (A1, A2, B1, B2) has 4 children.

            max_repeated_species
                Maximum number of species that may occur in more than one child.

                Example:
                    (A1, A2, B1, B2, C1)
                    repeated species = A, B
                    repeated_species_n = 2

            max_total_extra_species_copies
                Total extra copies across all repeated species.

                Example:
                    (A1, A2, B1, B2, C1)
                    A has 1 extra copy, B has 1 extra copy
                    total_extra_species_copies = 2

            max_copies_per_species
                Maximum copies allowed for any one species.

                Example:
                    (A1, A2, A3, B1, C1)
                    A appears 3 times
                    max_copies_per_species = 3

            max_overlap_pair_ratio
                Fraction of child-child pairs that share species.

                Example:
                    (A1, A2, B1, B2)
                    total pairs = 4 * 3 / 2 = 6
                    overlapping pairs = A1-A2, B1-B2
                    overlap_pair_ratio = 2 / 6 = 0.333

        Default values:

            max_children = min(2 * n_species_run, 40)
            max_repeated_species = n_species_run
            max_total_extra_species_copies = n_species_run
            max_copies_per_species = 2
            max_overlap_pair_ratio = 0.30

        Interpretation:

            By default, balanced cases like this are allowed:

                (A1, A2, B1, B2, C1, C2, D1, D2)

            because each species appears at most twice.

            But species-explosion cases like this are rejected:

                (A1, A2, A3, A4, B1, C1, D1)

            because A appears four times.
        """
        return {
            # Allow roughly two terminal copies per species,
            # but do not let very large datasets create huge accepted polytomies.
            "max_children": min(2 * n_species_run, 40),

            # Many species may be repeated once in a balanced terminal-duplication
            # ambiguity.
            "max_repeated_species": n_species_run,

            # Allow every species to have one extra terminal copy.
            "max_total_extra_species_copies": n_species_run,

            # allow A1,A2 but reject A1,A2,A3,A4 by default.
            "max_copies_per_species": 2,

            # Reject highly tangled child-overlap structures.
            "max_overlap_pair_ratio": 0.30,
        }


    def orthologs_and_suspect_from_gene_lists(
            self,
            genes0,
            genes1,
            suspect_genes,
            misplaced_genes=None,
        ):
        """
        Same return format as Orthologs_and_Suspect(), but works from
        explicit gene-name lists instead of tree nodes.

        Returns:
            d0, d1, d0_sus, d1_sus
        """
        if misplaced_genes is None:
            misplaced_genes = set()

        d = [defaultdict(list) for _ in range(2)]
        d_sus = [defaultdict(list) for _ in range(2)]

        for genes, di, d_susi in zip((genes0, genes1), d, d_sus):
            for g in genes:
                if g in misplaced_genes:
                    continue

                sp, seq = self.SpeciesAndGene(g)

                if g in suspect_genes:
                    d_susi[sp].append(seq)
                else:
                    di[sp].append(seq)

        return d[0], d[1], d_sus[0], d_sus[1]


    def classify(self, children):
        """
        Classify a polytomy according to species overlap between child branches.

        Levels:

            POLYTOMY_CLEAN:
                No species repeated across child branches.

            POLYTOMY_MINOR_COLLAPSIBLE:
                One species repeated exactly twice, and both copies are
                species-pure child branches.

            POLYTOMY_MODERATE_COLLAPSIBLE:
                Repeated species are still species-pure and pass scalable
                threshold controls.

            POLYTOMY_HEAVY_COLLAPSIBLE:
                Repeated species are species-pure, but the polytomy exceeds
                scalable threshold controls.

            POLYTOMY_NON_COLLAPSIBLE:
                At least one repeated species occurs inside a mixed-species
                child branch.
        """
        child_species = []
        species_to_children = defaultdict(set)

        for i, child in enumerate(children):
            species_set = {
                self.GeneToSpecies(g)
                for g in child.get_leaf_names()
            }

            child_species.append(species_set)

            for sp in species_set:
                species_to_children[sp].add(i)

        repeated_species = {
            sp: child_ids
            for sp, child_ids in species_to_children.items()
            if len(child_ids) > 1
        }

        repeated_species_n = len(repeated_species)

        if repeated_species_n == 0:
            return POLYTOMY_CLEAN, {
                "child_species": child_species,
                "species_to_children": species_to_children,
                "repeated_species": repeated_species,
                "repeated_species_n": 0,
                "total_extra_species_copies": 0,
                "max_copies_per_species": 1,
                "overlap_pair_n": 0,
                "overlap_pair_ratio": 0.0,
                "collapsible": True,
                "reason": "clean",
            }

        total_extra_species_copies = sum(
            len(child_ids) - 1
            for child_ids in repeated_species.values()
        )

        max_copies_per_species = max(
            len(child_ids)
            for child_ids in repeated_species.values()
        )

        total_child_pair_n = len(children) * (len(children) - 1) // 2

        overlap_pair_n = sum(
            1
            for s0, s1 in itertools.combinations(child_species, 2)
            if s0.intersection(s1)
        )

        overlap_pair_ratio = (
            float(overlap_pair_n) / float(total_child_pair_n)
            if total_child_pair_n
            else 0.0
        )

        collapsible = True

        for sp, child_ids in repeated_species.items():
            for child_id in child_ids:
                if child_species[child_id] != {sp}:
                    collapsible = False
                    break

            if not collapsible:
                break

        profile = {
            "child_species": child_species,
            "species_to_children": species_to_children,
            "repeated_species": repeated_species,
            "repeated_species_n": repeated_species_n,
            "total_extra_species_copies": total_extra_species_copies,
            "max_copies_per_species": max_copies_per_species,
            "overlap_pair_n": overlap_pair_n,
            "overlap_pair_ratio": overlap_pair_ratio,
            "collapsible": collapsible,
        }

        if not collapsible:
            profile["reason"] = "non_collapsible"
            return POLYTOMY_NON_COLLAPSIBLE, profile

        if (
                repeated_species_n == 1
                and total_extra_species_copies == 1
                and max_copies_per_species == 2
            ):
            profile["reason"] = "minor_collapsible"
            return POLYTOMY_MINOR_COLLAPSIBLE, profile

        if (
                len(children) <= self.max_children
                and repeated_species_n <= self.max_repeated_species
                and total_extra_species_copies <= self.max_total_extra_species_copies
                and max_copies_per_species <= self.max_copies_per_species
                and overlap_pair_ratio <= self.max_overlap_pair_ratio
            ):
            profile["reason"] = "moderate_collapsible"
            return POLYTOMY_MODERATE_COLLAPSIBLE, profile

        profile["reason"] = "heavy_collapsible"
        return POLYTOMY_HEAVY_COLLAPSIBLE, profile


    def build_gene_groups(self, children, level, profile):
        """
        Build collapsed gene groups for accepted polytomies.

        Clean:
            each child remains one group.

        Collapsible overlap:
            repeated species-specific children are collapsed by species.
        """
        child_species = profile["child_species"]

        if level == POLYTOMY_CLEAN:
            groups = []

            for child, species_set in zip(children, child_species):
                groups.append((
                    list(child.get_leaf_names()),
                    set(species_set)
                ))

            return groups

        if level not in (
            POLYTOMY_MINOR_COLLAPSIBLE,
            POLYTOMY_MODERATE_COLLAPSIBLE,
        ):
            return []

        repeated_species = profile["repeated_species"]

        groups = []
        used_children = set()

        # Collapse repeated species-specific children.
        for sp in sorted(repeated_species):
            genes = []

            for child_id in sorted(repeated_species[sp]):
                genes.extend(children[child_id].get_leaf_names())
                used_children.add(child_id)

            groups.append((genes, {sp}))

        # Keep all non-repeated child clades.
        for i, child in enumerate(children):
            if i in used_children:
                continue

            groups.append((
                list(child.get_leaf_names()),
                set(child_species[i])
            ))

        groups.sort(key=lambda x: (sorted(x[1]), sorted(x[0])))

        return groups

    def get_orthologues(self, children, suspect_genes):

        level, profile = self.classify(children)

        if level > self.max_keep_level:
            return [], level, profile

        groups = self.build_gene_groups(children, level, profile)

        orthologue_tuples = []

        for (genes0, species0), (genes1, species1) in itertools.combinations(groups, 2):
            if species0.intersection(species1):
                continue

            orthologue_tuples.append(
                self.orthologs_and_suspect_from_gene_lists(
                    genes0,
                    genes1,
                    suspect_genes,
                    misplaced_genes=set(),
                )
            )

        return orthologue_tuples, level, profile


