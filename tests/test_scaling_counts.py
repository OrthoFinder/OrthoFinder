"""
Counts made without orthogroups x species matrices (memory with 1000s of species):
each gives what the matrix gave.
"""
import os
import random
from types import SimpleNamespace

import numpy as np
import pytest

from orthofinder.orthology import stats
from orthofinder.species_tree import concatenated_alignment as ca


def _matrix(ogs, min_seq):
    """The orthogroup matrix the species tree selection was given before (orthogroups x species present)."""
    iogs4 = [i for i, og in enumerate(ogs) if len(og) >= min_seq]
    species = sorted({g.i_sp for i in iogs4 for g in ogs[i]})
    col = {sp: j for j, sp in enumerate(species)}
    m = np.zeros((len(iogs4), len(species)))
    for row, i in enumerate(iogs4):
        for g in ogs[i]:
            m[row, col[g.i_sp]] += 1
    return m, iogs4


def test_species_tree_occupancy_from_the_orthogroups():
    rng = random.Random(7)
    for _ in range(200):
        n_sp = rng.randint(2, 30)
        ogs = [[SimpleNamespace(i_sp=rng.randrange(n_sp)) for _ in range(rng.randint(1, 3 * n_sp))]
               for _ in range(rng.randint(1, 300))]
        m, iogs4 = _matrix(ogs, 4)
        occupancy, iogs4_b = ca.orthogroup_occupancy(ogs, 4)
        assert iogs4_b == iogs4
        if m.size == 0:
            continue
        assert ca.get_orthogroup_occupancy_info(occupancy) == ca.get_orthogroup_occupancy_info(m)
        assert (ca.determine_orthogroups_for_species_tree(occupancy, iogs4, n_ogs_min=5, n_sufficient=50)
                == ca.determine_orthogroups_for_species_tree(m, iogs4, n_ogs_min=5, n_sufficient=50))


def test_complete_and_single_copy_orthogroups_as_the_matrix_gives():
    rng = random.Random(2)
    i_species = [0, 1, 2, 3]
    for _ in range(100):
        ogs = [[[rng.choice(i_species), k] for k in range(rng.randint(2, 7))] for _ in range(50)]
        m = stats.orthogroups_matrix(i_species, ogs)
        presence = [set(g[0] for g in og) for og in ogs]
        complete = [i for i, sp in enumerate(presence) if len(sp) == len(i_species)]
        single = [i for i in complete if len(ogs[i]) == len(i_species)]
        assert complete == list(((m > 0).sum(1) == len(i_species)).nonzero()[0])
        assert single == list((m == 1).all(1).nonzero()[0])


def test_ortholog_counts_written_as_before(tmp_path):
    counts = stats.NOrthoSp(2)
    counts.add_pair(0, 1, 3, 1)
    assert counts.n.dtype == np.uint64
    fn = str(tmp_path / "t.tsv")
    stats.write_orthologs_matrix(fn, counts.n, [0, 1], {"0": "A", "1": "B"})
    assert open(fn).read() == "\tA\tB\nA\t0.0\t3.0\nB\t1.0\t0.0\n"
    # (past 2^32: uint32 would wrap round silently)
    big = stats.NOrthoSp(2)
    for _ in range(3):
        big.add_pair(0, 1, 2 ** 31, 1)
    big += counts
    assert int(big.n[0, 1]) == 3 * 2 ** 31 + 3


def test_mcl_gene_number_to_id():
    from orthofinder.mcl_orthogroups import mcl

    def linear(starts, single_id, species):
        for i, start in enumerate(starts):
            if start > single_id:
                return "%d_%d" % (species[i - 1], single_id - starts[i - 1])
        return "%d_%d" % (species[-1], single_id - starts[-1])
    rng = random.Random(5)
    for _ in range(100):
        sizes = [rng.choice([0, 1, 2, 7]) for _ in range(rng.randint(1, 30))]   # (species without genes too)
        starts = [sum(sizes[:i]) for i in range(len(sizes))]
        species = sorted(rng.sample(range(100), len(sizes)))
        for gene in range(sum(sizes)):
            assert mcl.get_id_pair(starts, gene, species) == linear(starts, gene, species)


def test_concatenated_alignment(tmp_path):
    # rows of the single-copy genes, all-gap rows for the species without one, the
    # columns with too many gaps removed
    alignments = {0: ">0_0\nAC-T\n>1_0\nA--T\n>2_0\nACGT\n", 1: ">0_1\nMK\n>1_1\nM-\n>1_2\nMK\n"}
    for i, text in alignments.items():
        (tmp_path / ("OG%d.fa" % i)).write_text(text)
    gene = lambda i_sp, i_seq: SimpleNamespace(i_sp=i_sp, to_string=lambda: "%d_%d" % (i_sp, i_seq))
    ogs = [[gene(0, 0), gene(1, 0), gene(2, 0)], [gene(0, 1), gene(1, 1), gene(1, 2)]]
    out = tmp_path / "concatenated.fa"
    ca.create_concatenated_alignment([0, 1], ogs, lambda i: str(tmp_path / ("OG%d.fa" % i)), str(out), 1.0)
    # 3 species, f_max_gap 0.5: columns with at most 1.5 gaps kept. OG1: species 1 has two
    # genes, so only species 0's row is used, its columns 2 gaps of 3: removed, as is OG0's "-,-,G"
    assert out.read_text() == ">0\nACT\n>1\nA-T\n>2\nACT\n"


def test_alignment_of_rows_of_different_lengths_not_used(tmp_path):
    fn = tmp_path / "bad.fa"
    fn.write_text(">0_0\nAAA\n>1_0\nAA\n>2_0\nAAA\n")
    # (the short row is not one of those used: the alignment is malformed all the same)
    assert ca.orthogroup_pieces(str(fn), {"0_0", "2_0"}, {"0", "1", "2"}) == ("different lengths", None)


def _fasta_records_line_by_line(fn):
    header, parts = None, []
    with open(fn, "r", encoding="utf-8") as f:
        for line in f:
            if line[0] == ">":
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].rstrip(), []
            elif header is not None:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def test_fasta_records_read_in_pieces(tmp_path, monkeypatch):
    # records and their boundaries across the pieces read; no final newline; "\r\n"; empty records
    from orthofinder.gene_tree import msa_trees
    rng = random.Random(9)
    fn = tmp_path / "t.fa"
    for _ in range(300):
        parts = [rng.choice(["", "junk\n"])]
        for r in range(rng.randint(0, 5)):
            parts.append(">%d_%d%s\n" % (rng.randint(0, 3), r, rng.choice(["", " ", "  x"])))
            parts += [rng.choice(["", "A", "MKV" * 7]) + "\n" for _ in range(rng.choice([0, 1, 3]))]
        text = "".join(parts)
        if rng.random() < 0.3:
            text = text.rstrip("\n")
        if rng.random() < 0.1:
            text = text.replace("\n", "\r\n")
        fn.write_bytes(text.encode())
        for size in (1, 3, 7, 64):
            monkeypatch.setattr(msa_trees, "FASTA_READ_CHARS", size)
            assert list(msa_trees._fasta_records(str(fn))) == list(_fasta_records_line_by_line(str(fn)))


def test_orthogroup_fasta_files_as_fasta_writer(tmp_path, capsys):
    # the files FastaWriter writes, with the sequences in memory or sorted into ranges of orthogroups
    import filecmp
    import os
    from orthofinder.gene_tree import msa_trees
    from orthofinder.mcl_orthogroups.orthogroups_set import Seq
    rng = random.Random(1)
    wd = tmp_path / "wd"
    wd.mkdir()
    genes = []
    for sp in range(5):
        with open(wd / ("Species%d.fa" % sp), "w") as f:
            for seq in range(40):
                f.write(">%d_%d\n%s\n%s\n" % (sp, seq, "MK" * rng.randint(1, 30), "A" * rng.randint(0, 5)))
                genes.append("%d_%d" % (sp, seq))
    rng.shuffle(genes)
    ogs = [[Seq(g) for g in genes[i:i + rng.randint(1, 9)]] for i in range(0, len(genes), 6)]
    ogs[2] = ogs[2] + [Seq("0_999")]                    # (a gene without a sequence)
    names = {g: "name_" + g for g in genes}
    species = list(range(5))
    out = {d: tmp_path / d for d in ("old_ids", "old_names", "mem_ids", "mem_names", "ranges_ids", "ranges_names",
                                     "parallel_ids", "parallel_names")}
    for d in out.values():
        d.mkdir()
    fw = msa_trees.FastaWriter([str(wd) + "/"], species)
    for i, og in enumerate(ogs):
        fw.write_seqs_to_fasta(og, str(out["old_ids"] / ("OG%d.fa" % i)))
        fw.write_seqs_to_fasta_with_new_accessions(og, str(out["old_names"] / ("OG%d.fa" % i)), names)
    for mode, kwargs in (("mem", {}), ("ranges", {"buffer_bytes": 200}),
                         ("parallel", {"buffer_bytes": 600, "n_processes": 3})):
        msa_trees.write_orthogroup_fasta_files([str(wd) + "/"], species, ogs, lambda i: str(out[mode + "_ids"] / ("OG%d.fa" % i)),
                                               tmp_dir=str(tmp_path), **kwargs)
        msa_trees.write_orthogroup_fasta_files([str(wd) + "/"], species, ogs, lambda i: str(out[mode + "_names"] / ("OG%d.fa" % i)),
                                               names=names, q_report_missing=False, tmp_dir=str(tmp_path), **kwargs)
    for kind in ("ids", "names"):
        files = sorted(os.listdir(out["old_" + kind]))
        for mode in ("mem", "ranges", "parallel"):
            assert sorted(os.listdir(out[mode + "_" + kind])) == files
            assert all(filecmp.cmp(out["old_" + kind] / f, out[mode + "_" + kind] / f, shallow=False) for f in files)
    assert capsys.readouterr().out.count("ERROR: 0_999 not found") == 4     # (FastaWriter, then each mode with IDs)
    assert not [f for f in os.listdir(tmp_path) if f.startswith("og_fasta_")]


def test_orthogroup_fasta_files_fail_not_hang_if_a_process_cannot_start(tmp_path):
    # The processes start with "spawn": in a script without `if __name__ == "__main__":`
    # each fails as it starts. With data beyond a pipe's buffer the main process
    # waited for them without end (the data then went through the pipe that starts
    # them): it must stop, with an error.
    import subprocess
    import sys
    wd = tmp_path / "wd"
    wd.mkdir()
    for sp in range(3):
        with open(wd / ("Species%d.fa" % sp), "w") as f:
            for seq in range(20000):
                f.write(">%d_%d\nMKV\n" % (sp, seq))
    (tmp_path / "out").mkdir()
    script = tmp_path / "unguarded.py"
    script.write_text(
        "from orthofinder.gene_tree import msa_trees\n"
        "ogs = [['%%d_%%d' %% (sp, q) for sp in range(3)] for q in range(20000)]\n"
        "msa_trees.write_orthogroup_fasta_files([%r], [0, 1, 2], ogs, lambda i: %r %% i, buffer_bytes=1000,\n"
        "                                       tmp_dir=%r, n_processes=2)\n"
        % (str(wd) + "/", str(tmp_path / "out" / "OG%d.fa"), str(tmp_path)))
    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=120)
    assert result.returncode != 0 and "BrokenProcessPool" in result.stderr
    assert not [f for f in os.listdir(tmp_path) if f.startswith("og_fasta_")]


@pytest.mark.skipif(not os.path.exists("/proc/self/smaps_rollup"), reason="needs Linux /proc")
@pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")
def test_forked_process_does_not_copy_inherited_objects_in_a_collection():
    # a garbage collection in a forked process touches every object it inherited (its copy
    # of their memory made), unless they were frozen when it started (frozen_for_fork)
    import gc
    import time
    from orthofinder.utils import parallel_task_manager

    def private_kb(pid):
        with open("/proc/%d/smaps_rollup" % pid) as f:
            return sum(int(line.split()[1]) for line in f if line.startswith(("Private_Dirty", "Private_Clean")))

    def run(r, w):
        pid = os.fork()
        if pid == 0:
            gc.collect()
            os.write(w, b"x")
            time.sleep(30)
            os._exit(0)
        gc.collect()           # (the parent collects too, while the process runs)
        os.read(r, 1)
        kb = private_kb(pid)
        os.kill(pid, 9)
        os.waitpid(pid, 0)
        return kb

    def child_copies(frozen):
        r, w = os.pipe()
        if frozen:
            with parallel_task_manager.frozen_for_fork():
                return run(r, w)
        return run(r, w)
    frozen_before = gc.get_freeze_count()
    objects = [[i, str(i)] for i in range(1_000_000)]    # (~100 MB of lists: GC-tracked)
    assert child_copies(frozen=True) < child_copies(frozen=False) / 4
    assert gc.get_freeze_count() == frozen_before      # (unfrozen again in the parent)
    del objects


def test_orthogroup_fasta_temporary_files_do_not_grow_with_species_and_ranges(tmp_path, monkeypatch):
    # the sequences sorted by range into one pack a process: not a file a range and species' file
    from orthofinder.gene_tree import msa_trees
    wd = tmp_path / "wd"
    wd.mkdir()
    for sp in range(30):
        with open(wd / ("Species%d.fa" % sp), "w") as f:
            for seq in range(40):
                f.write(">%d_%d\nMKVLA\n" % (sp, seq))
    ogs = [["%d_%d" % (sp, q) for sp in range(30)] for q in range(40)]
    (tmp_path / "out").mkdir()
    seen = []
    write_range = msa_trees._write_range

    def spy(blocks, *args):
        if blocks:
            seen.append(sorted(os.listdir(os.path.dirname(blocks[0][3]))))
        return write_range(blocks, *args)
    monkeypatch.setattr(msa_trees, "_write_range", spy)
    msa_trees.write_orthogroup_fasta_files([str(wd) + "/"], list(range(30)), ogs, lambda i: str(tmp_path / "out" / ("OG%d.fa" % i)),
                                           buffer_bytes=2000, tmp_dir=str(tmp_path))
    assert len(seen) >= 5                                  # (several ranges, of 30 species' files)
    assert all(len(files) == 1 and files[0].startswith("pack_") for files in seen)
    assert len(os.listdir(tmp_path / "out")) == 40
