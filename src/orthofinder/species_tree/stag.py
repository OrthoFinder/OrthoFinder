#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2018 David Emms
#
# This program (OrthoFinder) is distributed under the terms of the GNU General Public License v3
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# For any enquiries send an email to David Emms
# david_emms@hotmail.com

import os
import sys
import glob

import argparse
import subprocess
import numpy as np
from itertools import combinations
try:
    from rich import print
except ImportError:
    ...

from ..phylo import tree, newick
from ..utils import util, parallel_task_manager
from . import consensus_tree as cons

# import tree, newick, util, parallel_task_manager
# import consensus_tree as cons
        
# def can_run_command(command, qAllowStderr = False, qPrint = True):
#     if qPrint: sys.stdout.write("Test can run \"%s\"" % command)      
#     capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     stdout = [x for x in capture.stdout]
#     stderr = [x for x in capture.stderr]

#     if len(stdout) > 0 and (qAllowStderr or len(stderr) == 0):
#         if qPrint: print(" - ok")
#         return True
#     else:
#         if qPrint: print(" - failed")
#         return False
        
# def CheckFastME(workingDir):
#     testFN = workingDir + "SimpleTest.phy"
#     with open(testFN, 'w') as outfile:
#         outfile.write("4\n1_1 0 0 0.2 0.25\n0_2 0 0 0.21 0.28\n3_1 0.2 0.21 0 0\n4_1 0.25 0.28 0 0")
#     outFN = workingDir + "SimpleTest.tre"
#     if os.path.exists(outFN): os.remove(outFN)  
#     print("")      
#     if not can_run_command("fastme -i %s -o %s" % (testFN, outFN), qAllowStderr=False):
#         print("ERROR: Cannot run fastme")
#         print('Please check "fastme" is installed and that the executable is in the system path.\n')
#         return False
#     os.remove(testFN)
#     os.remove(outFN)
#     fastme_stat_fn = workingDir + "SimpleTest.phy_fastme_stat.txt"
#     if os.path.exists(fastme_stat_fn): os.remove(fastme_stat_fn)

def write_phylip_matrix(m, names, out_fn, max_og=1e6):
    """
    m - list of mp.Array  (so that each represents an nSeq x nSeq matrix
    """
    sliver = 1e-6
    with open(out_fn, 'w') as outfile:
        n = len(m)
        outfile.write("%d\n" % n)
        for i in range(n):
            outfile.write(names[i] + " ")
            # values could be -inf, these are the most distantly related so replace with max_og
            V = [0. + (m[i][j] if m[i][j] > -9e99 else max_og) for j in range(n)] # "0. +": hack to avoid printing out "-0"
            V = [sliver if v < sliver  else v for v in V]  # make sure scientific notation is not used (not accepted by fastme)
            values = " ".join(["%.6f" % v for v in V])   
            outfile.write(values + "\n")

class UnrecognisedGene(Exception):
    pass

class GeneToSpecies(object):
    def __init__(self, gene_map_fn):
        if not os.path.isfile(gene_map_fn):
            print(("ERROR: Could not open %s" % gene_map_fn))
            sys.exit()
        self.exact = dict()
        self.startswith = dict()
        if gene_map_fn.endswith("SpeciesIDs.txt"):
            with open(gene_map_fn, 'r') as infile:
                for line in infile:
                    _, sp = line.rstrip().rsplit(".", 1)[0].split()
                    sp = sp.replace(".", "_")
                    self.startswith[sp] = sp
        else:
            with open(gene_map_fn, 'r') as infile:
                for i, line in enumerate(infile):
                    t = line.rstrip().split()
                    if len(t) != 2: 
                        print(("ERROR: Invalid format in gene to species mapping file line %d" % (i+1))) 
                        print(line)
                        sys.exit()
                    g, sp = t
                    if g.endswith("*"):
                        self.startswith[g[:-1]] = sp
                    else:
                        self.exact[g] = sp
        self.species = sorted(list(set(list(self.exact.values()) + list(self.startswith.values()))))
        print(("%d species in mapping file:" % len(self.species)))
        for s in self.species:
            print(s)
        print("\nSTAG will infer a species tree from each gene tree with all species present") 
        self.sp_to_i = {s:i for i,s in enumerate(self.species)}
                
    def to_species(self, gene):
        if gene in self.exact: return self.exact[gene]
        else:
            for k, v in self.startswith.items():
                if gene.startswith(k): return v
            # if not found then raise an Exceptions
            raise UnrecognisedGene(gene)
                
    def species_to_index_dict(self):
        return self.sp_to_i
    
    def number_of_species(self):
        return len(self.species)

class GeneToSpeciesOrthoFinder(GeneToSpecies):
    def __init__(self, species_to_use):
        # no exact conversions
        self.exact = dict()
        self.startswith = {("%d_" % sp):str(sp) for sp in species_to_use}
        self.species = list(map(str,species_to_use))
        self.sp_to_i = {s:i for i,s in enumerate(self.species)}

def get_distances_fast(t, n_sp, g_to_i):
    D = np.ones((n_sp, n_sp)) * 9e99
    for n in t.traverse('postorder'):
        if n.is_leaf():
            n.add_feature('d', {g_to_i[n.name]:max(0.0, n.dist)})
        else:
            children = n.get_children()
            for ch0, ch1 in combinations(children,2):
                for sp0,dist0 in ch0.d.items():
                    for sp1,dist1 in ch1.d.items():
                        if sp0==sp1: continue
                        i = sp0 if sp0<sp1 else sp1
                        j = sp1 if sp0<sp1 else sp0
                        D[i,j] = min(D[i,j], dist0+dist1)
                spp = {k for ch in children for k in ch.d.keys()}
                d = {k:(min([ch.d[k] for ch in children if k in ch.d])+max(0.0, n.dist)) for k in spp}
                n.add_feature('d', d)
    for i in range(n_sp):
        for j in range(i):
            D[i,j] = D[j, i]
        D[i,i]=0.
    return D

def process_trees(dir_in, dir_matrices, dir_trees_out, GeneToSpecies, q_verbose=True, q_skip_single_copy=False, q_for_of=False):
    n_sp = GeneToSpecies.number_of_species()
    s_to_i = GeneToSpecies.species_to_index_dict()
    n_success = 0
    n_not_all_present = 0
    n_fail = 0
    if q_verbose: print("\nProcessing gene trees:")
    for fn in glob.glob(dir_in + "/*"):
        try:
            t = tree.Tree(fn)
        except newick.NewickError:
            print((os.path.split(fn)[1] + " - WARNING: ETE could not interpret tree file, it will be ignored"))
            n_fail += 1
            continue
        try:
            genes = t.get_leaf_names()
            species = list(map(GeneToSpecies.to_species, genes))
        except UnrecognisedGene as e:
            print((os.path.split(fn)[1] + " - WARNING: unrecognised gene, %s" % str(e)))
            n_fail += 1
            continue
        n_this = len(set(species))
        if n_this != n_sp:
#            print(os.path.split(fn)[1] + " - Only %d species, skipping" % nThis)
            n_not_all_present += 1
            continue
        if q_skip_single_copy and n_this == len(genes):
            # Single copy - don't recalculate the tree
            tree_out_fn = dir_trees_out + os.path.split(fn)[1] + ".tre"
            for n in t:
                n.name = s_to_i[GeneToSpecies.to_species(n.name)]
            t.write(outfile = tree_out_fn, format=5)
            if q_verbose: print((os.path.split(fn)[1] + " - Processed"))
            continue
        g_to_i = {g:s_to_i[s] for g,s in zip(genes, species)}
        D = get_distances_fast(t, n_sp, g_to_i)
        species_names_fastme = list(map(str,range(n_sp)))
        matrix_fn = dir_matrices + os.path.split(fn)[1] + ".dist.phylip"
        tree_out_fn = dir_trees_out + os.path.split(fn)[1] + ".tre"
        write_phylip_matrix(D, species_names_fastme, matrix_fn, max_og=1e6)
        command = "fastme -i %s -o %s -w O -s -n" % (matrix_fn, tree_out_fn)
        if q_for_of:
            parallel_task_manager.run_command(command, q_print_on_error=True, q_print_stderr=True, raise_on_error=True)
        else:
            popen = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            popen.communicate()
        n_success += 1
        if q_verbose: print((os.path.split(fn)[1] + " - Processed"))
    if q_verbose: print(("\nExamined %d trees" % (n_success + n_not_all_present + n_fail)))
    print(("%d trees had all species present and will be used by STAG to infer the species tree" % n_success))
      
def astral(tree_dir, astral_jar_file, q_for_of=False):
    trees_fn = tree_dir + "../TreesFile.txt"
    with open(trees_fn, 'w') as outfile:
        for fn in glob.glob(tree_dir + "/*"):
            t = tree.Tree(fn)    
            outfile.write(t.write(format=9) + "\n")
    species_tree_fn = tree_dir + "../SpeciesTree_ids.txt"
    command = " ".join(["java", "-Xmx6000M", "-jar", astral_jar_file, "-i", trees_fn, "-o", species_tree_fn])
    if q_for_of:
        parallel_task_manager.run_command(command, q_print_on_error=True, q_print_stderr=True, raise_on_error=True)
    else:
        subprocess.call(command, shell=True)
    return tree.Tree(species_tree_fn)     
      
def infer_species_tree(tree_dir, species, output_fn, astral_jar=None):
    if astral_jar == None:
        t = cons.consensus_tree(tree_dir)
    else:
        t = astral(tree_dir, astral_jar)
    for n in t:
        n.name = species[int(n.name)]
    t.write(outfile=output_fn)

def run_for_ortho_finder(dir_in, d_working, species_to_use, species_tree_ids_fn_out):
    dir_matrices = d_working + "Distances_SpeciesTree/"
    os.mkdir(dir_matrices)
    dir_trees_out = d_working + "SpeciesTrees_ids/"
    os.mkdir(dir_trees_out)
    gene_to_species = GeneToSpeciesOrthoFinder(species_to_use)
    process_trees(dir_in, dir_matrices, dir_trees_out, gene_to_species, q_verbose=False, q_for_of=True)
    infer_species_tree(dir_trees_out, gene_to_species.species, species_tree_ids_fn_out)
    
# def main(args):
#     dir_in = args.gene_trees
#     astral_jar = None
#     gene_to_species = GeneToSpecies(args.species_map)
#     dir_out = paths.create_new_working_directory(dir_in + "/../STAG_Results")
#     CheckFastME(dir_out)
#     dir_matrices = dir_out + "DistanceMatrices/"
#     os.mkdir(dir_matrices)
#     dir_trees_out = dir_out + "Trees/"
#     os.mkdir(dir_trees_out)
#     process_trees(dir_in, dir_matrices, dir_trees_out, gene_to_species, qVerbose=(not args.quiet))
#     outputFN = dir_out + "SpeciesTree.tre"
# #    if args.astral_jar == None:
#     if astral_jar == None:
#         infer_species_tree(dir_trees_out, gene_to_species.species, outputFN)
#     else:
#         infer_species_tree(dir_trees_out, gene_to_species.species, outputFN, astral_jar=astral_jar)
#     print(("STAG species tree: " + os.path.abspath(outputFN) + "\n"))

# if __name__ == "__main__":
#     text = """
# *********************************************************
# *                                                       *
# *      STAG: Species Tree inference from All Genes      *
# *                                                       *
# *********************************************************"""
#     print((text[1:] + "\n"))
#     parser = argparse.ArgumentParser()
#     parser.add_argument("species_map", help = "Map file from gene names to species names, or SpeciesIDs.txt file from OrthoFinder")
#     parser.add_argument("gene_trees", help = "Directory conaining gene trees")
#     parser.add_argument("-q", "--quiet", help = "Only print sparse output", action="store_true")
# #    parser.add_argument("-a", "--astral_jar", help = "ASTRAL jar file. Use ASTRAL to combine STAG species tree estimates instead of greedy consensus tree.")
#     args = parser.parse_args()
#     main(args)