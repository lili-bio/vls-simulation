#!/usr/bin/env python
# coding: utf-8
# simulates vls gene family evolution  
# combines birth-and-death model and converted evolution model


import argparse

# Initialize parser object
parser = argparse.ArgumentParser(description = 'vls gene family simulation')

# Add arguments
parser.add_argument('-u', '--mutation', type=float, default=1e-4,
                    help='gene mutation rate, per site per generation (default = 1e-4)')
parser.add_argument('-r', '--recombination', type=float, default=1e-3,
                    help='gene conversion rate, per gene per generation (default = 1e-3)')
parser.add_argument('-l', '--length', type=int, default=100,
                    help='gene length in amino acid (default = 100)')
parser.add_argument('-a', '--alpha', type=float, default=0.1,
                    help='gene gain rate, per gene per generation (default = 0.1)')
parser.add_argument('-b', '--beta', type=float, default=0.1,
                    help='gene loss rate, per gene per generation (default = 0.1)')
parser.add_argument('-n', '--number', type=int, default=10,
                    help='number of gene copies (default = 10)')
parser.add_argument('-L', '--conversionLength', type=int, default=20,
                    help='length of the fragment involved in gene conversion (default = 20)')
parser.add_argument('-t', '--tree', type=str, default="Bb_genome_reduced.tree",
                    help='the tree file in newick format (default = Bb_genome_reduced.tree)')
parser.add_argument('-s', '--sweep', type=bool, default=False,
                    help='whether perform selective sweep in the conversion events (default = False)')
parser.add_argument('-o', '--output', type=str, default="bb_vls_sim.fas",
                    help='the output alignment file (default = bb_vls_sim.fas)')

args = parser.parse_args()

import numpy as np
import math
import random
from numpy.random import default_rng
rng = default_rng() # random number generator
import copy
from Bio import Phylo

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child.name] = clade
    return parents

def random_protein_sequence(l, n):
    AA = ('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')
    copies = {}
    copy = ''.join(np.random.choice(AA) for _ in range(l))
    for i in range(n):
        gene_id = 'gene%i' %i
        copies[gene_id] = copy   
    return copies 

def mutation(genes, prob, l):
    for name, seq in genes.items(): # for each gene,
        exp = prob * l
        n_mut =  rng.poisson(exp) # determine how many sites to be mutated
        print("n_mut= ", n_mut)
        if n_mut > 0:
            seq = list(seq) # convert string to list in order to subscribe
            sites_mut = random.sample(range(l), n_mut) # choose the sites to mutate
            for site in sites_mut: # mutate the sites one by one
                AA = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
                ori = seq[site]
                AA.remove(ori) # remove the original letter to avoid sampling the same letter
                seq[site] = random.choice(AA) # random substitution
            genes[name] = ''.join(seq)
    return(genes) 

def conversion_sweep(genes, prob, l, L, n):
    exp = prob * n 
    n_con = rng.poisson(exp) # determine how many times of conversion will happen
    print("n_con= ", n_con)
    for j in range(n_con):
        names = list(genes)
        doner = random.choice(names) # pick a gene as doner
        names.remove(doner) # remove the doner to avoid getting the same gene
        start = random.choice(range(l-L)) # pick a start point on the sequence within 1-80
        end = start + L # conversion fragment is 20 aa long
        for recipient in names: # copy the fragment to every gene
            seq = list(genes[recipient])
            seq[start:end] = list(genes[doner])[start:end]
            genes[recipient] = ''.join(seq)
    return(genes)

def conversion_random(genes, prob, l, L, n):
    exp = prob * n 
    n_con = rng.poisson(exp) # determine how many times of conversion will happen
    print("n_con= ", n_con)
    for j in range(n_con):
        names = list(genes)
        doner = random.choice(names) # pick a gene as doner
        names.remove(doner) # remove the doner to avoid getting the same gene
        recipient = random.choice(names) # pick a gene as recipient       
        start = random.choice(range(l-L)) # pick a start point on the sequence within 1-80
        end = start + L # conversion fragment is 20 aa long
        seq = list(genes[recipient])
        seq[start:end] = list(genes[doner])[start:end]
        genes[recipient] = ''.join(seq)
    return(genes)

def duplication(genes, prob, n):
    exp = prob * n
    n_gain = rng.poisson(exp) # determine how many genes to be duplicated
    print("n_gain= ", n_gain)
    if n_gain > 0:
        gene_gain = random.choice(list(genes)) # choose a gene to duplicate 
        for j in range(n_gain):
            i = 1
            name = '%s_%i' %(gene_gain, i)
            while name in genes: # check if name is already there
                i += 1 
                name = '%s_%i' %(gene_gain, i) # create a new name
            genes[name] = genes[gene_gain] # add new copy to the dict
    return(genes)

def deletion(genes, prob, n):
    exp = prob * n
    n_loss = rng.poisson(exp) # determine how many genes to be deleted
    if n_loss > 0:
        if n_loss > (len(genes) - 2): 
            n_loss = len(genes) - 2 # extinction protection
        print("n_loss= ", n_loss)    
        gene_loss = random.sample(list(genes), n_loss) # choose the genes to delete
        for item in gene_loss: # remove items from the dict one by one
            genes.pop(item)
    return(genes)

import textwrap
def write_fasta(dictionary, filename):
    with open(filename, "w") as outfile:
        for key, value in dictionary.items():
            outfile.write(">" + key + "\n")
            outfile.write("\n".join(textwrap.wrap(value, 60)))
            outfile.write("\n")
            

### Main body of the simulator ###

# read the tree and get the clades
tree = Phylo.read(args.tree, "newick") # genome tree of B.burgdorferi
internal_clades = tree.get_nonterminals()
for i in range(len(internal_clades)):
    internal_clades[i].name = "node%i" %i # name internal nodes
#Phylo.draw(tree)

# initiation 
root = random_protein_sequence(args.length, args.number)
tree_dic = {}
tree_dic['node0'] = root # a nesting dictionary

# tree_walking and evolution
allclades = list(tree.find_clades(order="level"))
parents = all_parents(tree)


for clade in allclades[1:]: # except the root 
    time = clade.branch_length/1e-6 # convert branch length to time
    p_mut = 1 - math.exp(-(args.mutation*time)) # probability of mutation
    p_con = 1 - math.exp(-(args.recombination*time))
    p_gain = (args.alpha/(args.alpha+args.beta))*(1-math.exp(-((args.alpha+args.beta)*time))) # probability of gain
    p_loss = 1 - p_gain # probability of loss
    
    parent = parents[clade.name].name # retrieve the name of the parent node
    genes = tree_dic[parent] # retrieve parent's gene set
    genes = mutation(genes, p_mut, args.length) 
    genes = duplication(genes, p_gain, args.number) 
    genes = deletion(genes, p_loss, args.number)
    if args.sweep: 
        gene = conversion_sweep(genes, p_con, args.length, args.conversionLength, args.number)
    else:
        genes = conversion_random(genes, p_con, args.length, args.conversionLength, args.number)
    genes2 = genes.copy()
    tree_dic[clade.name] = genes2 # add new gene set to the dict
            
print("Simulation is done!")    

# extract terminal names and sequences
tips = [tip.name for tip in tree.get_terminals()]
seq_dic = {}
for taxon in tips:
    for name, seq in tree_dic[taxon].items():
        seq_dic["%s_%s" %(taxon, name)] = seq
    
# export simulated sequences to a fasta file
write_fasta(seq_dic, args.output)
print("Success! File written") 

