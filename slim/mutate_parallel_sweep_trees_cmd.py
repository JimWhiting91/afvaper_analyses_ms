# Keywords: Python, tree-sequence recording, tree sequence recording

# Use this python script to read trees back in, perform recapitation, and add back in neutral mutations

import msprime, pyslim, time, gzip, sys
import numpy as np
import matplotlib.pyplot as plt

# Test variables

# Set theta and other variables
tree_res=str(sys.argv[1])
#theta=float(sys.argv[2])
mut_rate=float(sys.argv[2])
#rho=float(sys.argv[3])
pop_size=float(sys.argv[3])
rec_rate=float(sys.argv[4])

# Make our recombination map
recombmap = msprime.RecombinationMap(positions=[0,2000000,4000000,6000000,8000000,10000000,12000000,14000000,16000000,18000000,20000000,22000000],
                                     rates=[1.00E-09,2.80E-09,4.60E-09,6.40E-09,8.20E-09,1.00E-08,2.80E-08,4.60E-08,6.40E-08,8.20E-08,1.00E-07,0.00E+00])


# Load the .trees file
ts1 = pyslim.load(tree_res)

# If we have a variable recombination map, ie rec=0, use our map
if(rec_rate == 0):
    print("Recapitating with recombination map")
    ts1_recap = ts1.recapitate(recombination_map=recombmap, Ne=pop_size, random_seed=1)
else:
    print("Recapitating with fixed recombination")
    ts1_recap = ts1.recapitate(recombination_rate=rec_rate, Ne=pop_size, random_seed=1)

# Mutate + simplify the output
mutated = msprime.mutate(ts1_recap, rate=mut_rate, random_seed=1, keep=True)
mutated2 =  mutated.simplify()

print(f"The tree sequence has {mutated2.num_trees} trees on a genome of length {mutated2.sequence_length},"
      f" {mutated2.num_individuals} individuals, {mutated2.num_samples} 'sample' genomes,"
      f" and {mutated2.num_mutations} mutations.")

# Convert to pyslim format
ts = pyslim.SlimTreeSequence(mutated2)

# Get all the individuals alive at present (aka, time 0)
inds = ts.individuals_alive_at(0,stage='late',remembered_stage='early')

# Determine all samples that are in a population
samples_of_pop = dict()
for i in inds:
    pop = ts.individual(i).population
    samples_of_pop.setdefault(pop, []).append(i)

# Sample the desired invividuals
# Create new sample IDs
num_samples = 20
sample_names = []
sampled_inds = []
for pop in sorted(samples_of_pop):
    samples = sorted(np.random.choice(samples_of_pop[pop], num_samples, replace=False))
    for sam in samples:
        sampled_inds.append(sam)
        new_sam_name = 'msp_{}'.format(sam)
        sample_names.append(new_sam_name)

# Write VCF
#with open(sys.argv[2], "w") as vcf_file:
with open(''.join([tree_res,".vcf"]), "w") as vcf_file:
    ts.write_vcf(output=vcf_file,
        individuals=sampled_inds,
        individual_names=sample_names,
        contig_id='chr1')
