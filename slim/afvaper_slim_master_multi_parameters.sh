#!/bin/bash
#SBATCH -D .
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=sim_run
#SBATCH --error=logs/%x_%j_%a
#SBATCH --output=logs/%x_%j_%a
#SBATCH --export=All
#SBATCH --mail-type=END
#SBATCH --mail-user=j.whiting2@exeter.ac.uk
#SBATCH --array=1-27
#SBATCH --mem=64G

# This script handles the running and processing of parallel selective sweep runs

#########################################################################################################
# Load modules we want
module load vcftools
module load bcftools/1.11
module load r/4.1.0
module load slim
#module load Python/3.7.2-GCCcore-8.2.0
source ~/afvaper_env/bin/activate

#########################################################################################################
# Set our environment
MAIN=/home/jimw91/workdir/afvaper_slim

# Get parameters from the metadata
METADATA=$MAIN/slim/sim_parameter_metadata_v3.txt

ITERATIONS=100
SELECTION=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | cut -f1)
POPN=10000
MUT_RATE=0.00000005
REC_RATE=0.00000001
REC_MAP=$MAIN/data/variable_recombination.map
CORES=8
DAUGHTER_POP_SIZE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | cut -f2)
MIG_RATE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | cut -f3)
EVOLVING_GENERATIONS=200
EIGEN_WINDOW_SIZE=(50 200 500)
EIGEN_PERMUTE=1000
RUN_NAME="21,09,07,Running_a_full_divergent_example"
RUN_TYPE="FULL_DIVERGENT"

# Make an output dir for this run
mkdir $MAIN/outputs/slim
mkdir $MAIN/outputs/slim/$RUN_NAME
cd $MAIN/outputs/slim/$RUN_NAME

# Make a log file with run info
echo "run_type=${RUN_TYPE}" > ${RUN_NAME}.settings.log
echo "selection=${SELECTION}" >> ${RUN_NAME}.settings.log
echo "popN=${POPN}" >> ${RUN_NAME}.settings.log
echo "theta=${THETA}" >> ${RUN_NAME}.settings.log
echo "rho=${RHO}" >> ${RUN_NAME}.settings.log
for eigen in "${EIGEN_WINDOW_SIZE[@]}"
do
echo "eigen_window_size=${eigen}" >> ${RUN_NAME}.settings.log
done
echo "daughter_pop_size=${DAUGHTER_POP_SIZE}" >> ${RUN_NAME}.settings.log
echo "mig_rate=${MIG_RATE}" >> ${RUN_NAME}.settings.log
#echo "final_freq=${FINAL_FREQ}" >> ${RUN_NAME}.settings.log
echo "evolving_generations=${EVOLVING_GENERATIONS}" >> ${RUN_NAME}.settings.log

#########################################################################################################
# First, run simulations and process outputs to filtered VCFs
echo "STARTING SIMULATAION STAGE"

# Run some runs of our parallel simulation - remember to set run_type
if [ "$RUN_TYPE" = "FULL_PARALLEL" ]
then
seq 1 $ITERATIONS | parallel -j $CORES "slim -t -d rec_rate=$REC_RATE -d selection_coef=$SELECTION -d popN=$POPN -d daughter_pop_size=$DAUGHTER_POP_SIZE -d mig_rate=$MIG_RATE -d evolving_generations=$EVOLVING_GENERATIONS $MAIN/slim/SGV_sweeps_all_parallel_trees_fullrecap.slim > run_{}.log"
elif [ "$RUN_TYPE" = "TWO_PARALLEL" ]
then
seq 1 $ITERATIONS | parallel -j $CORES "slim -t -d rec_rate=$REC_RATE -d selection_coef=$SELECTION -d popN=$POPN -d daughter_pop_size=$DAUGHTER_POP_SIZE -d mig_rate=$MIG_RATE -d evolving_generations=$EVOLVING_GENERATIONS $MAIN/slim/SGV_sweeps_two_parallel_lineages_trees_fullrecap.slim > run_{}.log"
elif [ "$RUN_TYPE" = "FULL_DIVERGENT" ]
then
seq 1 $ITERATIONS | parallel -j $CORES "slim -t -d rec_rate=$REC_RATE -d selection_coef=$SELECTION -d popN=$POPN -d daughter_pop_size=$DAUGHTER_POP_SIZE -d mig_rate=$MIG_RATE -d evolving_generations=$EVOLVING_GENERATIONS $MAIN/slim/SGV_sweeps_four_divergent_lineages_fullrecap.slim > run_{}.log"
elif [ "$RUN_TYPE" = "FULL_PARALLEL_RECOMBINATION" ]
then
seq 1 $ITERATIONS | parallel -j $CORES "slim -t -d sweep_locus=999999 -d selection_coef=$SELECTION -d popN=$POPN -d daughter_pop_size=$DAUGHTER_POP_SIZE -d mig_rate=$MIG_RATE -d evolving_generations=$EVOLVING_GENERATIONS $MAIN/slim/SGV_sweeps_all_parallel_trees_full_recap_variable_recombination.slim > run_{}_lowRecomb.log"
seq 1 $ITERATIONS | parallel -j $CORES "slim -t -d sweep_locus=10999999 -d selection_coef=$SELECTION -d popN=$POPN -d daughter_pop_size=$DAUGHTER_POP_SIZE -d mig_rate=$MIG_RATE -d evolving_generations=$EVOLVING_GENERATIONS $MAIN/slim/SGV_sweeps_all_parallel_trees_full_recap_variable_recombination.slim > run_{}_medRecomb.log"
seq 1 $ITERATIONS | parallel -j $CORES "slim -t -d sweep_locus=20999999 -d selection_coef=$SELECTION -d popN=$POPN -d daughter_pop_size=$DAUGHTER_POP_SIZE -d mig_rate=$MIG_RATE -d evolving_generations=$EVOLVING_GENERATIONS $MAIN/slim/SGV_sweeps_all_parallel_trees_full_recap_variable_recombination.slim > run_{}_highRecomb.log"
fi

# Clean up burns and established...
rm -f *burn* *daughter* *established*

# Now convert each of our final simulation states to VCF
if [ "$RUN_TYPE" = "FULL_PARALLEL_RECOMBINATION" ]
  then
    ls *.trees | parallel -j $CORES "python $MAIN/slim/mutate_parallel_sweep_trees_cmd.py {} $MUT_RATE $POPN 0"
  else
    ls *.trees | parallel -j $CORES "python $MAIN/slim/mutate_parallel_sweep_trees_cmd.py {} $MUT_RATE $POPN $REC_RATE"
fi

# Now MAF filter and invariant filter our VCFs
ls *.vcf | parallel -j $CORES "bcftools view -c 1 -c 1:nonmajor {} > {}_hets.vcf"
ls *hets.vcf | parallel -j $CORES "bcftools view -i 'MAF > 0.01' {} > {}_maf.vcf"

# Tidy up
rm -f *trees.vcf *hets.vcf

#########################################################################################################
# Now we'll do our calculations of windowed Fst with vcftools
echo "STARTING FST STAGE"
# Make all our popfiles
rm -f simple.popmap
for i in $(seq 1 5)
do
  for j in $(seq 1 20)
  do
    echo "pop${i}" >> simple.popmap
  done
done

# Merge with IDs
for vcf in $(ls *maf.vcf)
do
  bcftools query -l $vcf | paste - simple.popmap > $vcf.popmap
done

# Use vcftools to get Fst...
for vcf in $(ls *maf.vcf)
  do
  for i in {2..5}
  do
    grep "pop1" $vcf.popmap | cut -f1 > pop1.popmap
    grep "pop${i}" $vcf.popmap | cut -f1 > pop${i}.popmap
    vcftools --vcf $vcf --weir-fst-pop pop1.popmap --weir-fst-pop pop${i}.popmap --fst-window-size 10000 --out ${vcf}_pop1_vs_pop${i}
  done
done
#########################################################################################################
# Finally perform Eigen analysis and save results to a rds
echo "STARTING EIGEN STAGE"

# This is parallelised internally so we just loop
for vcf in $(ls *maf.vcf)
do
  for wind_size in "${EIGEN_WINDOW_SIZE[@]}"
  do
  Rscript $MAIN/slim/afvape_slim_vcf.R $vcf $wind_size $CORES $EIGEN_PERMUTE
  done
done

# Move final results
mkdir final_results
mv *fst final_results
mv *rds final_results

# Tidy up dregs
rm -f *finished.trees* *popmap
