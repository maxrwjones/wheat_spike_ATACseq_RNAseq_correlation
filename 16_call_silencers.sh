#!/bin/bash
#SBATCH --mem 32G
#SBATCH --partition=jic-short,nbi-short
#SBATCH --time=0-02:00
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=../logs/%x_%A.out
#SBATCH --error=../logs/%x_%A.err

source /jic/software/staging/RCSUPPORT-2912/stagingloader

### Change working directory
cd ../outputs

Rscript call_silencers.r

### Create directory for simulated silencers used in step 18
### (Better to do this now because step 18 involves array jobs.)
mkdir silencer_simulations