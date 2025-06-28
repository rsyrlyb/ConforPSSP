#!/bin/bash

nodes=1
cores=1
queue="gpu"
walltime="7-00:00:00"
mem="20G"
gpu=1

dir_pn=/home/raulia/software/miniconda3/envs/conforpssp/

sbatch << EOF
#!/bin/bash
#SBATCH -p $queue
#SBATCH -c $cores
#SBATCH --mem=$mem
#SBATCH --nodes=1
#SBATCH --gres=gpu:$gpu
#SBATCH -o "slurm"_%j
#SBATCH --array=1-1
#SBATCH -J mlfold

source activate $dir_pn
export NCCL_P2P_DISABLE=1
unset TMPDIR

python run_confor-pssp.py --fasta_file $1 --output_dir ./ --model 3

EOF
