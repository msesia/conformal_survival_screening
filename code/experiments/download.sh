DATA=$1

mkdir -p results_hpc

rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/conformal_survival_private/code/experiments/results/* results_hpc/

#rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/conformal_survival_private/code/experiments/figures/* results_hpc/figures/
