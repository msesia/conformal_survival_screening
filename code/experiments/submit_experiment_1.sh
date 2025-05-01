#!/bin/bash

# Parameters
SETUP=6

if [[ $SETUP == 1 ]]; then
  # Data distribution setting
  SETTING_LIST=(1 2 3)
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf")
  # List of training sample sizes
  N_TRAIN_LIST=(100 200 500 1000 2000 5000 10000)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(10000)
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of maximum number of features to use when fitting censoring model
  N_FEAT_CENS_LIST=(10)
  # List of distribution shift parameters for training data (between 0 and 1)
  SHIFT_LIST=(0)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G

elif [[ $SETUP == 2 ]]; then
  # Data distribution setting
  SETTING_LIST=(4)
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf")
  # List of training sample sizes
  N_TRAIN_LIST=(5000)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(10000)
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of maximum number of features to use when fitting censoring model
  N_FEAT_CENS_LIST=(10)
  # List of distribution shift parameters for training data (between 0 and 1)
  SHIFT_LIST=(0 1)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G

elif [[ $SETUP == 3 ]]; then
  # Data distribution setting
  SETTING_LIST=(1 2 3)
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf")
  # List of training sample sizes
  N_TRAIN_LIST=(5000)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(10000)
  # List of calibration sample sizes
  N_CAL_LIST=(10 20 50 100 200 500 1000)
  # List of maximum number of features to use when fitting censoring model
  N_FEAT_CENS_LIST=(10)
  # List of distribution shift parameters for training data (between 0 and 1)
  SHIFT_LIST=(0)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G

elif [[ $SETUP == 4 ]]; then
  # Data distribution setting
  SETTING_LIST=(1 2 3)
  SURV_MODEL_TYPE_LIST=("grf" "rf" "cox" "survreg")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf" "cox")
  # List of training sample sizes
  N_TRAIN_LIST=(100 1000 10000)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(5000)
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of maximum number of features to use when fitting censoring model
  N_FEAT_CENS_LIST=(10)
  # List of distribution shift parameters for training data (between 0 and 1)
  SHIFT_LIST=(0)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G

elif [[ $SETUP == 5 ]]; then
  # Data distribution setting
  SETTING_LIST=(1 2 3)
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf")
  # List of training sample sizes
  N_TRAIN_LIST=(5000)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(100 200 500 1000 2000 5000)
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of maximum number of features to use when fitting censoring model
  N_FEAT_CENS_LIST=(10)
  # List of distribution shift parameters for training data (between 0 and 1)
  SHIFT_LIST=(0)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G

elif [[ $SETUP == 6 ]]; then
  # Data distribution setting
  SETTING_LIST=(1 2 3)
  # Survival model types
  SURV_MODEL_TYPE_LIST=("grf")
  # Censoring model types
  CENS_MODEL_TYPE_LIST=("grf")
  # List of training sample sizes
  N_TRAIN_LIST=(200 500)
  # List of censoring training sample sizes
  N_TRAIN_CENS_LIST=(5000)
  # List of calibration sample sizes
  N_CAL_LIST=(500)
  # List of maximum number of features to use when fitting censoring model
  N_FEAT_CENS_LIST=(10 20 50 100)
  # List of distribution shift parameters for training data (between 0 and 1)
  SHIFT_LIST=(0)
  # Sequence of batches for parallel simulation
  BATCH_LIST=$(seq 1 5)
  # Memory
  MEMO=5G


fi

# Slurm parameters
TIME=00-00:20:00                    # Time required (20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/setup_"$SETUP

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/setup_"$SETUP

# Loop over configurations
for BATCH in $BATCH_LIST; do
  for SETTING in "${SETTING_LIST[@]}"; do
    for N_CAL in "${N_CAL_LIST[@]}"; do
      for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
        for N_TRAIN_CENS in "${N_TRAIN_CENS_LIST[@]}"; do
          for N_FEAT_CENS in "${N_FEAT_CENS_LIST[@]}"; do
            for SURV_MODEL_TYPE in "${SURV_MODEL_TYPE_LIST[@]}"; do
              for CENS_MODEL_TYPE in "${CENS_MODEL_TYPE_LIST[@]}"; do
                for SHIFT in "${SHIFT_LIST[@]}"; do

                  # Generate a unique and interpretable file name based on the input parameters
                  JOBN="setup_${SETUP}/setting${SETTING}_surv_${SURV_MODEL_TYPE}_cens_${CENS_MODEL_TYPE}_train${N_TRAIN}_trainc${N_TRAIN_CENS}_cal${N_CAL}_nfc${N_FEAT_CENS}_shift${SHIFT}_batch${BATCH}.txt"
                  OUT_FILE=$OUT_DIR"/"$JOBN
                  #ls $OUT_FILE
                  COMPLETE=0

                  if [[ -f $OUT_FILE ]]; then
                    COMPLETE=1
                  fi

                  if [[ $COMPLETE -eq 0 ]]; then
                    # R script to be run with command line arguments
                    SCRIPT="./experiment_1.sh $SETUP $SETTING $SURV_MODEL_TYPE $CENS_MODEL_TYPE $N_TRAIN $N_TRAIN_CENS $N_CAL $N_FEAT_CENS $SHIFT $BATCH"

                    # Define job name for this configuration
                    OUTF=$LOGS"/"$JOBN".out"
                    ERRF=$LOGS"/"$JOBN".err"

                    # Assemble slurm order for this job
                    ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" $SCRIPT"

                    # Print order
                    echo $ORD
                    # Submit order
                    $ORD
                    # Run command now
                    #./$SCRIPT
                    
                  fi

                done
              done
            done
          done
        done
      done
    done
  done
done
