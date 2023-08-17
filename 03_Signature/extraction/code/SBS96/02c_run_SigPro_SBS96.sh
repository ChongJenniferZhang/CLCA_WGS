#!/bin/bash

# User needs to set the environment variable CLCA_PROJECT_DIR before running the bash script
if [[ ! -v CLCA_PROJECT_DIR ]]; then
    echo "CLCA_PROJECT_DIR should be set to be the root directory containing the CLCA paper code"
    exit 
elif [[ -z "$CLCA_PROJECT_DIR" ]]; then
    echo "CLCA_PROJECT_DIR should be set to be the root directory containing the CLCA paper code"
    exit 
else
    echo "CLCA_PROJECT_DIR has the value: $CLCA_PROJECT_DIR"
fi

export LOG_PATH=$CLCA_PROJECT_DIR/03_Signature/extraction/code/SBS96/02a_run_SigPro_SBS96.log

nice singularity exec -e $CLCA_PROJECT_DIR/containers/signature.sif python3.8 $CLCA_PROJECT_DIR/03_Signature/extraction/code/SBS96/02a_run_SigPro_SBS96.py $CLCA_PROJECT_DIR > $LOG_PATH 2>&1 &