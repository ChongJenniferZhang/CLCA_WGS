This folder contains code and results of mutational signature analysis in the paper: Deep Whole Genome Analysis of 494 Hepatocellular Carcinomas.

### Running environment
The software needed is installed inside the singularity container `signature.sif`.

Download this container from [Zenodo](https://doi.org/10.5281/zenodo.7260221) and put 
it inside the `containers` folder in the top directory of a clone of this repository.
This is because 
the script that calls the container to run the analysis utilizes the folder name structure.

It is recommended to run the code on a Linux server with >= 100 CPU cores.

### Structure of this folder
#### `common_code`: Code needed for both mutational signature extraction and assignment analysis.


#### `extraction`: Mutational signature extraction.
- `input`: Mutational catalog of CLCA cohort for SBS96, DBS78 and ID83.
- `code`: R and Python code to run mutational signature extraction software.
Based on users' computing environment, they can choose either one of the following options to run the analysis:
  - `***.pbs`: Job submission script for the [Altair PBS Professional job scheduler](https://altair.com/pbs-professional/) on HPC clusters to call the singularity container. These scripts will need to be modified for other job schedulers.
  - `***.sh`: Bash script to call the singularity container on a Linux server.
- `raw_output`: Selected important raw outputs from mutational signature extraction software.
- `summary`: Mutational signatures extracted from CLCA cohort for SBS96, DBS78 and ID83.


#### `attribution`: Mutational signature assignment.
- `input`: Mutational catalog and signatures of CLCA cohort for SBS96, DBS78 and ID83.
- `code`: R code to run mutational signature assignment software. Based on users' computing environment, they can choose either one of the following options to run the analysis:
  - `***.pbs`: Job submission script for the [Altair PBS Professional job scheduler](https://altair.com/pbs-professional/) on HPC clusters to call the singularity container. These scripts will need to be modified for other job schedulers.
  - `***.sh`: Bash script to call the singularity container on a Linux server.
- `raw_output`: Selected important raw outputs from mutational signature assignment software.
- `summary`: Exposure matrix of mutational signatures in CLCA cohort for SBS96, DBS78 and ID83.
