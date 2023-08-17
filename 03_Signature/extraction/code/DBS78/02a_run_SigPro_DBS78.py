# This script is designed to run in a singularity/apptainer container,
# it is suggested to use a bash script to call this Python script

from pathlib import Path
from SigProfilerExtractor import sigpro as sig
import pandas as pd
import sys

if __name__ == "__main__":
  context_type = "DBS78"
  min_K = 2
  max_K = 5
  nmf_replicates = 1000
  random_seed = 1234  
  num_cpu = 100
  
  project_dir = sys.argv[1]
  
  output_dir = project_dir+"/03_Signature/extraction/raw_output/DBS78/SigProfilerExtractor/"
  path = Path(output_dir)
  if path.exists()==False:
      path.mkdir(parents=True, exist_ok=True)
  
  replicates = list(range(1, nmf_replicates+1))
  seed = [random_seed] * nmf_replicates
  seed_df = pd.DataFrame(list(zip(replicates, seed)), columns=["Replicates","Seeds"])
  seed_df = seed_df.set_index("Replicates")
  seed_df.to_csv(output_dir+"/Seeds.txt", sep="\t")
  
  input_catalog = project_dir+"/03_Signature/extraction/input/DBS78/clca_catalog_DBS78.txt"
  seeds = output_dir+"/Seeds.txt"
  
  sig.sigProfilerExtractor(
    input_type="matrix",
    output=output_dir,
    input_data=input_catalog,
    reference_genome="GRCh37",
    context_type=context_type,
    exome=False,
    minimum_signatures=min_K,
    maximum_signatures=max_K,
    nmf_replicates=nmf_replicates,
    resample=True,
    seeds=seeds,
    matrix_normalization="gmm",
    nmf_init="random",
    precision="single",
    cpu=num_cpu,
  )
