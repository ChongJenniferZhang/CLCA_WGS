from pathlib import Path
from SigProfilerExtractor import sigpro as sig
import pandas as pd
import sys
import os

if __name__ == "__main__":
  context_type = "ID83"
  cwd = os.getcwd()
  
  output_dir = cwd+"/example_output/SigPro/"
  path = Path(output_dir)
  if path.exists()==False:
      path.mkdir(parents=True, exist_ok=True)
  
  input_catalog = cwd+"/example_input/catalog.txt"
  
  sig.sigProfilerExtractor(
    input_type="matrix",
    output=output_dir,
    input_data=input_catalog,
    reference_genome="GRCh37",
    context_type=context_type,
    exome=False,
    minimum_signatures=2,
    maximum_signatures=4,
    nmf_replicates=5,
    resample=True,
    matrix_normalization="gmm",
    nmf_init="random",
    precision="single",
    cpu=10,
  )
