import itertools
import os
import pathlib
import re
import random

import pandas as pd

def grouper(iterable, n, fillvalue=None):
  "Collect data into fixed-length chunks or blocks"
  # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
  args = [iter(iterable)] * n
  return itertools.zip_longest(*args, fillvalue=fillvalue)


def write_batch_fst_computation(group_label, pred_fids, temp_dir):
  """
  Batch FST computations
  """
  formatted_fids = ' '.join([str(i) for i in pred_fids])
  script_path = temp_dir.joinpath(f'run_{group_label}.sh')
  script_string = (
    "#!/bin/bash\n"
    f"#SBATCH -J fst_{group_label}\n"
    f"#SBATCH -o fst_{group_label}.o%j\n"
    f"#SBATCH -e fst_{group_label}.o%j\n"
    "#SBATCH -p normal\n"
    "#SBATCH -N 1\n"
    "#SBATCH -n 4\n"
    "#SBATCH -t 36:00:00\n"
    "#SBATCH -A OTH21148\n"
    "#SBATCH --mail-user=joyce.wang@utexas.edu\n"
    "#SBATCH --mail-type=begin\n"
    "#SBATCH --mail-type=end\n\n"
    "plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'\n"
    "scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'\n"
    f"group={group_label}\n"
    f"for pred_fid in {formatted_fids}\n"
    "do\n"
    "  $plink \\\n"
    "    --bim ${scratch}/ukb_merged/ukb_merged.bim \\\n"
    "    --bed ${scratch}/ukb_merged/ukb_merged.bed \\\n"
    "    --fam ${scratch}/ukb_merged/ukb_merged_edited.fam \\\n"
    "    --family \\\n"
    "    --fst \\\n"
    "    --keep-cluster-names GWAS ${pred_fid} \\\n"
    "    --out ../data/fst/fst_${pred_fid}\n\n"
    "  rm ../data/fst/fst_${pred_fid}.fst\n"
    "  echo ${pred_fid} >> ../data/fst/fst_${group}.est\n"
    "  grep 'Fst estimate:' ../data/fst/fst_${pred_fid}.log >> ../data/fst/fst_${group}.est\n"
    "  rm ../data/fst/fst_${pred_fid}.log\n"
    "done\n"
    )
  with open(script_path, 'w') as f:
    f.write(script_string)

def main():
  temp_dir = pathlib.Path('temp_fst_path/')
  temp_dir.mkdir(exist_ok=True)

  my_file = open("data/ukb_populations/pc_futher_fst_id.txt", "r")
  all_pred = my_file.read()
  all_pred = all_pred.split("\n")
  all_pred.pop()
  my_file.close()

  group_size = 200

  for i in range(len(all_pred)):
    all_pred[i] = "pred" + str(all_pred[i])

  for i, pred_fid_group in enumerate(grouper(all_pred, group_size, '')):
    write_batch_fst_computation(i, pred_fid_group, temp_dir)

if __name__ == "__main__":
    main()
