#!/bin/bash

module load Python/3.8.6-GCCcore-10.2.0

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pygeneral

python3 ~/develop/crystalCF/scripts_python/core_program_export.py

