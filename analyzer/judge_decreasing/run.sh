#!/bin/bash

source ~/.pyenv/versions/miniconda3-3.9.1/envs/my_env/bin/activate network
script_dir=$(cd $(dirname $0); pwd)
python ${script_dir}/judge_decreasing.py $1

