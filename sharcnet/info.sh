#!/bin/bash

# Path variables
project_dir=$(realpath "$project_dir")
sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

input_dir="$project_dir/input"
step_00_dir="$project_dir/00_oneline"
step_01_dir="$project_dir/01_windows"
step_02_dir="$project_dir/02_clusters"
step_03_dir="$project_dir/03_screen"
step_04_dir="$project_dir/04_locarna"
step_05_dir="$project_dir/05_evaluation"

logs_dir="$project_dir/logs"

container_dir="$sharcnet_dir/containers"
rnatools="$container_dir/rnatools_v2.1.sif"

## Pipeline parameters
window_size=250
overlap=75

fold_temperature=21
