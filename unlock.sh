#! /bin/bash -login

## Exit if any command fails
set -e

## Load required modules
module load python/3.6.6

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecificy which workflow to unlock (i.e. RNApipeCore)'
            exit 2
            ;;
    
    'RNApipeCore' | 'RNAcore' | 'rnaCore')
            ## Unlock snakemake workflow
            snakemake -j 100 --unlock -s workflows/RNApipeCore.snakefile --cores 1
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun with sbatch RNApipeCore.sh"