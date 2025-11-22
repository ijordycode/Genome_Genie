#!/bin/bash

# Run the samplesheet script
python3 Scripts/samplesheet_generator.py

# Run the snakemake file
snakemake -j 2 --configfile config.yaml -n 