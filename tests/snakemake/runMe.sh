#!/bin/bash

snakemake --snakefile Snakefile --configfile tests/config/simulation.yml --scheduler greedy --cores 1
