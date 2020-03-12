import os
import pandas as pd
from pathlib import Path
import re
from snakemake.utils import validate, min_version

min_version('5.10.0')

configfile: 'config.yaml'
validate(config, 'schemas/config.schema.yaml')

read_metadata = pd.read_table(config['read-metadata']) \
    .set_index('filename', drop=False)
validate(read_metadata, 'schemas/read_metadata.schema.yaml')

wildcard_constraints:
    iteration=r'\d*'

localrules: all, cluster_config, link_contigs, bam_fofn, index_fasta,
    fasta_slices, fasta_slice_fofn, aggregate_pilon,
    genome_wide_coverage_plot

rule all:
    input: 'data/contigs_pilon_{iterations}.fasta' \
        .format(iterations=config['iterations'])

include: 'rules/data_management.smk'
include: 'rules/polishing.smk'
include: 'rules/alignment_stats.smk'
include: 'rules/alignment.smk'
include: 'rules/cluster_config.smk'
