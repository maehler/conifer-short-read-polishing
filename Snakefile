import os
import pandas as pd
from pathlib import Path
import re
from snakemake.utils import validate

configfile: 'config.yaml'
validate(config, 'schemas/config.schema.yaml')

read_metadata = pd.read_table(config['read_metadata']) \
    .set_index('filename', drop=False)
validate(read_metadata, 'schemas/read_metadata.schema.yaml')

wildcard_constraints:
    iteration=r'\d*'

localrules: all, cluster_config, link_contigs, bam_fofn

rule all:
    input:
        'results/pilon/contigs_pilon_{iterations}.fa'.format(iterations=config['iterations'])

rule link_contigs:
    '''
    Symlink the contigs that should be used as the base for the polishing.
    '''
    input: Path(config['contig_fasta']).resolve()
    output: 'data/contigs_pilon_0.fa'
    shell:
        '''
        cd data
        ln -s {input} $(basename {output})
        '''

rule pilon:
    '''
    Run one iteration of pilon.

    TODO:
        - Set up a conda environment.
        - Write an actual command.
        - Split the input and run pilon on each slice and merge the results.
    '''
    input:
        fasta=lambda wildcards: 'data/contigs_pilon_{iteration}.fa'.format(iteration=int(wildcards.iteration)-1),
        bams=lambda wildcards: 'results/alignments/read_alignments_{iteration}.fofn'.format(iteration=int(wildcards.iteration)-1)
    output:
        'results/pilon/contigs_pilon_{iteration}.fa',
        'data/contigs_pilon_{iteration}.fa'

rule bam_fofn:
    '''
    Create a file-of-filenames of short-read alignments.
    '''
    input:
        expand('results/alignments/read_alignments_{{iteration}}/{bam_name}_alignments.bam',
            bam_name=[Path(x).stem for x in read_metadata.filename])
    output: 'results/alignments/read_alignments_{iteration}.fofn'
    run:
        with open(output[0], 'w') as f:
            for fname in input[0]:
                f.write(fname)
        
def get_read_filename(wildcards):
    '''
    Get the filename(s) for a particular set of reads. If the read data is paired,
    two file names are returned.
    '''
    matches = [re.search(wildcards.bam_name, Path(x).stem) for x in read_metadata.filename]
    rows = [i for i, x in enumerate(matches) if x is not None]

    if len(rows) > 1:
        raise ValueError('more than one entry in the metadata matches {bam_name}' \
            .format(bam_name=bam_name))

    files = []

    files.append(read_metadata.filename[rows[0]])
    if read_metadata.pair[rows[0]]:
        files.append(read_metadata.pair[rows[0]])

    return files

rule bwa_align:
    '''
    Align short reads against the contigs using bwa mem.

    TODO:
        - Make sure appropriate output is created.
    '''
    input:
        reads=get_read_filename,
        index=expand('reference/contigs_pilon_{{iteration}}.{ext}',
            ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        bam='results/alignments/read_alignments_{iteration}/{bam_name}_alignments.bam',
        bai='results/alignments/read_alignments_{iteration}/{bam_name}_alignments.bam.bai'
    threads: 20
    params:
        index_prefix='reference/contigs_pilon_{iteration}'
    conda: 'envs/bwa.yaml'
    shell:
        '''
        bwa mem -t {threads} -o {output.bam}.sam {params.index_prefix} {input.reads}
        samtools view -bu {output.bam}.sam | \\
            samtools sort -m 5G -@ $(({threads} - 1)) -o {output.bam}
        samtools index -@ {threads} {output.bam}
        rm {output.bam}.sam
        '''

rule bwa_index:
    '''
    Create a bwa index.
    '''
    input:
        fasta='data/contigs_pilon_{iteration}.fa'
    output:
        expand('reference/contigs_pilon_{{iteration}}.{ext}',
            ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix='reference/contigs_pilon_{iteration}'
    conda: 'envs/bwa.yaml'
    shell:
        '''
        bwa index {input.fasta} -p {params.prefix}
        '''

rule cluster_config:
    '''
    Generate a cluster profile.
    '''
    output: directory( \
        '{home}/.config/snakemake/{profile_name}' \
            .format(home=Path.home(), \
                    profile_name=config['cluster']['profile_name']))
    params:
        url=config['cluster']['cookiecutter_url'],
        profile_name=config['cluster']['profile_name']
    conda: 'envs/cluster_config.yaml'
    shell: 'bash scripts/cluster_config.sh {params.url} {params.profile_name}'
