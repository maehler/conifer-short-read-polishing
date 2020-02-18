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

localrules: all, cluster_config, link_contigs, bam_fofn, index_fasta,
    fasta_slices, fasta_slice_fofn

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

checkpoint fasta_slices:
    '''
    Split the contigs. In order to save space, the slices only
    contain the names of the contigs, and not the sequences.
    '''
    input: 'data/contigs_pilon_{iteration}.fa.fai'
    output: directory('data/contigs_pilon_{iteration}_slices')
    run:
        outdir = Path(output[0])
        outdir.mkdir(exist_ok=True)
        current_slice = []
        current_slice_idx = 0
        current_slice_size = 0
        with open(input[0]) as f:
            for line in f:
                line = line.strip().split()
                current_slice.append(line[0])
                current_slice_size += int(line[1])
                if current_slice_size >= config['slice_size']:
                    with open('{dir}/contigs_pilon_{iteration}_{slice_idx}' \
                            .format(dir=outdir, iteration=wildcards.iteration, slice_idx=current_slice_idx), 'w') as of:
                        of.write('\n'.join(current_slice))
                    current_slice = []
                    current_slice_size = 0
                    current_slice_idx += 1
            if current_slice_size > 0:
                with open('{dir}/contigs_pilon_{iteration}_{slice_idx}' \
                        .format(dir=outdir, iteration=wildcards.iteration, slice_idx=current_slice_idx), 'w') as of:
                    of.write('\n'.join(current_slice))

def collect_fasta_slices(wildcards):
    '''
    Collect the fasta slices and return them sorted by their slice index.
    '''
    output_dir = checkpoints.fasta_slices.get(iteration=wildcards.iteration).output[0]
    fname_pattern = '{output_dir}/contigs_pilon_{iteration}_{{slice}}'.format(output_dir=output_dir, iteration=wildcards.iteration)
    gwc = glob_wildcards(fname_pattern)
    fnames = expand(fname_pattern, slice=gwc.slice)
    return sorted(fnames, key=lambda x: int(re.search(r'\d+$', x).group(0)))

rule fasta_slice_fofn:
    '''
    Create a file-of-filenames of the fasta slices.
    '''
    input: collect_fasta_slices
    output: 'data/contigs_pilon_{iteration}_slices.fofn'
    run:
        with open(output[0], 'w') as f:
            f.write('\n'.join(str(Path(x).resolve()) for x in input))

rule index_fasta:
    '''
    Index a fasta file.
    '''
    input: '{fastafile}'
    output: '{fastafile}.fai'
    wildcard_constraints:
        fastafile='.+\.(fa|fasta)(\.gz)?'
    conda: 'envs/samtools.yaml'
    shell: 'samtools faidx {input}'

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
            for fname in input:
                print(Path(fname).resolve(), file=f)
        
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
