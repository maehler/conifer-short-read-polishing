def get_pilon_slices(wildcards):
    output_dir = checkpoints.fasta_slices.get(iteration=int(wildcards.iteration)-1).output[0]
    fname_pattern = '{output_dir}/contigs_pilon_{iteration}_{{slice}}' \
        .format(output_dir=output_dir, iteration=int(wildcards.iteration)-1)
    gwc = glob_wildcards(fname_pattern)
    fnames = expand('results/pilon_{iteration}/polished_slices/contigs_pilon_{iteration}_{fix}_{slice}.fasta',
        iteration=wildcards.iteration, fix=wildcards.fix, slice=gwc.slice)
    return sorted(fnames, key=lambda x: int(re.search(r'(\d+)\.fasta$', x).group(1)))

rule aggregate_pilon:
    input: get_pilon_slices
    output:
        fasta=protected('results/pilon_{iteration}/contigs_pilon_{iteration}_{fix}.fasta'),
        linked_fasta='data/contigs_pilon_{iteration}_{fix}.fasta'
    shell:
        '''
        echo {input} | xargs -n100 cat > {output.fasta}
        cd data
        ln -s ../{output.fasta} $(basename {output.linked_fasta})
        '''

rule pilon:
    '''
    Run one iteration of pilon on a slice of the assembly.
    '''
    input:
        fasta=lambda wildcards: 'data/contigs_pilon_{iteration}.fasta' \
            .format(iteration=int(wildcards.iteration)-1),
        bams=lambda wildcards: 'results/alignments/read_alignments_{iteration}.fofn' \
            .format(iteration=int(wildcards.iteration)-1),
        fasta_slice=lambda wildcards: 'data/contigs_pilon_{iteration}_slices/contigs_pilon_{iteration}_{{slice}}' \
            .format(iteration=int(wildcards.iteration)-1)
    output:
        temp('results/pilon_{iteration}/polished_slices/contigs_pilon_{iteration}_{fix}_{slice}.fasta')
    params:
        strays=lambda wildcards: '--nostrays' \
            if wildcards.fix == 'snps,indels' or wildcards.fix == 'bases' \
            else ''
    threads: 10
    conda: '../envs/pilon.yaml'
    shell:
        '''
        bam_arg=$(cat {input.bams} | xargs -n1 -I{{}} echo '--frags {{}}')
        pilon \\
            -Xms59G \\
            -Xmx59G \\
            --threads {threads} \\
            --fix {wildcards.fix} {params.strays} \\
            --genome {input.fasta} \\
            ${{bam_arg}} \\
            --targets {input.fasta_slice} \\
            --output results/pilon_{wildcards.iteration}/polished_slices/contigs_pilon_{wildcards.iteration}_{wildcards.fix}_{wildcards.slice}
        '''

checkpoint fasta_slices:
    '''
    Split the contigs. In order to save space, the slices only
    contain the names of the contigs, and not the sequences.
    '''
    input: 'data/contigs_pilon_{iteration}.fasta.fai'
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
                        print('\n'.join(current_slice), file=of)
                    current_slice = []
                    current_slice_size = 0
                    current_slice_idx += 1
            if current_slice_size > 0:
                with open('{dir}/contigs_pilon_{iteration}_{slice_idx}' \
                        .format(dir=outdir, iteration=wildcards.iteration, slice_idx=current_slice_idx), 'w') as of:
                    print('\n'.join(current_slice), file=of)

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
    conda: '../envs/samtools.yaml'
    shell: 'samtools faidx {input}'
