rule bam_fofn:
    """
    Create a file-of-filenames of short-read alignments.
    """
    input:
        expand('results/alignments/read_alignments_{{iteration}}/{bam_name}_alignments.bam',
            bam_name=[Path(x).stem for x in read_metadata.filename])
    output: 'results/alignments/read_alignments_{iteration}.fofn'
    run:
        with open(output[0], 'w') as f:
            for fname in input:
                print(Path(fname).resolve(), file=f)

def get_read_filename(wildcards):
    """
    Get the filename(s) for a particular set of reads. If the read data is paired,
    two file names are returned.
    """
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
    """
    Align short reads against the contigs using bwa mem.
    """
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
    conda: '../envs/bwa.yaml'
    envmodules: 'bioinfo-tools', 'bwa/0.7.17', 'samtools/1.10'
    shell:
        """
        bwa mem -t {threads} -o {output.bam}.sam {params.index_prefix} {input.reads}
        samtools view -bu {output.bam}.sam | \\
            samtools sort -m 5G -@ $(({threads} - 1)) -o {output.bam}
        samtools index -@ {threads} {output.bam}
        rm {output.bam}.sam
        """

rule bwa_index:
    """
    Create a bwa index.
    """
    input:
        fasta='data/contigs_pilon_{iteration}.fasta'
    output:
        expand('reference/contigs_pilon_{{iteration}}.{ext}',
            ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix='reference/contigs_pilon_{iteration}'
    conda: '../envs/bwa.yaml'
    envmodules: 'bioinfo-tools', 'bwa/0.7.17'
    shell:
        """
        bwa index {input.fasta} -p {params.prefix}
        """
