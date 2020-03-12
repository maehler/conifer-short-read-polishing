rule link_contigs:
    '''
    Symlink the contigs that should be used as the base for the polishing.
    '''
    input: Path(config['contig-fasta']).resolve()
    output: 'data/contigs_pilon_0.fasta'
    shell:
        '''
        cd data
        ln -s {input} $(basename {output})
        '''
