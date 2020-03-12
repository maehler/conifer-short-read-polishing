rule alignment_window_coverage:
    input:
        fasta_index='data/contigs_pilon_{iteration}.fasta.fai',
        bams='results/alignments/read_alignments_{iteration}.fofn'
    output:
        coverage='results/alignments/read_alignments_{iteration}_window_coverage.tsv',
        genome_size=temp('results/alignments/read_alignments_{iteration}_sizes.txt'),
        windows=temp('results/alignments/read_alignments_{iteration}_windows.bed')
    params:
        window_size=5000,
        step_size=2500
    conda: '../envs/bedtools.yaml'
    envmodules: 'bioinfo-tools', 'BEDTools/2.27.1'
    shell:
        """
        cut -f1,2 {input.fasta_index} > {output.genome_size}

        bedtools makewindows \\
            -g {output.genome_size} \\
            -w {params.window_size} \\
            -s {params.step_size} \\
            > {output.windows}

        while read -r bamfile; do
            bedtools coverage \\
                -a {output.windows} \\
                -b ${{bamfile}} \\
                -sorted \\
                >> {output.coverage}
        done < {input.bams}
        """

rule coverage_histograms:
    input: 'results/alignments/read_alignments_{iteration}.fofn'
    output: 'results/alignments/read_alignments_{iteration}_genomecov.tsv'
    conda: '../envs/bedtools.yaml'
    envmodules: 'bioinfo-tools', 'BEDTools/2.27.1'
    shell:
        """
        while read -r bamfile; do
            bedtools genomecov \\
                -ibam ${{bamfile}} \\
                >> {output}
        done < {input}
        """

rule genome_wide_coverage:
    input: 'results/alignments/read_alignments_{iteration}_genomecov.tsv'
    output: 'results/alignments/read_alignments_{iteration}_global_coverage.tsv'
    shell:
        """
        awk '{{
            coverage[$2] = $3
        }} END {{
            for (depth in coverage) {{
                printf("%d\\t%d\\n", depth, coverage[depth])
            }}
        }}' {input} | \\
        sort -k1,1n \\
        > {output}
        """

rule genome_wide_coverage_plot:
    input: 'results/alignments/read_alignments_{iteration}_global_coverage.tsv'
    output: 'results/alignments/read_alignments_{iteration}_plots/global_coverage_{maxdepth}.png'
    conda: '../envs/r.yaml'
    envmodules: 'R/3.6.1', 'R_packages/3.6.1'
    script: '../scripts/plot_global_coverage.R'

rule flagstats:
    input: 'results/alignments/read_alignments_{iteration}.fofn'
    output: 'results/alignments/read_alignments_{iteration}_flagstats.tsv'
    conda: '../envs/samtools.yaml'
    envmodules: 'bioinfo-tools', 'samtools/1.10'
    shell:
        """
        while read -r bamfile; do
            samtools flagstat -O tsv ${{bamfile}} | \\
                cat <(echo "#" ${{bamfile}}) - \\
                >> {output}
        done < {input}
        """

rule alignment_mapq:
    input: 'results/alignments/read_alignments_{iteration}.fofn'
    output: 'results/alignments/read_alignments_{iteration}_mapq_freq.tsv'
    conda: '../envs/samtools.yaml'
    envmodules: 'bioinfo-tools', 'samtools/1.10'
    shell:
        """
        while read -r bamfile; do
            samtools view ${{bamfile}} | \\
                awk '{{mapqs[$5] += 1; unmapped[$5] += and($2, 4) == 4;}}
                    END {{for (i in mapqs) printf("%d\\t%d\\t%d\\n", i, mapqs[i], unmapped[i])}}' | \\
                sort -k1,1n | \\
                cat <(echo -e '#' ${{bamfile}} '\\n' 'quality\\tcount\\tunmapped') - \\
                >> {output}
        done < {input}
        """
