def get_pilon_logs(wildcards):
    fname_pattern = 'results/pilon_{iteration}/polished_slices/contigs_pilon_{iteration}_{{slice}}.log' \
        .format(iteration=wildcards.iteration)
    gwc = glob_wildcards(fname_pattern)
    return expand(fname_pattern, slice=gwc.slice)

rule pilon_stats:
    input: get_pilon_logs
    output: 'results/pilon_{iteration}/stats.tsv'
    script: '../scripts/parse_pilon_log.py'
