include: "common.smk"
include: "identify.smk"
include: "interactions.smk"
include: "map_to_genome.smk"
include: "map_to_viewpoints.smk"
include: "peaks.smk"
include: "pileup.smk"
include: "compare.smk"
include: "stats.smk"

ruleorder: bam_sort_viewpoints > bam_move_to_final_location
ruleorder: bam_index_viewpoints > bam_move_to_final_location
