
#vclust prefilter -i PIGEON3.0.fa  -o pigeon3_prefilteredforvclust.fasta

vclust align -i PIGEON3.0.fa  -o pigeon3_vclust_aligned.fasta --filter pigeon3_prefilteredforvclust.fasta --outfmt lite

vclust cluster -i pigeon3_vclust_aligned.fasta -o pigeon3_votus_clustered.tsv --ids pigeon3_vclust_aligned.ids.fasta  --algorithm leiden --metric ani --ani 0.95 --qcov 0.85
