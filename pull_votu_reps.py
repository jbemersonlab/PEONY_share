import pandas as pd

       # Load data
clusters = pd.read_csv("pigeon3_votus_clustered.tsv", sep="\t")
lengths = pd.read_csv("lengths.tsv", sep="\t", header=None, names=["contig_id", "length"])

        # Merge cluster + length tables
df = clusters.merge(lengths, on="contig_id", how="left")

        # For each vOTU (cluster), pick the longest contig
reps = df.loc[df.groupby("cluster")["length"].idxmax()]

        # Save table with cluster_id + representative contig
reps.to_csv("votu_reps.tsv", sep="\t", index=False)

        # Save plain list of representatives for FASTA extraction
reps["contig_id"].to_csv("votureps_list.tsv", sep="\t", index=False, header=False)

