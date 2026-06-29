# Run GeNomad on s contigs
rule GeNomadIndividual:
    input:
        ContigsIn = "{study}.fa"
##        CheckRenameAssemblyContigs = "Checks/3.2-RenameContigs_{study}.done"
    output:
        CheckGenomad = "Checks/4.1-Genomad_all_{study}.done",
        Contigs1kbp = temp("genomad/{study}/{study}_1000bp_contigs.fna")
    conda:
        "genomad"
    threads: 12
    params:
        tag = "{study}",
        GenomadDB = "/group/jbemersogrp/databases/genomad/genomad_db_1.9",
        GenomadFolder = "genomad/{study}"
    resources:
        mem_mb = "100gb",
        partition = "low",
        time = "7-00:00:00"
    shell:'''
        # Filter contigs to >= 1000 bp
        seqtk seq -L 1000 {input.ContigsIn} > {output.Contigs1kbp} && \
        mkdir -p {params.GenomadFolder} && \
        mamba run -n genomad genomad end-to-end \
                    --cleanup \
                    --composition virome \
                    --enable-score-calibration \
                    -t {threads} {output.Contigs1kbp} \
        {params.GenomadFolder} {params.GenomadDB} && \
        touch {output.CheckGenomad}
        '''

rule FilterGenomad:
    input:
        CheckGenomad = "Checks/4.1-Genomad_all_{study}.done"
    output:
        FilteredGenomadResults = "genomad/{study}/{study}_FilteredGenomadResults.tsv",
        CheckFilteredGenomadResults = "Checks/4.3-FilterGenomadResults_{study}.done",
        FilteredContigs = "genomad/{study}_FilteredContigs_all.fna"
    conda:
        "genomad"
    resources:
        mem_mb = "100gb",
        partition = "low",
        time = "7-00:00:00"
    params:
        tag = "{study}",
        Contigs = "genomad/{study}/{study}_1000bp_contigs_summary/{study}_1000bp_contigs_virus.fna",
        GenomadResults = "genomad/{study}/{study}_1000bp_contigs_summary/{study}_1000bp_contigs_virus_summary.tsv"
    shell:'''
        awk -F'\t' 'NR==1 || $2 >= 10000 || ($2 >= 1000 && $11 ~ /Monodnaviria/)' {params.GenomadResults} | \
        cut -f1 > {output.FilteredGenomadResults} && \
        seqtk subseq {params.Contigs} {output.FilteredGenomadResults} > {output.FilteredContigs} && \
        touch {output.CheckFilteredGenomadResults}
        '''
#concatenate all contigs together for next step (clustering)
rule CatGenomadResults:
    input:
        Contigs = expand("genomad/{study}_FilteredContigs_all.fna", study = all_studies)
    output:
        CheckCatGenomadResults = "Checks/4.4-CatGenomadResults_all.done",
        CatGenomadContigs = "genomad/individualFilteredViralContigs.fna"
    params:
        tag = "individual"
    shell:'''
        cat {input.Contigs} > {output.CatGenomadContigs} && \
        touch {output.CheckCatGenomadResults}
        '''

#extract contig lengths, will need for vOTU clustering

rule contig_lengths:
    input:
        fasta="genomad/individualFilteredViralContigs.fna"
    output:
        tsv="genomad/contig_lengths.tsv",
        lengthcheck="Checks/lengthcheck.done"
    conda:
        "seqkit"
    shell:
        """
        seqkit fx2tab -nl {input.fasta} | awk '{{print $1"\t"$2}}' \
            > {output.tsv} && \
        touch {output.lengthcheck}
        """
################################################################################
#### 4.2 Run CheckV
################################################################################

rule checkv:
    """Run CheckV on GeNomad viral contigs to assess quality and completeness."""
    input:
        CheckGenomad = "Checks/4.1-Genomad_all_{study}.done"
    output:
        Checkcheckv = "Checks/{study}_checkv.done"
    params:
        genomad_out = "genomad/{study}",
        checkv_out = "checkv/{study}",
        checkv_db = "/group/jbemersogrp/databases/checkv/checkv-db-v1.5",
        tag = "{study}"
    threads: 16
    resources:
        mem_mb = 20000,
        time = "08:00:00",
        partition = "bmh"
    message: "🧫 Running CheckV on GeNomad output for {wildcards.sample}"
    shell: """
        mkdir -p {params.checkv_out} && \
        checkv end_to_end \
        {params.genomad_out}/{wildcards.study}.contigs_RENAMED_summary/{wildcards.study}.contigs_RENAMED_virus.fna \
            {params.checkv_out} \
            -d {params.checkv_db} \
            -t {threads} && \
        touch {output.Checkcheckv}
    """

################################################################################
#### 4.3 Add Sample ID to FASTA Headers & TSVs
################################################################################

#rule add_sample_headers:
#    """Append sample names to FASTA headers and TSV identifiers for merged analysis."""
#    input:
#        Checkcheckv = "Checks/4-viral_identification/{sample}_checkv.done"
#    output:
#        Check_checkv = "Checks/4-viral_identification/{sample}_headeradd.done"
#    params:
#        checkv_dir = "4-viral_identification/checkv/{sample}",
#        genomad_summary = "4-viral_identification/genomad/{sample}/{sample}.contigs_summary",
#        tag = "{sample}_add_sample_headers"
#    threads: 4
#    resources:
#        mem_mb = 10000,
#        time = "06:00:00",
#        partition = "med2"
#    message: "🧾 Adding sample name to FASTA and TSV headers for {wildcards.sample}"
#    shell: """
#        # Add sample prefix to virus FASTA headers
#        sed "s/^>\(.*\)/>{wildcards.sample}_\\1/" \
#            "{params.checkv_dir}/viruses.fna" > "{params.checkv_dir}/viruses_header.fna"

        # Add sample prefix to provirus FASTA headers
#        sed "s/^>\(.*\)/>{wildcards.sample}_\\1/" \
#            "{params.checkv_dir}/proviruses.fna" > "{params.checkv_dir}/proviruses_header.fna"

        # Add sample prefix to contamination table
 #       awk -v s={wildcards.sample}_ -F'\\t' '{{OFS="\\t"; $1=s $1}}1' \
 #           "{params.checkv_dir}/contamination.tsv" > "{params.checkv_dir}/contamination_header.tsv"

        # Add sample prefix to GeNomad virus summary table
 #       awk -v s={wildcards.sample}_ -F'\\t' '{{OFS="\\t"; $1=s $1}}1' \
 #           "{params.genomad_summary}/{wildcards.sample}.contigs_virus_summary.tsv" \
 #           > "{params.genomad_summary}/{wildcards.sample}.contigs_virus_summary_header.tsv"

  #      touch {output.Check_checkv}
  #  """

################################################################################
#### 4.4 Combine CheckV + GeNomad Outputs Across All Samples
################################################################################

# NOTE: This rule assumes you have a list of samples as a global variable, e.g. `samples = [...]`
#rule combine_checkv_results:
#    """
#    Merge all per-sample CheckV and GeNomad results into combined FASTA and summary tables.
#    Includes sanity checks for total sequences and lines.
#    """
#    input:
#        headeradd_done = expand("Checks/4-viral_identification/{sample}_headeradd.done", sample=samples)
#    output:
#        Check_combinecheckv = "Checks/4-viral_identification/genomad_checkv_merge.done",
#        combined_fasta = "4-viral_identification/checkv/combined_viruses.fna",
#        combined_summary = "4-viral_identification/checkv/all_genomad_files.tsv",
#        combined_checkv = "4-viral_identification/checkv/checkv_contamination.tsv"
#    params:
#        provirus_files = expand("4-viral_identification/checkv/{sample}/proviruses_header.fna", sample=samples),
#        virus_files = expand("4-viral_identification/checkv/{sample}/viruses_header.fna", sample=samples),
#        contamination_files = expand("4-viral_identification/checkv/{sample}/contamination_header.tsv", sample=samples),
#        summary_files = expand(
#            "4-viral_identification/genomad/{sample}/{sample}.contigs_summary/{sample}.contigs_virus_summary_header.tsv",
#            sample=samples
#        ),
#        tag = "combine_checkv_results"
#    threads: 2
#    resources:
#        mem_mb = 10000,
#        time = "06:00:00",
#        partition = "med2"
#    message: "📦 Combining CheckV + GeNomad results across all samples"
#    shell: r"""
#        set -euo pipefail

#        echo "Merging provirus + virus FASTA..."
#        cat {params.provirus_files} {params.virus_files} > {output.combined_fasta}
#        sed -i 's/ /_/g; s|/|_|g; s/|/__/g' {output.combined_fasta}
#        echo "Total viral sequences in combined FASTA:"
 #       grep -c ">" {output.combined_fasta}

 #       echo "Merging GeNomad summary tables..."
 #       cat {params.summary_files} > {output.combined_summary}
 #       echo "Total lines in combined summary:"
 #       wc -l < {output.combined_summary}

#        echo "Merging CheckV contamination tables..."
#        cat {params.contamination_files} > {output.combined_checkv}
#        echo "Total lines in combined contamination table:"
#        wc -l < {output.combined_checkv}

#        touch {output.Check_combinecheckv}
#    """
