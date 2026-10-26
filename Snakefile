all_samples = ['shCTL_cyto-1_R1', 'shCTL_cyto-1_R2', 'shCTL_cyto-2_R1', 'shCTL_cyto-2_R2', 'shCTL_nucleus-1_R1', 'shCTL_nucleus-1_R2', 'shCTL_nucleus-2_R1', 'shCTL_nucleus-2_R2',
               'shPABPN1_cyto-1_R1', 'shPABPN1_cyto-1_R2', 'shPABPN1_cyto-2_R1', 'shPABPN1_cyto-2_R2', 'shPABPN1_nucleus-1_R1', 'shPABPN1_nucleus-1_R2', 'shPABPN1_nucleus-2_R1', 'shPABPN1_nucleus-2_R2']
paired_samples = ['shCTL_cyto-1', 'shCTL_cyto-2', 'shCTL_nucleus-1', 'shCTL_nucleus-2', 'shPABPN1_cyto-1', 'shPABPN1_cyto-2', 'shPABPN1_nucleus-1', 'shPABPN1_nucleus-2']
merged_reps = ['shCTL_nucleus', 'shCTL_cyto', 'shPABPN1_nucleus', 'shPABPN1_cyto']
irfinder_shCTL_nucleus_reps = ['shCTL_nucleus', 'shCTL_nucleus-1', 'shCTL_nucleus-2']
irfinder_shCTL_cyto_reps = ['shCTL_cyto', 'shCTL_cyto-1', 'shCTL_cyto-2']
irfinder_shPABPN1_nucleus_reps = ['shPABPN1_nucleus', 'shPABPN1_nucleus-1', 'shPABPN1_nucleus-2']
irfinder_shPABPN1_cyto_reps = ['shPABPN1_cyto', 'shPABPN1_cyto-1', 'shPABPN1_cyto-2']

rule finale :
    input :
        fastqc_before_trim = expand("results/fastqc_before_trim/{sample}_fastqc.html", sample = all_samples),
        fastqc_after_trim_R1 = expand("results/fastqc_after_trim/{paired_sample}_R1_val_1_fastqc.html", paired_sample = paired_samples),
        fastqc_after_trim_R2= expand("results/fastqc_after_trim/{paired_sample}_R2_val_2_fastqc.html", paired_sample = paired_samples),
        read_counts = "results/read_counts/read_counts.out",
        irfinder_diff_nuc = "results/irfinder/irfinder_diff_nucleus/shPABPN1_vs_shCTL_nucleus.txt",
        irfinder_diff_cyto = "results/irfinder/irfinder_diff_nucleus/shPABPN1_vs_shCTL_cyto.txt"

rule fastqc_before_trim :
    threads : 1
    input:
        "/home/sbache/projects/def-kchoquet/sequencing_data_repo/NGS/250508_SB1_shCTL_shPABPN1_nucleus_cyto/{sample}.fastq.gz"
    output:
        "results/fastqc_before_trim/{sample}_fastqc.html"
    params :
        outdir = "results/fastqc_before_trim"
    log:
        error='logs_error/fastqc_before_trim{sample}.e'
    shell :
        """
        module load fastqc/0.12.1
        fastqc -o {params.outdir} {input} 2>{log.error}
        """

rule trimgalore :
    threads : 1
    input:
        input_R1 = "/home/sbache/projects/def-kchoquet/sequencing_data_repo/NGS/250508_SB1_shCTL_shPABPN1_nucleus_cyto/{paired_sample}_R1.fastq.gz",
        input_R2= "/home/sbache/projects/def-kchoquet/sequencing_data_repo/NGS/250508_SB1_shCTL_shPABPN1_nucleus_cyto/{paired_sample}_R2.fastq.gz"
    output :
        output_R1 = "results/trimgalore/{paired_sample}_R1_val_1.fq.gz",
        output_R2= "results/trimgalore/{paired_sample}_R2_val_2.fq.gz"
    params :
        outdir = "results/trimgalore"
    log:
        error='logs_error/trimgalore_{paired_sample}.e'
    conda :
        "envs/trimgalore.yaml"
    shell :
        """
        trim_galore --paired -o {params.outdir} {input.input_R1} {input.input_R2} 2>{log.error}
        """

rule fastqc_after_trim :
    threads : 1
    input:
        input_R1 = "results/trimgalore/{paired_sample}_R1_val_1.fq.gz",
        input_R2 = "results/trimgalore/{paired_sample}_R2_val_2.fq.gz"
    output:
        output_R1 = "results/fastqc_after_trim/{paired_sample}_R1_val_1_fastqc.html",
        output_R2 = "results/fastqc_after_trim/{paired_sample}_R2_val_2_fastqc.html"
    params :
        outdir = "results/fastqc_after_trim"
    log:
        error='logs_error/fastqc_after_trim_{paired_sample}.e'
    shell :
        """
        module load fastqc/0.12.1
        fastqc -o {params.outdir} {input.input_R1} 2>{log.error}
        fastqc -o {params.outdir} {input.input_R2}
        """

rule star_build_index :
    threads : 8
    input :
        fasta_ref = "references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        gtf = "references/Homo_sapiens.GRCh38.114.gtf"
    output :
        "results/star_ref/Genome"
    params :
        outdir = "results/star_ref"
    log:
        error='logs_error/star_build_index.e'
    shell :
        """
        module load star/2.7.11b
        STAR --runMode genomeGenerate --genomeDir {params.outdir} --genomeFastaFiles {input.fasta_ref} --sjdbGTFfile {input.gtf} --runThreadN {threads} 2>{log.error}
        """

rule star_alignment :
    threads : 8
    input :
        star_genome = "results/star_ref/Genome",
        R1 = "results/trimgalore/{paired_sample}_R1_val_1.fq.gz",
        R2= "results/trimgalore/{paired_sample}_R2_val_2.fq.gz"
    output :
        "results/alignment/{paired_sample}Aligned.sortedByCoord.out.bam"
    params :
        gunzip_R1 = "results/trimgalore/{paired_sample}_R1_val_1.fq",
        gunzip_R2= "results/trimgalore/{paired_sample}_R2_val_2.fq",
        star_ref_dir = "results/star_ref/",
        prefix = "results/alignment/{paired_sample}"
    log:
        error='logs_error/star_alignment_{paired_sample}.e'
    shell :
        """
        gunzip {input.R1} -c > {params.gunzip_R1}
        gunzip {input.R2} -c > {params.gunzip_R2}
        module load star/2.7.11b
        STAR --runMode alignReads --genomeDir {params.star_ref_dir} --outSAMtype BAM SortedByCoordinate --readFilesIn {params.gunzip_R1} {params.gunzip_R2} --runThreadN {threads} --outFileNamePrefix {params.prefix} 2>{log.error}
        """

rule featurecounts :
    threads: 4
    input :
        bamfiles = expand("results/alignment/{paired_sample}Aligned.sortedByCoord.out.bam", paired_sample = paired_samples)
    output :
        "results/read_counts/read_counts.out"
    params:
        gtf = "references/Homo_sapiens.GRCh38.114.gtf"
    log :
        error = 'logs_error/featurecounts.e'
    conda :
        "envs/featurecounts.yaml"
    shell :
        """
        featureCounts -a {params.gtf} -o {output} {input.bamfiles} -T {threads} -B -C -p -s 2 2>{log.error}
        """

rule irfinder_build_ref :
    threads: 5
    input:
        star_genome = "results/star_ref/Genome",
        fasta_ref= "references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        gtf="references/Homo_sapiens.GRCh38.114.gtf"
    output :
        "results/irfinder_ref/IRFinder/introns.unique.bed"
    params :
        outdir = "results/irfinder_ref",
        star_ref = "results/star_ref"
    log:
        error='logs_error/irfinder_build_ref.e'
    shell :
        """
        module load apptainer/1.4.5
        rm -rf {params.outdir}
        singularity run /home/sbache/projects/def-kchoquet/sbache/scripts/IRFinder BuildRefFromSTARRef -r {params.outdir} -x {params.star_ref} -f {input.fasta_ref} -g {input.gtf} -t {threads} 2>{log.error}
        """

rule irfinder_individual :
    threads : 3
    input :
        irfinder_introns = "results/irfinder_ref/IRFinder/introns.unique.bed",
        bam = "results/alignment/{paired_sample}Aligned.sortedByCoord.out.bam",
    output:
        "results/irfinder/{paired_sample}/IRFinder-IR-dir.txt"
    params:
        outdir = "results/irfinder/{paired_sample}",
        irfinder_ref = "results/irfinder_ref"
    log:
        error='logs_error/irfinder_individual_{paired_sample}.e'
    shell :
        """
        module load apptainer/1.4.5
        singularity run /home/sbache/projects/def-kchoquet/sbache/scripts/IRFinder BAM -r {params.irfinder_ref} -d {params.outdir} -t {threads} {input.bam} 2>{log.error}
        """

rule merge_bam_files :
    threads : 1
    input :
        rep1 = "results/alignment/{merged_rep}-1Aligned.sortedByCoord.out.bam",
        rep2= "results/alignment/{merged_rep}-2Aligned.sortedByCoord.out.bam"
    output :
        "results/alignment/merged_bam/{merged_rep}Aligned.sortedByCoord.out.bam"
    log:
        error='logs_error/merge_bam_files_{merged_rep}.e'
    shell :
        """
        module load StdEnv/2020
        module load samtools/1.17
        samtools cat {input.rep1} {input.rep2} > {output} 2>{log.error}
        """

rule irfinder_merged :
    threads : 3
    input :
        irfinder_introns = "results/irfinder_ref/IRFinder/introns.unique.bed",
        bam = "results/alignment/merged_bam/{merged_rep}Aligned.sortedByCoord.out.bam",
    output:
        "results/irfinder/{merged_rep}/IRFinder-IR-dir.txt"
    params:
        outdir = "results/irfinder/{merged_rep}",
        irfinder_ref = "results/irfinder_ref"
    log:
        error='logs_error/irfinder_merged_{merged_rep}.e'
    shell :
        """
        module load apptainer/1.4.5
        singularity run /home/sbache/projects/def-kchoquet/sbache/scripts/IRFinder BAM -r {params.irfinder_ref} -d {params.outdir} -t {threads} {input.bam} 2>{log.error}
        """

rule irfinder_diff_nucleus :
    threads : 3
    input :
        shCTL = expand("results/irfinder/{irfinder_shCTL_nucleus_rep}/IRFinder-IR-dir.txt", irfinder_shCTL_nucleus_rep = irfinder_shCTL_nucleus_reps),
        shPABPN1 = expand("results/irfinder/{irfinder_shPABPN1_nucleus_rep}/IRFinder-IR-dir.txt", irfinder_shPABPN1_nucleus_rep = irfinder_shPABPN1_nucleus_reps)
    output:
        "results/irfinder/irfinder_diff_nucleus/shPABPN1_vs_shCTL_nucleus.txt"
    log :
        error = 'logs_error/irfinder_diff_nucleus.e'
    shell :
        """
        module load perl
        /home/sbache/projects/def-kchoquet/sbache/scripts/IRFinder_dir/IRFinder-2.0-beta/bin/analysisWithLowReplicates.pl -A {input.shCTL} -B {input.shPABPN1} > {output} 2>{log.error}
        """

rule irfinder_diff_cytoplasm :
    threads : 3
    input :
        shCTL = expand("results/irfinder/{irfinder_shCTL_cyto_rep}/IRFinder-IR-dir.txt", irfinder_shCTL_cyto_rep = irfinder_shCTL_cyto_reps),
        shPABPN1 = expand("results/irfinder/{irfinder_shPABPN1_cyto_rep}/IRFinder-IR-dir.txt", irfinder_shPABPN1_cyto_rep = irfinder_shPABPN1_cyto_reps)
    output:
        "results/irfinder/irfinder_diff_nucleus/shPABPN1_vs_shCTL_cyto.txt"
    log :
        error = 'logs_error/irfinder_diff_cytoplasm.e'
    shell :
        """
        module load perl
        /home/sbache/projects/def-kchoquet/sbache/scripts/IRFinder_dir/IRFinder-2.0-beta/bin/analysisWithLowReplicates.pl -A {input.shCTL} -B {input.shPABPN1} > {output} 2>{log.error}
        """

