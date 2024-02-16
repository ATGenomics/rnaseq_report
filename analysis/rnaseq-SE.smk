from pathlib import Path

wd = str(Path.home() / "workshop/rnaseq_report/analysis")

sample_links = {
    "ERR458493": "https://osf.io/5daup/download",
    "ERR458494": "https://osf.io/8rvh5/download",
    "ERR458495": "https://osf.io/2wvn3/download",
    "ERR458500": "https://osf.io/xju4a/download",
    "ERR458501": "https://osf.io/nmqe6/download",
    "ERR458502": "https://osf.io/qfsze/download",
}

SAMPLES = sample_links.keys()

rule all:
    input:
        expand(["rnaseqSE_snakemake/02_mapping/salmon/{name}_quant/quant.sf",
                "rnaseqSE_snakemake/02_mapping/kallisto/{name}_quant/abundance.h5"], name=SAMPLES),
        "rnaseqSE_snakemake/01_qc/multiqc_report.html"


rule download_reads:
    output:
        "rnaseqSE_snakemake/00_raw/{sample}.fq.gz"
    params:
        download_link=lambda wildcards: sample_links[wildcards.sample]
    shell:
        """
        curl -L {params.download_link} -o {output}
        """


rule fastqc_raw:
    input:
        "rnaseqSE_snakemake/00_raw/{sample}.fq.gz"
    output:
        "rnaseqSE_snakemake/01_qc/a_raw_fastqc/{sample}_fastqc.html"
    params:
        outdir="rnaseqSE_snakemake/01_qc/a_raw_fastqc"
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        fastqc {input} -f fastq --nogroup --outdir {params.outdir}
        """


rule adapter_file:
    output:
        temp("TruSeq2-SE.fa")
    params:
        url="https://raw.githubusercontent.com/ATGenomics/adapters/master"
    shell:
        """ 
        curl -L {params.url}/{output} -o {output}
        """


rule quality_trim:
    input:
        reads="rnaseqSE_snakemake/00_raw/{sample}.fq.gz",
        adapters="TruSeq2-SE.fa",
    output:
        "rnaseqSE_snakemake/01_qc/b_trimming/{sample}.qc.fq.gz",
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        trimmomatic SE {input.reads} {output} \
                ILLUMINACLIP:{input.adapters}:2:0:15 \
                LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25    
        """


rule fastqc_trimmed:
    input:
        "rnaseqSE_snakemake/01_qc/b_trimming/{sample}.qc.fq.gz"
    output:
        "rnaseqSE_snakemake/01_qc/c_trim_fastqc/{sample}.qc_fastqc.html"
    params:
        outdir="rnaseqSE_snakemake/01_qc/c_trim_fastqc"
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        fastqc {input} -f fastq --outdir {params.outdir}
        """


rule multiqc:
    input:
        raw=expand("rnaseqSE_snakemake/01_qc/a_raw_fastqc/{sample}_fastqc.html", sample=SAMPLES),
        trimmed=expand("rnaseqSE_snakemake/01_qc/c_trim_fastqc/{sample}.qc_fastqc.html", sample=SAMPLES),
    output:
        "rnaseqSE_snakemake/01_qc/multiqc_report.html"
    params:
        raw_dir="rnaseqSE_snakemake/01_qc/a_raw_fastqc",
        trimmed_dir="rnaseqSE_snakemake/01_qc/c_trim_fastqc"
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        multiqc {params.raw_dir} {params.trimmed_dir} --filename {output}
        """


rule download_yeast_ref:
    output:
        "rnaseqSE_snakemake/02_mapping/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
    params:
        url="ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/cdna"
    shell:
        """
        curl -L {params.url}/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -o {output}
        """


rule salmon_index:
    input:
        "rnaseqSE_snakemake/02_mapping/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
    output:
        directory("rnaseqSE_snakemake/02_mapping/reference/salmon/salmon_index")
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        salmon index --index {output} --transcripts {input} # --type quasi
        """


rule salmon_quantify:
    input:
        reads="rnaseqSE_snakemake/01_qc/b_trimming/{sample}.qc.fq.gz",
        index_dir="rnaseqSE_snakemake/02_mapping/reference/salmon/salmon_index",
    output:
        "rnaseqSE_snakemake/02_mapping/salmon/{sample}_quant/quant.sf"
    params:
        outdir=lambda wildcards: "rnaseqSE_snakemake/02_mapping/salmon/"
        + wildcards.sample
        + "_quant"
    threads: 4
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} \
                --threads {threads} --seqBias --gcBias --validateMappings
        """


rule kallisto_index:
    input:
        "rnaseqSE_snakemake/02_mapping/reference/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
    output:
        "rnaseqSE_snakemake/02_mapping/reference/kallisto/kallisto_index"
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        kallisto index --index {output} {input}
        """


rule kallisto_quantify:
    input:
        index_file="rnaseqSE_snakemake/02_mapping/reference/kallisto/kallisto_index",
        reads="rnaseqSE_snakemake/01_qc/b_trimming/{sample}.qc.fq.gz",
    output:
        "rnaseqSE_snakemake/02_mapping/kallisto/{sample}_quant/abundance.h5"
    params:
        outdir=lambda wildcards: "rnaseqSE_snakemake/02_mapping/kallisto/"
        + wildcards.sample
        + "_quant",
    threads: 4
    conda:
        f"{wd}/env_files/rnaseqSE.yaml"
    shell:
        """
        kallisto quant --single -l 200 -s 30 -i {input.index_file} -o {params.outdir} \
            -t {threads} {input.reads}
        """
