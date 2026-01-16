import os

wd = config["working_dir"]
database_dir = config["database_dir"]
card_version = config["CARD_version"]

fasta_suffix = config["assembly_suffix"]
fasta_dir = config["assemblies_folder"]
fastq_suffix = config["read_suffix"][0]
fastq_dir = config["read_folder"]

def output_combined_amr_table():
    result = []
    if config["run_assemblies"]:
        result.extend(expand("{wd}/output/assembly_card_{version}_cov{coverage}.tsv",
        wd = wd,
        version = card_version,
        coverage = config["min_template_coverage"]))
    if config["run_reads"]:
        result.extend(expand("{wd}/output/reads_card_{version}_cov{coverage}.fpkm.tsv",
        wd = wd,
        version = card_version,
        coverage = config["min_template_coverage"]))
    return(result)

rule all:
    input:
        output_combined_amr_table()

####################################
# Database download and indexing
####################################

rule card_download:
    output:
        "{database_dir}/card_db_{version}/nucleotide_fasta_protein_homolog_model.fasta"
    params:
        prefix=os.path.join(database_dir, f"card_db_{card_version}")
    localrule:
        True
    threads:
        1
    shell:
        """
        wget -P {params.prefix} https://card.mcmaster.ca/download/0/broadstreet-v{card_version}.tar.bz2
        cd {params.prefix}
        tar -xjf broadstreet-v{card_version}.tar.bz2
        """

rule kma_index_phm:
    input:
        "{database_dir}/card_db_{version}/nucleotide_fasta_protein_homolog_model.fasta"
    output:
        "{database_dir}/card_db_{version}/nucleotide_fasta_protein_homolog_model.fasta.name"
    log:
        "{database_dir}/card_db_{version}/kma_index_{version}.log"
    params:
        prefix=os.path.join(database_dir, f"card_db_{card_version}")
    localrule:
        True
    conda:
        "envs/aln.yml"
    threads:
        3
    resources:
        mem_mb = 5000
    shell:
        """
        time (
            kma index -i {input}
        ) >& {log}
        """

################
# Assemblies
################

def get_assembly_names():
    from glob import glob
    assembly_files = glob(f"{fasta_dir}/*.{fasta_suffix}")
    assembly_names = [os.path.basename(x).replace(f".{fasta_suffix}", "") for x in assembly_files]
    return assembly_names

rule run_kma_fasta:
    input:
        index_file = expand("{database_dir}/card_db_{version}/nucleotide_fasta_protein_homolog_model.fasta.name",
          database_dir = database_dir,
          version = card_version),
        fasta=f"{fasta_dir}/{{assembly_name}}.{fasta_suffix}"
    output:
        "{wd}/kma_fasta/{assembly_name}.res"
    log:
        "{wd}/logs/kma_{assembly_name}.log"
    params:
        output_suffix= lambda wildcards: f"{wildcards.wd}/kma_fasta/{wildcards.assembly_name}"
    localrule:
        True
    conda:
        "envs/aln.yml"
    threads:
        3
    resources:
        mem_mb = 5000
    shell:
        """
        time (
            kma_db_prefix=$(dirname {input.index_file})
            kma -t_db $kma_db_prefix/nucleotide_fasta_protein_homolog_model.fasta -t {threads} -cge -na -nc -nf -o {params.output_suffix} -i {input.fasta}
        ) >& {log}
        """

rule process_kma_fasta:
    input:
        res = expand("{wd}/kma_fasta/{assembly_name}.res",
                    wd = wd,
                    assembly_name = get_assembly_names())
    output:
        "{workdir}/output/assembly_card_{version}_cov{coverage}.tsv"
    log:
        "{workdir}/logs/assembly_card_{version}_cov{coverage}.log"
    params:
        input_prefix= lambda wildcards: f"{wildcards.workdir}/kma_fasta",
        max_id_diff = config["max_template_identity_difference"]
    localrule:
        True
    conda:
        "envs/aln.yml"
    threads:
        3
    resources:
        mem_mb = 5000
    shell:
        """
        time (
            Rscript {workflow.basedir}/scripts/combine_normalized_amr.R --folder {params.input_prefix} --output {output} --type assembly --min_cov {wildcards.coverage} --max_id_diff {params.max_id_diff}
        ) >& {log}
        """

################
# Trimmed reads
################

def get_sample_names():
    samples = []
    for obj in os.scandir(fastq_dir):
        if not obj.name.startswith('.') and obj.is_dir():
            samples.append(obj.name)
    return(samples)
    
def get_sample_runs(wildcards):
    runs = []
    sample_dir = os.path.join(fastq_dir, wildcards.sample)
    for obj in os.scandir(sample_dir):
        if obj.name.endswith(fastq_suffix) and obj.is_file():
            runs.append(obj.path)
    return(runs)

rule run_kma_fastq:
    input:
        index_file = expand("{database_dir}/card_db_{version}/nucleotide_fasta_protein_homolog_model.fasta.name",
          database_dir = database_dir,
          version = card_version),
        fwd=get_sample_runs
    output:
        res = "{wd}/kma_fastq/{sample}.res",
        mst = "{wd}/kma_fastq/{sample}.mapstat"
    log:
        "{wd}/logs/kma_{sample}.log"
    params:
        piped = lambda wildcards, input: "yes" if len(input.fwd) > 1 else "no",
        rev = lambda wildcards, input: [x.replace(fastq_suffix, fastq_suffix.replace("1", "2")) for x in input.fwd]
    shadow:
        "minimal"
    conda:
        "envs/aln.yml"
    threads:
        3
    resources:
        mem_mb = 15000
    shell:
        """
        time (
            if [ "{params.piped}" == "yes" ]; then
                mkfifo {wildcards.sample}.1.fq.gz
                mkfifo {wildcards.sample}.2.fq.gz
                cat {input.fwd} > {wildcards.sample}.1.fq.gz &
                cat {params.rev} > {wildcards.sample}.2.fq.gz &
                input_files="{wildcards.sample}.1.fq.gz {wildcards.sample}.2.fq.gz"
            else
                input_files="{input.fwd} {params.rev}"
            fi
            
            mkdir temp
            export TMPDIR="$PWD/temp"
            kma_db_prefix=$(dirname {input.index_file})
            kma -t_db $kma_db_prefix/nucleotide_fasta_protein_homolog_model.fasta -t {threads} -tmp $TMPDIR/ -1t1 -ID 50.0 -mrs 0.85 -ml 75 -cge -na -nc -nf -ef -o {wildcards.sample} -ipe $input_files
            rsync -a {wildcards.sample}.res {output.res}
            rsync -a {wildcards.sample}.mapstat {output.mst}
        ) >& {log}
        """

rule process_kma_fastq:
    input:
        res = expand("{wd}/kma_fastq/{sample}.res",
                    wd = wd,
                    sample = get_sample_names()),
        mst = expand("{wd}/kma_fastq/{sample}.mapstat",
                    wd = wd,
                    sample = get_sample_names())
    output:
        "{wd}/output/reads_card_{version}_cov{coverage}.fpkm.tsv"
    log:
        "{wd}/logs/reads_card_{version}_cov{coverage}.log"
    params:
        input_prefix = lambda wildcards: f"{wildcards.wd}/kma_fastq",
        min_count = config["min_fragment_no"],
        max_id_diff = config["max_template_identity_difference"]
    localrule:
        True
    conda:
        "envs/aln.yml"
    threads:
        3
    resources:
        mem_mb = 15000
    shell:
        """
        time (
            Rscript {workflow.basedir}/scripts/combine_normalized_amr.R --folder {params.input_prefix} --output {output} --type reads --min_cov {wildcards.coverage} --min_frag_count {params.min_count} --max_id_diff {params.max_id_diff}
        ) >& {log}
        """