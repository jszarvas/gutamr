# gutamr

Snakemake pipeline to predict and quantify acquired antimicrobial-resistance genes in sequencing reads or assembled data, using [CARD's]("https://card.mcmaster.ca/") `nucleotide_fasta_protein_homolog_model.fasta`.  

Integrated also into [MIntO](https://github.com/arumugamlab/MIntO). 

The pipeline needs to be configured by filling out the `gutamr.yaml` file. It expects the assemblies in one directory, and the reads in one sub-directory per sample under the specified folder. See [the config file for the MIntO tutorial](https://github.com/arumugamlab/MIntO/configuration/module_arg.yaml) for an example.  
It outputs tab-separated files with the filtered KMA results for both types of inputs, and an FPKM table for the reads.  
**NOTE:** Outputs are not filtered to contain only clinically relevant AROs. That is left to the user's discretion!

Example command:
```
snakemake --use-conda --keep-going --latency-wait 60 --cores 20 --resources mem=100 --conda-prefix /path/gutamr/conda_env --shadow-prefix /scratch/$USER --default-resources mem=4 --snakefile /path/gutamr/gutamr.smk --configfile gutamr.yaml -np --verbose | less
```