"""AgroLD_ETL pipeline, in progress: one rule per data source.

Each source's parser lives in riceKB/ and is invoked as a script (not
imported), so this Snakefile only needs to know its command line, not its
internals. Paths come from config.yaml.

    snakemake --cores 1 --configfile config.yaml
    snakemake --cores 1 -n                       # dry run

Only interpro and plantTFDB are wired in so far (stringDB has no parser yet
to call; uniprotToRDF is handled separately). Add a source by giving it its
own `rule <name>:` below and listing it in `rule all`'s inputs.
"""
configfile: "config.yaml"

OUTPUT_DIR = config["output_dir"]


rule all:
    input:
        f"{OUTPUT_DIR}/interpro.ttl",
        f"{OUTPUT_DIR}/planttfdb.ttl",
        f"{OUTPUT_DIR}/planttfdb_regulation.ttl",


rule interpro:
    """InterPro domain hierarchy, labels and GO mapping."""
    input:
        tree_file=config["interpro"]["tree_file"],
        entry_list=config["interpro"]["entry_list"],
        go_mapping=config["interpro"]["go_mapping"],
    output:
        f"{OUTPUT_DIR}/interpro.ttl",
    shell:
        "python3 riceKB/interpro.py {input.tree_file} {input.entry_list} "
        "{input.go_mapping} {output}"


rule planttfdb:
    """PlantTFDB gene families and the regulatory relations between them.

    One parser invocation writes both outputs; declaring them together as
    Snakemake output means a resume treats them as produced atomically (one
    missing file reruns the whole rule, not a partial one).
    """
    input:
        tf_list=config["planttfdb"]["tf_list"],
        regulation_file=config["planttfdb"]["regulation_file"],
    output:
        tf=f"{OUTPUT_DIR}/planttfdb.ttl",
        regulation=f"{OUTPUT_DIR}/planttfdb_regulation.ttl",
    shell:
        "python3 riceKB/plantTFDB.py {input.tf_list} {output.tf} "
        "{input.regulation_file} {output.regulation}"
