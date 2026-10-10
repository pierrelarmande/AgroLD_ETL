"""AgroLD_ETL pipeline, in progress: one rule per data source.

Each source's parser lives in riceKB/ and is invoked as a script (not
imported), so this Snakefile only needs to know its command line, not its
internals. Paths come from config.yaml.

    snakemake --cores 1 --configfile config.yaml
    snakemake --cores 1 -n                       # dry run
    snakemake --profile profiles/slurm           # on the cluster, one sbatch job per rule

Only interpro, plantTFDB and stringDB are wired in so far (uniprotToRDF is
handled separately). Add a source by giving it its own `rule <name>:` below
and listing it in `rule all`'s inputs.

Both current rules are light enough that profiles/slurm/config.yaml's
default resources cover them; give a heavier future rule its own
`resources:` block (e.g. a bigger mem_mb_per_cpu or a different
slurm_partition) rather than raising that shared default.
"""
import os
from snakemake.exceptions import WorkflowError

configfile: "config.yaml"

OUTPUT_DIR = config["output_dir"]


def load_agrold_taxon_ids():
    """{production_name: taxon_id} for the species list ensembl_rdf and this
    pipeline share (config.yaml's ensembl_rdf section) -- the single source
    of truth for which species exist and what each one's NCBI taxon is, so
    the two pipelines never disagree. species_agrold.yaml names the species;
    species_EnsemblPlants_63.txt (tab-separated, #name/species/.../
    taxonomy_id/...) is where the taxon id actually comes from.
    """
    import yaml

    species_file = config["ensembl_rdf"]["species_file"]
    species_table = config["ensembl_rdf"]["species_table"]
    for path in (species_file, species_table):
        if not os.path.exists(path):
            raise WorkflowError(
                f"ensembl_rdf config not found at {path} -- ensembl_rdf is "
                "expected as a sibling checkout; edit config.yaml's "
                "ensembl_rdf section if yours lives elsewhere.")

    with open(species_file) as f:
        wanted = set(yaml.safe_load(f)["species"])

    taxon_id = {}
    with open(species_table) as f:
        next(f)  # header
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name = fields[1]
            if name in wanted:
                taxon_id[name] = fields[3]

    missing = wanted - set(taxon_id)
    if missing:
        raise WorkflowError(
            f"{species_table} has no taxon id for: {sorted(missing)}")
    return taxon_id


AGROLD_TAXON_ID = load_agrold_taxon_ids()


def planttfdb_species():
    """Species plantTFDB actually has input data for: one subdirectory per
    species under planttfdb.data_dir, named by its ensembl_rdf production
    name. A subdirectory for a species outside the AgroLD list is an error,
    not a silent skip -- it is more likely a typo than an intentional extra
    species the shared config doesn't know about yet.
    """
    data_dir = config["planttfdb"]["data_dir"]
    if not os.path.isdir(data_dir):
        return []
    found = sorted(d for d in os.listdir(data_dir)
                   if os.path.isdir(os.path.join(data_dir, d)))
    unknown = [d for d in found if d not in AGROLD_TAXON_ID]
    if unknown:
        raise WorkflowError(
            f"{data_dir} has data for species not in {config['ensembl_rdf']['species_file']}: "
            f"{unknown}")
    return found


PLANTTFDB_SPECIES = planttfdb_species()


def stringdb_species():
    """Species StringDB actually has input data for: one subdirectory per
    species under stringdb.data_dir, same convention and same unknown-
    species error as planttfdb_species() above.
    """
    data_dir = config["stringdb"]["data_dir"]
    if not os.path.isdir(data_dir):
        return []
    found = sorted(d for d in os.listdir(data_dir)
                   if os.path.isdir(os.path.join(data_dir, d)))
    unknown = [d for d in found if d not in AGROLD_TAXON_ID]
    if unknown:
        raise WorkflowError(
            f"{data_dir} has data for species not in {config['ensembl_rdf']['species_file']}: "
            f"{unknown}")
    return found


STRINGDB_SPECIES = stringdb_species()


rule all:
    input:
        f"{OUTPUT_DIR}/interpro.ttl",
        expand(f"{OUTPUT_DIR}/planttfdb_{{species}}.ttl", species=PLANTTFDB_SPECIES),
        expand(f"{OUTPUT_DIR}/planttfdb_{{species}}_regulation.ttl", species=PLANTTFDB_SPECIES),
        expand(f"{OUTPUT_DIR}/stringdb_{{species}}.ttl", species=STRINGDB_SPECIES),


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
    """PlantTFDB gene families and the regulatory relations between them, for
    one species. The taxon id comes from the shared ensembl_rdf species
    list (load_agrold_taxon_ids above), not from the input file, which
    carries no species information of its own.

    One parser invocation writes both outputs; declaring them together as
    Snakemake output means a resume treats them as produced atomically (one
    missing file reruns the whole rule, not a partial one).
    """
    input:
        tf_list=lambda wc: os.path.join(
            config["planttfdb"]["data_dir"], wc.species, config["planttfdb"]["tf_list"]),
        regulation_file=lambda wc: os.path.join(
            config["planttfdb"]["data_dir"], wc.species, config["planttfdb"]["regulation_file"]),
    output:
        tf=f"{OUTPUT_DIR}/planttfdb_{{species}}.ttl",
        regulation=f"{OUTPUT_DIR}/planttfdb_{{species}}_regulation.ttl",
    params:
        taxon_id=lambda wc: AGROLD_TAXON_ID[wc.species],
    shell:
        "python3 riceKB/plantTFDB.py {input.tf_list} {output.tf} "
        "{input.regulation_file} {output.regulation} --taxon-id {params.taxon_id}"


rule stringdb:
    """STRING protein-protein association network for one species, with
    STRING's own protein ids resolved to native gene locus ids via the
    aliases file (see riceKB/stringDB.py). Output is large -- for
    Arabidopsis thaliana alone the links file's ~15M rows produce a ~5GB
    Turtle file -- so this rule may need its own heavier `resources:` for
    bigger genomes, per the module docstring above.
    """
    input:
        aliases_file=lambda wc: os.path.join(
            config["stringdb"]["data_dir"], wc.species,
            f"{AGROLD_TAXON_ID[wc.species]}.{config['stringdb']['aliases_suffix']}"),
        links_file=lambda wc: os.path.join(
            config["stringdb"]["data_dir"], wc.species,
            f"{AGROLD_TAXON_ID[wc.species]}.{config['stringdb']['links_suffix']}"),
    output:
        f"{OUTPUT_DIR}/stringdb_{{species}}.ttl",
    params:
        taxon_id=lambda wc: AGROLD_TAXON_ID[wc.species],
    shell:
        "python3 riceKB/stringDB.py {input.aliases_file} {input.links_file} "
        "{output} --taxon-id {params.taxon_id}"
