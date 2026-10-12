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

plantTFDB's and stringdb's raw input is fetched on demand (download_*
rules) and deleted again as soon as it's converted (temp() outputs) --
unlike interpro, which still expects data already placed under data_dir by
hand. Run with `snakemake --keep-going` so a species one of these sources
has no data for doesn't block the others.

Most rules are light enough that profiles/slurm/config.yaml's default
resources cover them; give a heavier future rule its own `resources:` block
(e.g. a bigger mem_mb_per_cpu or a different slurm_partition) rather than
raising that shared default.
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


# Like stringdb below, plantTFDB/PlantRegMap's raw input is downloaded on
# demand by the download_planttfdb rule, so every AgroLD species is a
# candidate; a species neither site covers simply fails its own
# download_planttfdb job (see bin/download_planttfdb.py) -- run with
# `snakemake --keep-going` so one such failure doesn't block the others.
PLANTTFDB_SPECIES = sorted(AGROLD_TAXON_ID)


# Unlike planttfdb (no fetch script wired in yet), stringdb's raw input is
# downloaded on demand by the download_stringdb rule below, so there is no
# pre-existing data_dir to scan: every AgroLD species is a candidate. A
# species STRING has no data for simply fails its own download_stringdb job
# (see bin/download_stringdb.py) -- run with `snakemake --keep-going` so one
# such failure doesn't block the others.
STRINGDB_SPECIES = sorted(AGROLD_TAXON_ID)


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


rule download_planttfdb:
    """Fetch one species' PlantTFDB TF list + PlantRegMap merged regulation
    network on demand, via bin/download_planttfdb.py. Declared temp() so
    Snakemake deletes them as soon as the planttfdb rule below has consumed
    them -- same reasoning as download_stringdb below: no permanent local
    copy once converted.

    Unlike download_stringdb, output filenames here don't embed a taxon id
    (PlantTFDB/PlantRegMap files are always named the same within a species'
    subdirectory), so {species} is the only wildcard needed and no
    taxon-matching trick is required.
    """
    output:
        tf_list=temp(os.path.join(
            config["planttfdb"]["data_dir"], "{species}", config["planttfdb"]["tf_list"])),
        regulation_file=temp(os.path.join(
            config["planttfdb"]["data_dir"], "{species}", config["planttfdb"]["regulation_file"])),
    shell:
        "python3 bin/download_planttfdb.py -s {wildcards.species} "
        "-o {config[planttfdb][data_dir]} "
        "--tf-list-filename {config[planttfdb][tf_list]} "
        "--regulation-filename {config[planttfdb][regulation_file]} "
        "-f {config[ensembl_rdf][species_file]}"


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


rule download_stringdb:
    """Fetch one species' STRING aliases + links files on demand, via
    bin/download_stringdb.py. Declared temp() so Snakemake deletes them as
    soon as the stringdb rule below has consumed them: STRING's raw,
    unfiltered links file is several GB per species on disk, far more than
    the filtered RDF it gets converted to, so there is no reason to keep a
    permanent local copy per species once converted.

    {taxon} is a free wildcard here (not looked up from AGROLD_TAXON_ID
    directly, since output paths can't be computed by a function in
    Snakemake) -- it gets bound by matching this rule's output pattern
    against the concrete path the stringdb rule's input asks for, which
    already embeds the right taxon id.
    """
    output:
        aliases_file=temp(os.path.join(
            config["stringdb"]["data_dir"], "{species}",
            "{taxon}." + config["stringdb"]["aliases_suffix"])),
        links_file=temp(os.path.join(
            config["stringdb"]["data_dir"], "{species}",
            "{taxon}." + config["stringdb"]["links_suffix"])),
    wildcard_constraints:
        taxon=r"\d+",
    shell:
        "python3 bin/download_stringdb.py -s {wildcards.species} "
        "-o {config[stringdb][data_dir]} "
        "-f {config[ensembl_rdf][species_file]} "
        "-t {config[ensembl_rdf][species_table]}"


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
