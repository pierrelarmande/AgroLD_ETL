"""AgroLD_ETL pipeline, in progress: one rule per data source.

Each source's parser lives in riceKB/ and is invoked as a script (not
imported), so this Snakefile only needs to know its command line, not its
internals. Paths come from config.yaml.

    snakemake --cores 1 --configfile config.yaml
    snakemake --cores 1 -n                       # dry run
    snakemake --profile profiles/slurm           # on the cluster, one sbatch job per rule

interpro, oryzaBaseParser, plantTFDB, stringDB and uniprotToRDF are wired in
so far. Add a source by giving it its own `rule <name>:` below and listing
it in `rule all`'s inputs.

plantTFDB's, stringdb's and uniprot's raw input is fetched on demand
(download_* rules) and deleted again as soon as it's converted (temp()
outputs) -- unlike interpro, which still expects data already placed under
data_dir by hand. Run with `snakemake --keep-going` so a species one of
these sources has no data for doesn't block the others.

Most rules are light enough that profiles/slurm/config.yaml's default
resources cover them; give a heavier future rule its own `resources:` block
(e.g. a bigger mem_mb_per_cpu or a different slurm_partition) rather than
raising that shared default.

Every source's Turtle output is validated with rapper before rule all ever
sees it, via the generic validate_turtle rule below -- a rule that produces
Turtle names its own output *.raw.ttl (wrapped in temp(), since
validate_turtle is its only consumer) instead of *.ttl, and validate_turtle
turns it into the final *.ttl by the naming convention alone, no per-rule
wiring needed. Catches malformed syntax from a parser bug (escaping,
encoding, ...) in one pass over the whole finished file, the same mechanism
ensembl_rdf's bin/rapper_postprocess.sh uses.
"""
import os
from snakemake.exceptions import WorkflowError

configfile: "config.yaml"

# riceKB/*.py scripts do `from riceKB.globalVars import *`, which needs the
# repo root (not riceKB/ itself, the directory `python3 riceKB/x.py` puts on
# sys.path[0] by default) on PYTHONPATH. Don't rely on the submitting
# shell's own PYTHONPATH already having it (it won't, e.g. on a cluster
# compute node) -- set it explicitly for every rule's shell command instead.
shell.prefix(f'export PYTHONPATH="{workflow.basedir}:$PYTHONPATH"; ')

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


# UniProt publishes ONE flat-file per taxonomic division (e.g. "plants"),
# covering every species in it -- unlike planttfdb/stringdb above, there is
# no per-species variant to fetch, so this is keyed by taxon id rather than
# species name: riceKB/uniprotToRDF.py itself filters the species it wants
# out of the shared file (see download_uniprot/uniprot rules below).
# uniprot_trembl_plants.dat.gz is ~14GB -- both banks stay gzip-compressed
# on disk (download_uniprot doesn't decompress them, and uniprotToRDF.py
# streams/decompresses on the fly), so there's no decompressed copy eating
# scratch space, at the cost of re-decompressing the file in every one of
# UNIPROT_TAXA's jobs rather than once.
UNIPROT_BANKS = ["sprot", "trembl"]
UNIPROT_TAXA = sorted(set(AGROLD_TAXON_ID.values()))


rule all:
    input:
        f"{OUTPUT_DIR}/interpro.ttl",
        f"{OUTPUT_DIR}/oryzabase_genes.ttl",
        expand(f"{OUTPUT_DIR}/planttfdb_{{species}}.ttl", species=PLANTTFDB_SPECIES),
        expand(f"{OUTPUT_DIR}/planttfdb_{{species}}_regulation.ttl", species=PLANTTFDB_SPECIES),
        expand(f"{OUTPUT_DIR}/stringdb_{{species}}.ttl", species=STRINGDB_SPECIES),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.plants.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.keyword.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.pubmed.ttl", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),
        expand(f"{OUTPUT_DIR}/mappings_{{taxon}}_{{bank}}_gene_prot.tsv", taxon=UNIPROT_TAXA, bank=UNIPROT_BANKS),


rule validate_turtle:
    """Validate and normalize a pipeline's raw Turtle output with rapper --
    the same mechanism ensembl_rdf's bin/rapper_postprocess.sh uses, applied
    here as one generic rule instead of a per-source step. {name} is a free
    wildcard (can contain dots, e.g. "uniprot_4530_trembl.plants") matched
    against whatever *.raw.ttl a producing rule's output already names.

    rapper's own exit code is 1 on a genuine syntax error (verified against
    a deliberately malformed file), so Snakemake's normal "remove a failed
    job's declared outputs" behaviour already keeps a bad {output} from
    lingering. {output}.rapper.log is kept only when something failed, for
    debugging -- on success its content is already visible in Snakemake's
    own job log via the cat below, so it's removed rather than left as
    permanent clutter next to every one of the pipeline's Turtle files.
    """
    input:
        f"{OUTPUT_DIR}/{{name}}.raw.ttl",
    output:
        f"{OUTPUT_DIR}/{{name}}.ttl",
    wildcard_constraints:
        # Without this, {name} is free to end in ".raw" itself, so Snakemake
        # considers applying this very rule to its own output recursively
        # (name="X.raw" producing "X.raw.ttl" as an input to name="X.raw.raw"
        # producing "X.raw.raw.ttl", ...) -- confirmed: Snakemake rejects the
        # DAG with a PeriodicWildcardError/AmbiguousRuleException without
        # this constraint.
        name=r".+(?<!\.raw)",
    shell:
        "rapper -i turtle -o turtle {input} > {output} 2> {output}.rapper.log; "
        "status=$?; cat {output}.rapper.log >&2; "
        "if [ $status -eq 0 ]; then rm -f {output}.rapper.log; fi; "
        "exit $status"


rule interpro:
    """InterPro domain hierarchy, labels and GO mapping."""
    input:
        tree_file=config["interpro"]["tree_file"],
        entry_list=config["interpro"]["entry_list"],
        go_mapping=config["interpro"]["go_mapping"],
    output:
        temp(f"{OUTPUT_DIR}/interpro.raw.ttl"),
    shell:
        "python3 riceKB/interpro.py {input.tree_file} {input.entry_list} "
        "{input.go_mapping} {output}"


rule oryzabase:
    """OryzaBase gene list -- names, symbols, RAP-DB/MSU/Gramene cross-refs,
    GO/TO/PO classifications -- plus the side mapping files other
    AgroLD_ETL parsers (see riceKB/rapdbModel.py) join against.

    OryzaBase's gene list is exported on demand via a web form
    (https://shigen.nig.ac.jp/rice/oryzabase/gene/download), with a
    timestamp baked into the downloaded filename -- no stable URL to fetch
    automatically, so like interpro above, it's placed by hand under
    data_dir (rename or symlink the download to the path config.yaml
    expects).
    """
    input:
        config["oryzabase"]["gene_list"],
    output:
        genes=temp(f"{OUTPUT_DIR}/oryzabase_genes.raw.ttl"),
        pub_dict=f"{OUTPUT_DIR}/pub_dictionnary.txt",
        pub_dict_rapdb=f"{OUTPUT_DIR}/pub_dictionnary_rapdb.txt",
        pub_dict_msu=f"{OUTPUT_DIR}/pub_dictionnary_msu.txt",
        pub_dict_hub=f"{OUTPUT_DIR}/pub_genome_hub.txt",
        dict_po=f"{OUTPUT_DIR}/dict_PO_hub.txt",
        dict_to=f"{OUTPUT_DIR}/dict_TO_hub.txt",
    shell:
        "python3 riceKB/oryzaBaseParser.py {input} {output.genes} "
        "--pub-dict {output.pub_dict} --pub-dict-rapdb {output.pub_dict_rapdb} "
        "--pub-dict-msu {output.pub_dict_msu} --pub-dict-hub {output.pub_dict_hub} "
        "--dict-po {output.dict_po} --dict-to {output.dict_to}"


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
        tf=temp(f"{OUTPUT_DIR}/planttfdb_{{species}}.raw.ttl"),
        regulation=temp(f"{OUTPUT_DIR}/planttfdb_{{species}}_regulation.raw.ttl"),
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
        temp(f"{OUTPUT_DIR}/stringdb_{{species}}.raw.ttl"),
    params:
        taxon_id=lambda wc: AGROLD_TAXON_ID[wc.species],
    shell:
        "python3 riceKB/stringDB.py {input.aliases_file} {input.links_file} "
        "{output} --taxon-id {params.taxon_id}"


rule download_uniprot:
    """Fetch one shared UniProt taxonomic division flat-file (one per bank:
    sprot or trembl), still gzip-compressed, via bin/download_uniprot.py --
    not per species, unlike stringdb/planttfdb above: UniProt publishes ONE
    file per division, covering every species in it, and
    riceKB/uniprotToRDF.py filters the species it wants out of this shared
    file itself (see the uniprot rule below). Kept gzip-compressed on disk
    (not decompressed here) because uniprot_trembl_plants.dat.gz is ~14GB
    compressed -- decompressing it to a second on-disk copy isn't viable on
    a quota-limited cluster; riceKB/uniprotToRDF.py decompresses on the fly
    as it reads instead. Declared temp() so Snakemake deletes it once every
    uniprot job for this bank, across every taxon, has consumed it.
    """
    output:
        temp(os.path.join(config["uniprot"]["data_dir"], "uniprot_{bank}_plants.dat.gz")),
    wildcard_constraints:
        bank=r"sprot|trembl",
    shell:
        "python3 bin/download_uniprot.py --bank {wildcards.bank} --division plants "
        "-o {config[uniprot][data_dir]}"


rule uniprot:
    """UniProt protein records -- plus the keywords and PubMed citations
    they reference -- for one taxon and one bank (sprot or trembl),
    filtered out of the shared, still-gzipped plants division dump by
    riceKB/uniprotToRDF.py itself.

    {taxon} is a free wildcard here (not looked up via AGROLD_TAXON_ID,
    since output paths can't be computed by a function in Snakemake); it's
    rule all's expand() over AGROLD_TAXON_ID's values (UNIPROT_TAXA above)
    that actually drives which taxa get requested. A taxon with no entry
    for this bank in the plants division fails this job (see
    riceKB/uniprotToRDF.py's exit(1) on zero matches) -- run with
    `snakemake --keep-going` so that doesn't block the other taxa/banks.

    riceKB/uniprotToRDF.py builds its own output filenames internally from
    bank/taxon (uniprot_{taxon}_{bank}.plants.ttl etc. -- its CLI doesn't
    take explicit output paths, deliberately left unchanged), so unlike the
    other sources this rule can't just rename its own `output:` to
    *.raw.ttl and let the script write there directly; it renames the
    script's fixed-name files to *.raw.ttl itself afterwards instead, so
    validate_turtle still picks them up the same way as everywhere else.
    """
    input:
        os.path.join(config["uniprot"]["data_dir"], "uniprot_{bank}_plants.dat.gz"),
    output:
        protein=temp(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.plants.raw.ttl"),
        keyword=temp(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.keyword.raw.ttl"),
        pubmed=temp(f"{OUTPUT_DIR}/uniprot_{{taxon}}_{{bank}}.pubmed.raw.ttl"),
        gene_mapping=f"{OUTPUT_DIR}/mappings_{{taxon}}_{{bank}}_gene_prot.tsv",
    wildcard_constraints:
        taxon=r"\d+",
        bank=r"sprot|trembl",
    params:
        output_dir=OUTPUT_DIR,
    shell:
        "python3 riceKB/uniprotToRDF.py {wildcards.bank} {wildcards.taxon} "
        "{params.output_dir} {input} && "
        "mv {params.output_dir}/uniprot_{wildcards.taxon}_{wildcards.bank}.plants.ttl {output.protein} && "
        "mv {params.output_dir}/uniprot_{wildcards.taxon}_{wildcards.bank}.keyword.ttl {output.keyword} && "
        "mv {params.output_dir}/uniprot_{wildcards.taxon}_{wildcards.bank}.pubmed.ttl {output.pubmed}"
