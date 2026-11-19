-- Materializes <gene> agrold_vocabulary:encodes <protein>, joining data two
-- other converters already produce -- no new RDF source, just a SPARQL
-- INSERT/WHERE over what's already loaded:
--   - riceKB/uniprotToRDF.py emits <protein> agrold_vocabulary:orderedLocusName
--     "<symbol>" for every UniProt OrderedLocusNames entry.
--   - ensembl_rdf/rdf_converter_ensembl_db.py's rdfize_gene() emits
--     <gene> dcterms:identifier "<stable_id>" for every gene (line 349).
-- The UniProt symbol is (case-insensitively) the real Ensembl gene
-- stable_id -- see TODO.md for the empirical check (rice/Arabidopsis,
-- against real Ensembl core dumps) behind this join and behind choosing a
-- plain SPARQL INSERT/WHERE over SHACL-AF or a new Python converter.
--
-- bin/prepare_virtuoso_import.py's collect_agrold_etl_files() puts a
-- species' uniprot_<taxon>_*.ttl in the same named graph
-- (http://purl.agrold.org/graph/species/<taxon_id>) as that species' own
-- ensembl_rdf-derived genes, so the join and the INSERT both stay within a
-- single graph -- one query covers every species, none hardcoded; ?g is
-- bound per species graph by the WHERE clause itself.
--
-- `?gene a agrold_vocabulary:Gene` guards against ensembl_rdf's
-- rdfize_translation() also emitting dcterms:identifier (on proteins, with
-- their own distinct stable_id format) -- without this, a protein's own
-- identifier could in principle be mistaken for a gene's.
--
-- Run once, after bin/prepare_virtuoso_import.py's load.sql has fully
-- loaded both the species and uniprot graphs. Safe to re-run: inserting an
-- already-present triple is a no-op.
--
-- Usage: isql < bin/sparql/encodes_from_orderedLocusName.sql

SPARQL
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX agrold_vocabulary: <http://purl.agrold.org/vocabulary/>
INSERT {
  GRAPH ?g { ?gene agrold_vocabulary:encodes ?protein . }
}
WHERE {
  GRAPH ?g {
    ?gene a agrold_vocabulary:Gene ;
          dcterms:identifier ?stableId .
    ?protein agrold_vocabulary:orderedLocusName ?symbol .
    FILTER(STRSTARTS(STR(?g), "http://purl.agrold.org/graph/species/"))
    FILTER(UCASE(?stableId) = UCASE(?symbol))
  }
};

checkpoint;
