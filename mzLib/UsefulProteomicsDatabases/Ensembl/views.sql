-- Views over an ensembl-orthology-snapshot (format 1), as DuckDB table macros. Nothing here is
-- stored: each view is a filter or join over the snapshot's Parquet files, and each matches the
-- member of mzLib's OrthologySnapshot with the same name.
--
-- Load once per connection, then pass the snapshot directory as the first argument:
--   .read compara-116/views.sql
--   SELECT status, count(*) FROM pair_status('compara-116', 'homo_sapiens', 'mus_musculus') GROUP BY status;

-- The pair file holding species a and b, whichever order they are given in.
CREATE OR REPLACE MACRO orthology_pair_file(root, a, b) AS
    root || '/pairs/' || least(a, b) || '__' || greatest(a, b) || '.parquet';

-- The ortholog rows between two different species, oriented so gene is species a's.
CREATE OR REPLACE MACRO orthologs(root, a, b) AS TABLE
    SELECT homology_id, relationship_type,
           CASE WHEN species_a = a THEN gene_a ELSE gene_b END AS gene,
           CASE WHEN species_a = a THEN protein_a ELSE protein_b END AS protein,
           CASE WHEN species_a = a THEN identity_a ELSE identity_b END AS identity,
           CASE WHEN species_a = a THEN gene_b ELSE gene_a END AS partner_gene,
           CASE WHEN species_a = a THEN protein_b ELSE protein_a END AS partner_protein,
           CASE WHEN species_a = a THEN identity_b ELSE identity_a END AS partner_identity,
           dn, ds, goc_score, wga_coverage, is_high_confidence
    FROM read_parquet(orthology_pair_file(root, a, b))
    WHERE relationship_class = 'ortholog' AND a <> b;

-- Every gene of species a, all biotypes, with exactly one status in species b. Which genes are the
-- denominator is the caller's choice (filter on gene_biotype); it is not the same as protein_coding.
CREATE OR REPLACE MACRO pair_status(root, a, b) AS TABLE
    WITH o AS (SELECT DISTINCT gene FROM orthologs(root, a, b)),
         ma AS (SELECT gene_id, group_id FROM read_parquet(root || '/members/' || a || '.parquet')),
         tb AS (SELECT DISTINCT group_id FROM read_parquet(root || '/members/' || b || '.parquet'))
    SELECT g.gene_id, g.gene_biotype,
           CASE WHEN o.gene IS NOT NULL THEN 'has_ortholog'
                WHEN ma.group_id IS NULL THEN 'not_in_any_tree'
                WHEN tb.group_id IS NOT NULL THEN 'no_edge_in_shared_tree'
                ELSE 'tree_lacks_target_species' END AS status
    FROM read_parquet(root || '/genes/' || a || '.parquet') g
    LEFT JOIN o ON o.gene = g.gene_id
    LEFT JOIN ma ON ma.gene_id = g.gene_id
    LEFT JOIN tb ON tb.group_id = ma.group_id
    WHERE a <> b;

-- One gene per species where all three pairs are orthologs; never chained. For more species, use
-- OrthologySnapshot.SpeciesSet in mzLib, or join one more orthologs() per added species the same way.
CREATE OR REPLACE MACRO species_set3(root, a, b, c) AS TABLE
    SELECT ab.gene AS gene_1, ab.partner_gene AS gene_2, ac.partner_gene AS gene_3,
           (ab.relationship_type = 'ortholog_one2one' AND ac.relationship_type = 'ortholog_one2one'
            AND bc.relationship_type = 'ortholog_one2one') AS all_one2one
    FROM orthologs(root, a, b) ab
    JOIN orthologs(root, a, c) ac ON ac.gene = ab.gene
    JOIN orthologs(root, b, c) bc ON bc.gene = ab.partner_gene AND bc.partner_gene = ac.partner_gene;
