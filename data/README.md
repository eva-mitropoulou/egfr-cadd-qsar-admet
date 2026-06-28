# Local Data

Raw and processed EGFR activity tables are intentionally kept local.

The scripts fetch public ChEMBL EGFR IC50 records and build the processed
tables used for descriptors, fingerprints, validation, and ranking. I keep the
generated tables out of git because they include record-level molecule
representations and can be regenerated.

The repository does keep the source code, final reports, metric JSON files,
figures, tests, and small structure/redocking artifacts needed to inspect the
workflow without publishing every intermediate table.
