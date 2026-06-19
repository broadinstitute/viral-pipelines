# LucaVirus Capsid Test Panel

Small UniProt-derived protein FASTA panel for manual LucaVirus
`viral_capsid` smoke testing.

Files:

- `capsid_positive.uniprot.faa`: reviewed viral capsid proteins expected to be
  positive controls for the `viral_capsid` profile.
- `viral_non_capsid.uniprot.faa`: reviewed viral polymerase/replication
  proteins used as viral negative controls for the `viral_capsid` profile.
- `nonviral.uniprot.faa`: reviewed human, yeast, and bacterial proteins used as
  non-viral negative controls.
- `lucavirus_capsid_test_panel.uniprot.faa`: combined FASTA for direct WDL
  input.
- `metadata.tsv`: accession-level category labels and UniProt source URLs.

The FASTA records keep original UniProt headers. Use `metadata.tsv` to map
accessions back to expected test categories.

Run with:

```bash
cromwell -Dconfig.file=pipes/cromwell/cromwell.local-gpu.conf \
  run pipes/WDL/workflows/classify_lucavirus.wdl \
  -i test/input/WDL/cromwell-gpu/test_inputs-classify_lucavirus-capsid_panel-local_gpu.json \
  -m lucavirus_capsid_panel.metadata.json
```
