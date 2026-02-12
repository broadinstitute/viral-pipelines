# Genomad Database Test Fixture

The genomad_single workflow requires a genomad database for testing. Due to size constraints (706MB compressed), this fixture is not stored in git.

## Quick Setup

If you have access to `/data/databases/genomad/genomad_db`:

```bash
cd /path/to/viral-pipelines
tar -cf - -C /data/databases/genomad genomad_db | zstd -19 -o test/input/genomad_db-test.tar.zst
```

## Alternative: Download from Zenodo

1. Download the official genomad database from https://zenodo.org/records/10594875
2. Extract and re-compress:

```bash
# Download
wget https://zenodo.org/records/10594875/files/genomad_db_v1.9.tar.gz

# Extract
tar -xzf genomad_db_v1.9.tar.gz

# Re-compress with zstd
tar -cf - genomad_db | zstd -19 -o test/input/genomad_db-test.tar.zst

# Clean up
rm -rf genomad_db genomad_db_v1.9.tar.gz
```

## Verify Setup

```bash
# Should show ~706MB file
ls -lh test/input/genomad_db-test.tar.zst

# Test the workflow
miniwdl run pipes/WDL/workflows/genomad_single.wdl \
  -i test/input/WDL/test_inputs-genomad_single-local.json
```

## File Details

- **Uncompressed size:** 1.4GB
- **Compressed size (zstd -19):** 706MB
- **Database version:** 1.9
- **Contains:** MMseqs2 databases, taxonomy files, marker metadata
