# Setup

## 1. Environment

```bash
conda env create -f environment.yml
conda activate resko-network
pip install -e ".[dev]"
```

## 2. Verify

```bash
pytest                     # 70 tests, ~4 s. All should pass.
```

If they don't, stop here — the tests pin the three corrections, so a failure
means the corrections aren't active.

## 3. VS Code

```bash
mv dot-vscode .vscode
```

(The directory ships as `dot-vscode` because some tooling refuses to write
dot-directories.)

Then Cmd-Shift-P → *Python: Select Interpreter* → `resko-network`. The Testing
panel picks up `tests/` automatically.

Recommended extensions: `ms-python.python`, `charliermarsh.ruff`.

## 4. Point it at your existing data

Nothing needs refetching. The ported stages read the R pipeline's outputs:

```bash
resko rank      /path/to/RESKO-Network/results/A12B_compound_protein_evidence_corrected.csv
resko audit-cache /path/to/RESKO-Network/results/A17D_api_cache
resko coverage  /path/to/RESKO-Network/eef1a_master_seed_list.csv \
                /path/to/RESKO-Network/SIDER/drug_names.tsv
```

`resko rank --compare` shows the corrected ranking against the R original's.
