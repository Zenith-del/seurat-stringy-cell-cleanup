# seurat-stringy-cell-cleanup

Codex skill and R templates for diagnosing Seurat UMAP stringy, bridge-like, trailing, or fragmented cell structures, plus an interactive multi-region selector for manual cleanup.

This repository is a clean extraction from a larger plotting skill. It does not modify the parent `r-sci-plot-template` skill and does not include real study data, run outputs, server paths, or project-specific results.

## What Is Included

- `skills/seurat-stringy-cell-cleanup/`: the standalone Codex skill.
- `skills/seurat-stringy-cell-cleanup/scripts/seurat_multi_region_selector.R`: Shiny/Plotly multi-region selector helper, version `2026-09-26.3`.
- `skills/seurat-stringy-cell-cleanup/references/seurat-umap-stringy-cell-removal.md`: detailed diagnostic and manual-selection guidance.
- `templates/stringy_cell_cleanup_copy_first.R`: generic copy-first R template.
- `tests/test_multi_region_selector.R`: synthetic non-real-data test script.
- `docs/VALIDATION.md`: validation notes and release checks.

## Safety Model

UMAP bridge geometry is not, by itself, evidence of doublets, low quality cells, ambient RNA, or a biological transition. Diagnostic or automatic filtering should audit doublets, QC metrics, ambient RNA, PC/UMAP sensitivity, and biological continuity before removing cells.

Manual removal is different: when a user explicitly selects regions and clicks the delete button, the recorded reason is `user_selected_region`. The helper subsets a separate working object, writes audit tables, and leaves the original object and original file untouched.

For the direct deletion-effect figure, the template preserves original UMAP coordinates and redraws retained cells only. It does not recompute normalization, PCA, neighbors, clusters, or UMAP unless a user explicitly asks for a separate full reanalysis branch.

## Requirements

R packages:

- `Seurat`
- `SeuratObject`
- `ggplot2`
- `shiny`
- `plotly`
- `htmlwidgets`
- `jsonlite`

Optional diagnostic packages depend on the requested workflow, for example `scDblFinder`, `SingleCellExperiment`, `BiocParallel`, or `decontX`.

## Use The Template

Set at least the input path, then run the template in RStudio or another interactive R session:

```r
Sys.setenv(STRINGY_INPUT_RDS = "/path/to/seurat_object.rds")
Sys.setenv(STRINGY_OUT_DIR = "/path/to/stringy_cleanup_runs")
source("templates/stringy_cell_cleanup_copy_first.R")
```

Optional overrides:

```r
Sys.setenv(STRINGY_REDUCTION = "umap")
Sys.setenv(STRINGY_GROUP_BY = "seurat_clusters")
Sys.setenv(STRINGY_HELPER_R = "skills/seurat-stringy-cell-cleanup/scripts/seurat_multi_region_selector.R")
Sys.setenv(STRINGY_SHINY_HOST = "127.0.0.1")
Sys.setenv(STRINGY_SHINY_PORT = "4877")
```

The selector lets you draw several lasso or rectangle regions without holding Shift. Each completed draw creates a numbered saved region. Checked regions are highlighted and removed only when you click the delete button.

## Install As A Codex Skill

Copy `skills/seurat-stringy-cell-cleanup` into your Codex skills directory:

```powershell
Copy-Item -Recurse -Force .\skills\seurat-stringy-cell-cleanup $env:USERPROFILE\.codex\skills\
```

Then invoke it as:

```text
Use $seurat-stringy-cell-cleanup to diagnose and handle UMAP stringy cells in my Seurat object.
```

## Validate Locally

Validate the skill metadata:

```powershell
$env:PYTHONUTF8 = "1"
python -X utf8 "$env:USERPROFILE\.codex\skills\.system\skill-creator\scripts\quick_validate.py" ".\skills\seurat-stringy-cell-cleanup"
```

Parse R scripts:

```powershell
Rscript -e "parse('skills/seurat-stringy-cell-cleanup/scripts/seurat_multi_region_selector.R'); parse('templates/stringy_cell_cleanup_copy_first.R'); parse('tests/test_multi_region_selector.R')"
```

Run synthetic tests when the R dependencies are installed:

```powershell
Rscript tests/test_multi_region_selector.R
```

## Release Hygiene

Before pushing to GitHub, keep real `.rds`, `.csv`, `.png`, `.pdf`, logs, run directories, server sync scripts, and credentials out of the repository. The `.gitignore` is intentionally strict for analysis outputs.

