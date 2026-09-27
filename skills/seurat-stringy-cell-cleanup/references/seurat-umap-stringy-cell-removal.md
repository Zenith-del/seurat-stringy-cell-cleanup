# Seurat UMAP Stringy Or Bridge-Like Cell Cleanup

Use this reference for Seurat UMAP structures described as stringy, bridge-like, trailing, stretched, fragmented, or 拉丝. Use the diagnostic stages for diagnosis or automatic filtering. For an explicit manual-deletion request, go directly to the multi-region interface after copy-first setup; diagnostics remain available but do not block the user's chosen manual action.

## Scope And Copy-First Setup

The original RDS, in-memory object, and existing outputs are read-only. Copy first, work in a fresh run directory, and record input path, checksum, package versions, assay, reduction, capture column, and parameters. Do not assign a cleaned object over the user's original variable.

For file-backed or on-disk assays, inspect backing paths before any operation that might write to storage. An RDS copy alone may still reference an external matrix store.

Minimal copy-first pattern:

```r
suppressPackageStartupMessages(library(Seurat))
input_rds <- Sys.getenv("STRINGY_INPUT_RDS")
base_dir <- Sys.getenv("STRINGY_OUT_DIR", unset = "stringy_cell_cleanup_runs")
stopifnot(nzchar(input_rds), file.exists(input_rds))
run_dir <- file.path(base_dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_stringy_cleanup"))
dir.create(run_dir, recursive = TRUE)
for (part in c("scripts", "rds", "plots", "tables", "logs")) {
  dir.create(file.path(run_dir, part), recursive = TRUE)
}
new_path <- function(part, filename) {
  path <- file.path(run_dir, part, filename)
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  path
}
input_md5 <- unname(tools::md5sum(input_rds))
working_rds <- new_path("rds", "input_working_copy.rds")
stopifnot(file.copy(input_rds, working_rds, overwrite = FALSE))
stopifnot(identical(input_md5, unname(tools::md5sum(working_rds))))
obj_work <- readRDS(working_rds)
stopifnot(inherits(obj_work, "Seurat"))
```

## Diagnostic Stages

### 1. Audit Doublets

Inspect existing doublet results and provenance before rerunning detectors. Reuse valid calls on a singlet-filtered object. Do not repeatedly run a detector on successively cleaned subsets until bridges disappear.

If rerunning `scDblFinder`, group by independent capture/library, not by treatment group. Use observed raw counts, not integrated/scaled/log-normalized values, and do not densify the matrix. Record version, seed, grouping column, and rate assumptions.

Compare candidate and non-candidate cells within capture and major cell type for doublet class, doublet score, nCount/nFeature, and co-expression of incompatible markers. Mixed markers support a hypothesis but can be mimicked by ambient RNA.

### 2. Overlay QC Metrics

Plot available QC features such as `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`, and any project-specific hemoglobin or ribosomal metrics. No single direction is universal; nuclei and cells can differ. Report missing metrics instead of fabricating zero-valued fields.

Reuse saved QC cutoffs when they exist. If they do not exist, state which input population defines the cutoff and save it before filtering.

### 3. Assess Ambient RNA Or Marker Contamination

Inspect suspected ambient genes, contamination metrics, and mixed marker expression using RNA counts or normalized RNA expression. Do not use integrated or scaled values for this audit. Preserve raw counts and store corrected counts in a new assay/object when correction is requested and justified.

### 4. Compare PCs Before Changing UMAP Geometry

Compare PC choices on identical cells, normalization, features, reduction, and seeds. Display comparisons are not deletion evidence. If evaluating clustering, build separate neighbors/graphs and cluster columns rather than overwriting current results.

### 5. Tune UMAP As Display Only

UMAP parameters can make structures look more or less compact. Changing `min.dist`, `n.neighbors`, spread, metric, or seeds does not repair doublets, contamination, or low-quality cells. Do not choose a final biological analysis solely because a UMAP looks cleaner.

### 6. Evaluate Biological Continuity

Retain plausible transitions when QC is normal, doublet evidence is weak, endpoints share a major cell type, markers change coherently, and the pattern is reproducible across captures. UMAP alone does not establish direction or pseudotime.

## Manual Multi-Region Selection And Direct Removal

Use this mode when the user explicitly asks to select and delete several regions. Copy or source `scripts/seurat_multi_region_selector.R`. Required packages are `SeuratObject`, `shiny`, `plotly`, `htmlwidgets`, and `jsonlite`.

The interface must support this exact interaction:

1. Each lasso or rectangle release appends a numbered region. Users continue drawing without Shift or per-region confirmation.
2. Saved outlines and IDs remain visible. New regions start checked; previous check choices are preserved.
3. Checked regions are highlighted and their barcode union is counted. Overlapping cells count once. Unchecked regions contribute no cells unless the same cell is in another checked region.
4. Select-all, select-none, undo-last-region, and clear-region edit selection only.
5. Clicking the delete button subsets a separate object and exports audit files. Empty and all-cell selections are refused. Cancel or closing the app returns no result.

Use full-resolution coordinates and barcode-backed Plotly `customdata`. Never infer barcodes from trace indices or `pointNumber` across traces. The saved barcode union is authoritative; geometry is provenance.

The helper exports a unique `manual_regions_*` directory containing region summaries, region membership, checked states, removed barcodes, removed metadata, per-cell decisions, counts, `sessionInfo()`, and a subset object. Use reason `user_selected_region`.

## Post-Removal Visualization

For a direct deletion-effect UMAP:

- preserve original UMAP coordinates;
- remove selected barcodes from the plotted data;
- label the figure as using original coordinates;
- save a plot source table with barcode, original UMAP coordinates, removal status, and coloring metadata;
- do not save or name outputs as recomputed unless a separate reanalysis branch actually reruns the reduction workflow.

Only run normalization, PCA, neighbors, clusters, or UMAP again when the user explicitly asks for a full reanalysis branch. Label that branch separately.

## Deliverables

- Fresh `scripts/`, `plots/`, `tables/`, `rds/`, and `logs` directories.
- Input provenance and unchanged checksum.
- Candidate, retained, review, and removed barcode tables with reasons.
- Before/selected/retained UMAP plots and source data.
- Diagnostic plots and tables when requested.
- Explicit notes for any missing packages, missing metrics, or unrun validation steps.

## Validation

Before touching real data, use synthetic Seurat cells. Verify rapid consecutive regions, retries, delayed checkbox echoes, overlap de-duplication, all-unchecked and all-cell guards, undo/clear/cancel behavior, audit-table consistency, and unchanged input serialization/checksum.

## Primary References

- Shiny checkbox groups, `stopApp`, and `testServer`.
- Plotly selection events and barcode-backed `customdata`.
- `scDblFinder` vignette for capture-aware doublet detection.
- `decontX` vignette for contamination correction.
- Seurat essential commands and `RunUMAP` documentation.

