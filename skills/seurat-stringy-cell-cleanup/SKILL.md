---
name: seurat-stringy-cell-cleanup
description: Diagnose Seurat UMAP stringy, bridge-like, trailing, or fragmented cell structures and prepare copy-first manual multi-region cleanup workflows. Use for Seurat/UMAP 拉丝 or bridge-like cell-removal requests, not for general plotting or unrelated single-cell analysis.
metadata:
  short-description: Seurat UMAP stringy-cell cleanup
---

# Seurat Stringy Cell Cleanup

Use this skill when a user asks about Seurat UMAP 拉丝, stringy, bridge-like, trailing, fragmented, or suspicious connecting cell structures, especially when they want diagnostic evidence or manual multi-region deletion.

This skill is extracted from `r-sci-plot-template` as a standalone workflow. Do not edit or rely on the parent skill when responding through this skill.

## Core Rules

- A UMAP bridge is not proof of a doublet, low-quality cell, ambient RNA, or a biological transition.
- Keep the original Seurat object, original RDS, and existing outputs read-only. Work from a copied object in a fresh run directory.
- For diagnosis or automatic filtering, read [references/seurat-umap-stringy-cell-removal.md](references/seurat-umap-stringy-cell-removal.md) before writing code.
- For explicit manual region deletion, use [scripts/seurat_multi_region_selector.R](scripts/seurat_multi_region_selector.R). Never substitute `CellSelector()` or a single `plotly_selected` result.
- Manual deletion records the reason as `user_selected_region`; do not claim a QC/doublet diagnosis unless diagnostic evidence was actually generated.
- For direct post-removal UMAP figures, preserve original UMAP coordinates and draw retained cells only. Do not recompute normalization, PCA, neighbors, clusters, or UMAP unless the user explicitly asks for a separate full reanalysis branch.
- Display generated plots during execution with `print()` and save audit-ready source tables.

## Manual Multi-Region Workflow

When the user explicitly wants to select and remove multiple regions:

1. Copy the input RDS into a fresh run directory and verify its checksum.
2. Copy `scripts/seurat_multi_region_selector.R` into the run's `scripts/` directory, source that copy, and log `MULTI_REGION_SELECTOR_VERSION`.
3. Launch the selector from an interactive R session or server RStudio browser/viewer. A headless SSH session alone cannot receive mouse selections.
4. Each lasso or rectangle release must append a numbered region without Shift or per-region confirmation. Keep numbered outlines visible, give each region a checkbox, and highlight/count the union of checked regions.
5. Delete only when the user clicks the delete button. Checkbox changes alone never perform deletion. Refuse empty and all-cell selections.
6. Save region membership, checked states, removed barcodes, removed metadata, per-cell decisions, counts, package versions, and a subset object.
7. Generate the direct deletion-effect UMAP from the original coordinates with selected cells removed from the plotted data.

Use [templates/stringy_cell_cleanup_copy_first.R](../../templates/stringy_cell_cleanup_copy_first.R) from this repository as the generic starter script when creating a project-specific cleanup script.

## Diagnostic Workflow

For diagnosis or evidence-supported filtering, follow the staged reference:

- copy-first setup and provenance;
- doublet audit by capture/library;
- QC overlays and project-specific cutoff provenance;
- ambient RNA or marker-contamination assessment;
- PC/UMAP sensitivity comparison as display diagnostics;
- biological-continuum review;
- candidate decision/export tables.

Automatic removal must stay within the user's requested scope and must preserve uncertain or plausible transitional cells unless the user chooses manual removal.

## Validation

For helper or template changes, use synthetic Seurat cells before touching real data. Verify rapid multi-region events, retries, delayed checkbox echoes, overlap de-duplication, empty/all-cell guards, undo/clear/cancel behavior, audit exports, and unchanged input serialization.

