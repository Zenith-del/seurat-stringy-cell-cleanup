# Validation Notes

Validated release baseline:

- Skill name: `seurat-stringy-cell-cleanup`
- Helper version: `2026-09-26.3`
- Parent skill left unchanged: `r-sci-plot-template`
- Real data included: no
- Local/server-specific paths included in packaged files: no intended defaults

Checks to run before publishing:

1. Run Codex skill validation on `skills/seurat-stringy-cell-cleanup`.
2. Parse all R scripts.
3. Run `tests/test_multi_region_selector.R` in an R environment with Seurat, Shiny, Plotly, htmlwidgets, and jsonlite.
4. Search for private paths, credentials, and large outputs:

```powershell
Select-String -Path .\**\* -Pattern "/home/|D:/ICH|C:\\|password|token|secret|credential|\.rds|\.csv|\.png|\.pdf" -CaseSensitive:$false
```

Expected synthetic-test behaviors:

- Three or more rapid region events can be retained before checkbox echoes arrive.
- Retry events are deduplicated.
- User unchecks persist when a new region arrives.
- Overlapping regions remove each barcode once.
- Empty selections and all-cell selections cannot delete.
- Undo, clear, and cancel do not delete cells.
- The original object serialization remains unchanged.
- Audit exports match the returned subset object.

