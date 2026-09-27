# Generic copy-first Seurat UMAP stringy-cell cleanup template.
# It preserves the original input file and draws direct deletion-effect UMAPs
# using original coordinates. It does not recompute UMAP unless you create a
# separate, explicitly requested reanalysis branch.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

# 00 Setup ----

truthy <- function(x) tolower(x) %in% c("1", "true", "yes", "y")

input_rds <- Sys.getenv("STRINGY_INPUT_RDS", unset = "")
out_root <- Sys.getenv("STRINGY_OUT_DIR", unset = "stringy_cell_cleanup_runs")
reduction_name <- Sys.getenv("STRINGY_REDUCTION", unset = "umap")
group_by_env <- Sys.getenv("STRINGY_GROUP_BY", unset = "")
shiny_host <- Sys.getenv("STRINGY_SHINY_HOST", unset = "127.0.0.1")
shiny_port_env <- Sys.getenv("STRINGY_SHINY_PORT", unset = "")
launch_browser <- truthy(Sys.getenv("STRINGY_LAUNCH_BROWSER", unset = if (interactive()) "true" else "false"))

if (!nzchar(input_rds) || !file.exists(input_rds)) {
  stop("Set STRINGY_INPUT_RDS to an existing Seurat RDS file.", call. = FALSE)
}
shiny_port <- if (nzchar(shiny_port_env)) as.integer(shiny_port_env) else NULL
if (length(shiny_port) && (is.na(shiny_port) || shiny_port <= 0)) {
  stop("STRINGY_SHINY_PORT must be a positive integer.", call. = FALSE)
}

script_path <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE) else {
    frame_files <- vapply(sys.frames(), function(frame) {
      if (!is.null(frame$ofile)) as.character(frame$ofile)[1] else NA_character_
    }, character(1))
    frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
    if (length(frame_files)) normalizePath(tail(frame_files, 1), winslash = "/", mustWork = FALSE) else NA_character_
  }
}
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_helper <- normalizePath(
  file.path(script_dir, "..", "skills", "seurat-stringy-cell-cleanup", "scripts", "seurat_multi_region_selector.R"),
  winslash = "/",
  mustWork = FALSE
)
same_dir_helper <- file.path(script_dir, "seurat_multi_region_selector.R")
helper_r <- Sys.getenv(
  "STRINGY_HELPER_R",
  unset = if (file.exists(same_dir_helper)) same_dir_helper else repo_helper
)
if (!file.exists(helper_r)) {
  stop("Cannot find selector helper. Set STRINGY_HELPER_R or copy seurat_multi_region_selector.R beside this script.", call. = FALSE)
}

make_fresh_dir <- function(root, prefix) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(root, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", prefix))
  candidate <- stem
  i <- 1L
  while (dir.exists(candidate)) {
    candidate <- paste0(stem, "_", sprintf("%02d", i))
    i <- i + 1L
  }
  if (!dir.create(candidate, recursive = FALSE, showWarnings = FALSE)) {
    stop("Cannot create run directory: ", candidate, call. = FALSE)
  }
  candidate
}

run_dir <- make_fresh_dir(file.path(out_root, "runs"), "manual_stringy_cleanup")
for (part in c("scripts", "rds", "plots", "tables", "logs")) {
  dir.create(file.path(run_dir, part), recursive = TRUE, showWarnings = FALSE)
}
new_path <- function(part, filename) {
  path <- file.path(run_dir, part, filename)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  path
}
write_csv <- function(data, filename, row.names = FALSE) {
  path <- new_path("tables", filename)
  utils::write.csv(data, path, row.names = row.names)
  message("Saved table: ", path)
  invisible(path)
}
save_plot <- function(plot, stem, width = 7.5, height = 6, dpi = 320) {
  print(plot)
  pdf_path <- new_path("plots", paste0(stem, ".pdf"))
  png_path <- new_path("plots", paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  invisible(c(pdf = pdf_path, png = png_path))
}

log_file <- file.path(run_dir, "logs", "stringy_cleanup.log")
log_con <- file(log_file, open = "at")
sink(log_file, split = TRUE)
sink(log_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

message("Run directory: ", run_dir)

# 01 Copy input and helper ----

input_md5 <- unname(tools::md5sum(input_rds))
working_rds <- new_path("rds", paste0(tools::file_path_sans_ext(basename(input_rds)), "_working_copy.rds"))
if (!file.copy(input_rds, working_rds, overwrite = FALSE)) {
  stop("Failed to copy input RDS.", call. = FALSE)
}
stopifnot(identical(input_md5, unname(tools::md5sum(working_rds))))

helper_run <- new_path("scripts", "seurat_multi_region_selector.R")
if (!file.copy(helper_r, helper_run, overwrite = FALSE)) {
  stop("Failed to copy selector helper.", call. = FALSE)
}
if (!is.na(script_path) && file.exists(script_path)) {
  file.copy(script_path, new_path("scripts", basename(script_path)), overwrite = FALSE)
}

selector_env <- new.env(parent = environment())
sys.source(helper_run, envir = selector_env)
selector_version <- selector_env$MULTI_REGION_SELECTOR_VERSION
if (is.null(selector_version) || utils::compareVersion(selector_version, "2026-09-26.3") < 0L) {
  stop("Selector helper is outdated; require version 2026-09-26.3 or newer.", call. = FALSE)
}

write_csv(data.frame(
  input_rds_read_only = input_rds,
  input_md5_before = input_md5,
  working_copy_rds = working_rds,
  helper_source = helper_r,
  helper_run_copy = helper_run,
  selector_version = selector_version,
  run_dir = run_dir,
  stringsAsFactors = FALSE
), "input_provenance.csv")

# 02 Load object and plot before selection ----

obj_work <- readRDS(working_rds)
if (!inherits(obj_work, "Seurat")) stop("Input RDS is not a Seurat object.", call. = FALSE)
if (!reduction_name %in% SeuratObject::Reductions(obj_work)) {
  stop("Reduction not found in object: ", reduction_name, call. = FALSE)
}
if (anyDuplicated(colnames(obj_work))) stop("Cell barcodes are duplicated.", call. = FALSE)

meta <- obj_work[[]]
preferred_group <- c("seurat_clusters", grep("_snn_res\\.", colnames(meta), value = TRUE), "sample_id", "orig.ident")
preferred_group <- preferred_group[preferred_group %in% colnames(meta)]
group_by <- if (nzchar(group_by_env) && !identical(tolower(group_by_env), "none")) {
  group_by_env
} else if (length(preferred_group)) {
  preferred_group[1]
} else {
  NULL
}
if (!is.null(group_by) && !group_by %in% colnames(meta)) {
  stop("STRINGY_GROUP_BY column not found: ", group_by, call. = FALSE)
}

p_before <- DimPlot(obj_work, reduction = reduction_name, group.by = group_by, pt.size = 0.2) +
  ggplot2::ggtitle("Before manual stringy-cell selection") +
  ggplot2::theme_classic(base_size = 11)
save_plot(p_before, "01_umap_before_manual_selection")

# 03 Manual multi-region selection ----

message("Launching selector. Draw lasso/rectangle regions; checked regions are removed only after clicking the delete button.")
selection <- selector_env$multi_region_cell_selector(
  object = obj_work,
  out_dir = file.path(run_dir, "tables"),
  reduction = reduction_name,
  group.by = group_by,
  launch.browser = launch_browser,
  host = shiny_host,
  port = shiny_port
)
if (is.null(selection)) stop("Selection cancelled; no cells removed.", call. = FALSE)

removed_cells <- selection$removed_cells
obj_subset <- selection$object
stopifnot(length(removed_cells) > 0, length(removed_cells) < ncol(obj_work))
stopifnot(!any(removed_cells %in% colnames(obj_subset)))
stopifnot(ncol(obj_subset) == ncol(obj_work) - length(removed_cells))

input_md5_after <- unname(tools::md5sum(input_rds))
if (!identical(input_md5, input_md5_after)) {
  stop("Input RDS checksum changed during workflow.", call. = FALSE)
}

subset_rds <- new_path("rds", "after_manual_region_removal_subset_only_original_coordinates.rds")
saveRDS(obj_subset, subset_rds)
writeLines(removed_cells, new_path("tables", "removed_barcodes_manual_regions.txt"))

decision_table <- data.frame(
  barcode = colnames(obj_work),
  decision = ifelse(colnames(obj_work) %in% removed_cells, "remove", "retain"),
  reason = ifelse(colnames(obj_work) %in% removed_cells, "user_selected_region", "outside_checked_regions"),
  stringsAsFactors = FALSE
)
write_csv(decision_table, "manual_region_cell_decisions.csv")

# 04 Direct deletion-effect UMAP using original coordinates ----

embedding <- SeuratObject::Embeddings(obj_work, reduction = reduction_name)[colnames(obj_work), 1:2, drop = FALSE]
plot_df <- data.frame(
  barcode = colnames(obj_work),
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  removed_status = ifelse(colnames(obj_work) %in% removed_cells, "removed", "retained"),
  stringsAsFactors = FALSE
)
if (!is.null(group_by)) {
  plot_df$group <- as.character(obj_work[[]][plot_df$barcode, group_by, drop = TRUE])
} else {
  plot_df$group <- as.character(SeuratObject::Idents(obj_work)[plot_df$barcode])
}
plot_df$group[is.na(plot_df$group)] <- "unknown"
write_csv(plot_df, "umap_original_coordinates_with_removal_status.csv")

retained_df <- plot_df[plot_df$removed_status == "retained", , drop = FALSE]
p_retained <- ggplot2::ggplot(retained_df, ggplot2::aes(UMAP_1, UMAP_2, color = group)) +
  ggplot2::geom_point(size = 0.2, alpha = 0.8) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "After manual region removal (original UMAP coordinates)",
    subtitle = "Selected cells are omitted from the plotted data; UMAP was not recomputed.",
    color = ifelse(is.null(group_by), "identity", group_by)
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot(p_retained, "02_umap_after_manual_removal_original_coordinates")

highlight_df <- plot_df
p_highlight <- ggplot2::ggplot(highlight_df, ggplot2::aes(UMAP_1, UMAP_2)) +
  ggplot2::geom_point(ggplot2::aes(color = removed_status), size = 0.2, alpha = 0.85) +
  ggplot2::scale_color_manual(values = c(retained = "grey75", removed = "#F26D4F")) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Manual regions selected for removal",
    subtitle = "Orange cells are user-selected; this is not an automatic QC diagnosis.",
    color = "Status"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot(p_highlight, "03_umap_selected_cells_highlighted_original_coordinates")

write_csv(data.frame(
  input_rds_read_only = input_rds,
  input_md5_before = input_md5,
  input_md5_after = input_md5_after,
  subset_rds = subset_rds,
  helper_selection_export_dir = selection$export_dir,
  cells_before = ncol(obj_work),
  removed_cells = length(removed_cells),
  cells_after_subset = ncol(obj_subset),
  original_checksum_unchanged = identical(input_md5, input_md5_after),
  umap_policy = "original coordinates retained; selected cells omitted from plotted data",
  reason = "user_selected_region",
  stringsAsFactors = FALSE
), "manual_region_removal_final_summary.csv")

writeLines(capture.output(utils::sessionInfo()), new_path("logs", "sessionInfo.txt"))
message("Done. Run directory: ", run_dir)

