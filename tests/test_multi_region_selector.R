suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(shiny)
  library(jsonlite)
})

command_line <- commandArgs(FALSE)
file_arg <- grep("^--file=", command_line, value = TRUE)
test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path("tests", "test_multi_region_selector.R")
repo_root <- normalizePath(file.path(dirname(test_file), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(repo_root) || basename(repo_root) != "seurat-stringy-cell-cleanup") {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
helper_r <- file.path(repo_root, "skills", "seurat-stringy-cell-cleanup", "scripts", "seurat_multi_region_selector.R")
if (!file.exists(helper_r)) stop("Cannot find helper: ", helper_r)
source(helper_r, encoding = "UTF-8")
if (utils::compareVersion(MULTI_REGION_SELECTOR_VERSION, "2026-09-26.3") < 0L) {
  stop("Helper version is too old")
}

region <- function(ids) list(barcodes = paste0("cell", ids), vertices = NULL)
stroke <- function(id, ids) list(event_id = id, barcodes = paste0("cell", ids), vertices = NULL)
batch <- function(...) as.character(jsonlite::toJSON(list(...), auto_unbox = TRUE))
expect_error <- function(expr) {
  if (!inherits(tryCatch(force(expr), error = identity), "error")) {
    stop("Expected error was not raised")
  }
}

set.seed(20260927)
counts <- matrix(rpois(200 * 24, lambda = 5), nrow = 200, ncol = 24)
rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
obj <- CreateSeuratObject(counts = counts, project = "synthetic")
obj$group <- rep(c("A", "B", "C"), length.out = ncol(obj))
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 100, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:5, n.neighbors = 6, min.dist = 0.2, seed.use = 20260927, verbose = FALSE)

before <- serialize(obj, NULL)
valid_cells <- colnames(obj)
regions_static <- list(R001 = region(1:4), R002 = region(3:6), R003 = region(20:22))
stopifnot(identical(multi_region_cells(regions_static, c("R001", "R002"), valid_cells), paste0("cell", 1:6)))
expect_error(multi_region_cells(list(R001 = region(99)), "R001", valid_cells))
expect_error(multi_region_cells(regions_static, "missing", valid_cells))
expect_error(multi_region_commit(obj, regions_static, character(), tempdir(), "umap"))
expect_error(multi_region_commit(obj, list(R001 = region(1:24)), "R001", tempdir(), "umap"))

done <- list()
out_dir <- file.path(tempdir(), paste0("selector_test_", Sys.getpid()))
app <- build_multi_region_cell_app(obj, out_dir, group.by = "group",
  on_done = function(value) { done[[length(done) + 1L]] <<- list(value) })

shiny::testServer(app, {
  session$flushReact()
  session$setInputs(region_events = batch(stroke("a", 1:4)))
  stopifnot(length(regions()) == 1L, identical(checked_ids(), "R001"), length(selected()) == 4L)

  session$setInputs(region_events = batch(stroke("a", 1:4), stroke("b", 3:6), stroke("c", 20:22)))
  stopifnot(length(regions()) == 3L, length(checked_ids()) == 3L, length(selected()) == 9L)

  session$setInputs(checked_regions = "R001")
  stopifnot(length(checked_ids()) == 3L, length(selected()) == 9L)

  session$setInputs(region_checks = list(visible = c("R001", "R002"), checked = "R002"))
  stopifnot(setequal(checked_ids(), c("R002", "R003")), length(selected()) == 7L)

  session$setInputs(region_events = batch(stroke("d", 8:9)))
  stopifnot(length(regions()) == 4L, !"R001" %in% checked_ids(), "R004" %in% checked_ids())

  session$setInputs(region_events = batch(stroke("e", integer()), stroke("f", 99)))
  stopifnot(length(regions()) == 4L)

  session$setInputs(undo_region = 1)
  stopifnot(length(regions()) == 3L, length(selected()) == 7L)

  session$setInputs(select_none = 1)
  stopifnot(length(selected()) == 0L)
  session$setInputs(delete_checked = 1)
  stopifnot(!finished)

  session$setInputs(select_all = 1)
  stopifnot(length(selected()) == 9L)
  session$setInputs(region_checks = list(visible = names(regions()), checked = c("R001", "R002")))
  stopifnot(length(selected()) == 6L)
  session$setInputs(delete_checked = 2)
})

stopifnot(length(done) == 1L)
result <- done[[1]][[1]]
stopifnot(identical(result$removed_cells, paste0("cell", 1:6)))
stopifnot(ncol(result$object) == 18L)
stopifnot(identical(before, serialize(obj, NULL)))
stopifnot(isTRUE(result$recompute_required))

decision <- read.csv(file.path(result$export_dir, "cell_decisions.csv"), stringsAsFactors = FALSE)
stopifnot(identical(decision$barcode[decision$decision == "remove"], result$removed_cells))
stopifnot(identical(readRDS(result$subset_rds), result$object))

cancelled <- FALSE
app_cancel <- build_multi_region_cell_app(obj, file.path(tempdir(), paste0("selector_cancel_", Sys.getpid())),
  on_done = function(value) { stopifnot(is.null(value)); cancelled <<- TRUE })
shiny::testServer(app_cancel, {
  session$flushReact()
  session$setInputs(region_events = batch(stroke("all", 1:24)))
  session$setInputs(delete_checked = 1)
  stopifnot(!finished)
  session$setInputs(clear_regions = 1)
  stopifnot(length(regions()) == 0L, length(selected()) == 0L)
  session$setInputs(region_events = batch(stroke("all", 1:24)))
  stopifnot(length(regions()) == 0L)
  session$setInputs(cancel = 1)
})
stopifnot(cancelled)
stopifnot(identical(before, serialize(obj, NULL)))

cat("ALL_SYNTHETIC_MULTI_REGION_SELECTOR_TESTS_PASSED\n")
