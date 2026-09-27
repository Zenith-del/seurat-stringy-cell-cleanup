# Source this helper into a local environment; it does not run an analysis.
# Requires SeuratObject, shiny, plotly, htmlwidgets, and jsonlite. Never installs packages.
MULTI_REGION_SELECTOR_VERSION <- "2026-09-26.3"

multi_region_cells <- function(regions, checked, valid_cells) {
  checked <- as.character(checked)
  if (anyNA(checked) || !all(checked %in% names(regions))) {
    stop("Unknown region ID")
  }
  cells <- unique(as.character(unlist(lapply(regions[checked], `[[`, "barcodes"),
                                      use.names = FALSE)))
  if (anyNA(cells) || !all(cells %in% valid_cells)) stop("Unknown cell barcode")
  valid_cells[valid_cells %in% cells]
}

multi_region_commit <- function(object, regions, checked, out_dir, reduction) {
  all_cells <- colnames(object)
  removed <- multi_region_cells(regions, checked, all_cells)
  if (!length(removed)) stop("Select at least one region containing cells.")
  if (length(removed) == length(all_cells)) stop("The current selection contains all cells; narrow the deletion range.")
  kept <- all_cells[!all_cells %in% removed]
  clean <- subset(object, cells = kept)
  stopifnot(identical(colnames(clean), kept), !any(removed %in% colnames(clean)))
  clean@misc$multi_region_removal <- list(
    status = "subset_only_recompute_required", reduction = reduction,
    removed_cells = removed, checked_regions = checked,
    reason = "user_selected_region", timestamp = as.character(Sys.time())
  )
  # Existing reductions/cluster labels are pre-removal context until recomputed.
  region_summary <- data.frame(
    region_id = names(regions),
    checked = names(regions) %in% checked,
    n_cells = vapply(regions, function(r) length(r$barcodes), integer(1)),
    stringsAsFactors = FALSE
  )
  membership <- do.call(rbind, lapply(names(regions), function(id) {
    data.frame(region_id = id, barcode = regions[[id]]$barcodes,
               checked = id %in% checked, stringsAsFactors = FALSE)
  }))
  decisions <- data.frame(
    barcode = all_cells,
    decision = ifelse(all_cells %in% removed, "remove", "retain"),
    reason = ifelse(all_cells %in% removed, "user_selected_region", "outside_checked_regions"),
    stringsAsFactors = FALSE
  )
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
    stop("Cannot create output directory: ", out_dir)
  }
  export_dir <- tempfile("manual_regions_", tmpdir = out_dir)
  if (!dir.create(export_dir)) stop("Cannot create a fresh selection directory")
  utils::write.csv(region_summary, file.path(export_dir, "regions.csv"), row.names = FALSE)
  utils::write.csv(membership, file.path(export_dir, "region_membership.csv"), row.names = FALSE)
  utils::write.csv(decisions, file.path(export_dir, "cell_decisions.csv"), row.names = FALSE)
  utils::write.csv(object[[]][removed, , drop = FALSE],
                   file.path(export_dir, "removed_cell_metadata.csv"), row.names = TRUE)
  writeLines(removed, file.path(export_dir, "removed_barcodes.txt"))
  saveRDS(regions, file.path(export_dir, "regions_with_geometry.rds"))
  utils::write.csv(data.frame(before = length(all_cells), removed = length(removed),
                             remaining = length(kept)),
                   file.path(export_dir, "cell_counts.csv"), row.names = FALSE)
  clean_path <- file.path(export_dir, "after_manual_region_removal_subset_only.rds")
  saveRDS(clean, clean_path)
  writeLines(c("Manual regions selected by the user; no automatic doublet/QC diagnosis.",
               "Overlapping cells are removed once. Only checked regions are applied.",
               "This object is subset-only. Its reductions and cluster labels are pre-removal context.",
               "Recompute the applicable analysis workflow before final interpretation."),
             file.path(export_dir, "selection_status.txt"))
  writeLines(capture.output(utils::sessionInfo()), file.path(export_dir, "sessionInfo.txt"))
  list(object = clean, removed_cells = removed, kept_cells = kept,
       regions = regions, checked_regions = checked, export_dir = export_dir,
       subset_rds = clean_path, recompute_required = TRUE)
}

build_multi_region_cell_app <- function(object, out_dir, reduction = "umap",
                                        group.by = NULL, on_done = shiny::stopApp) {
  packages <- c("SeuratObject", "shiny", "plotly", "htmlwidgets", "jsonlite")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  cells <- colnames(object)
  if (!length(cells) || anyNA(cells) || anyDuplicated(cells)) stop("Invalid object cell names")
  embedding <- SeuratObject::Embeddings(object, reduction = reduction)
  if (ncol(embedding) < 2L || anyDuplicated(rownames(embedding)) ||
      !all(cells %in% rownames(embedding))) stop("Missing or invalid embedding cell names")
  embedding <- embedding[cells, 1:2, drop = FALSE]
  if (any(!is.finite(embedding))) stop("Embedding contains non-finite coordinates")
  metadata <- object[[]]
  if (!is.null(group.by) && !group.by %in% names(metadata)) stop("Unknown group.by column")
  groups <- if (is.null(group.by)) as.character(SeuratObject::Idents(object)[cells]) else
    as.character(metadata[cells, group.by])
  groups[is.na(groups)] <- "unknown"
  palette <- stats::setNames(grDevices::hcl.colors(length(unique(groups)), "Dark 3"), unique(groups))
  base_colors <- unname(palette[groups])
  frame <- data.frame(barcode = cells, x = embedding[, 1], y = embedding[, 2],
                      group = groups, stringsAsFactors = FALSE)
  padded_range <- function(values) {
    limits <- range(values)
    padding <- if (diff(limits) > 0) diff(limits) * 0.07 else max(abs(limits[1]) * 0.07, 1)
    limits + c(-padding, padding)
  }
  ui <- shiny::fluidPage(
    shiny::titlePanel("UMAP multi-region cell removal"),
    shiny::tags$head(shiny::tags$style(shiny::HTML(
      "#umap_regions .shapelayer path { pointer-events: none; }"
    )), shiny::tags$script(shiny::HTML("
      $(document).on('change.multiRegions', '#checked_regions input[type=checkbox]', function(e) {
        // Only actual user clicks; ignore Shiny's programmatic change echoes.
        if (!e.originalEvent) return;
        var boxes = Array.from(document.querySelectorAll('#checked_regions input[type=checkbox]'));
        Shiny.setInputValue('region_checks', {
          visible: boxes.map(function(b) { return b.value; }),
          checked: boxes.filter(function(b) { return b.checked; }).map(function(b) { return b.value; })
        }, {priority:'event'});
      });
      Shiny.addCustomMessageHandler('multi_region_ack', function(message) {
        var el = document.getElementById('umap_regions');
        if (el && el.multiRegionAcknowledge) el.multiRegionAcknowledge(message.ids);
      });
      Shiny.addCustomMessageHandler('toggle_delete_checked', function(message) {
        var button = document.getElementById('delete_checked');
        if (!button) return;
        var disabled = !!message.disabled;
        button.disabled = disabled;
        if (disabled) button.classList.add('disabled');
        else button.classList.remove('disabled');
      });
    "))),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 3,
        shiny::strong("Accumulation mode: each completed draw adds one region"),
        shiny::p("Release the mouse to save a region. Continue drawing directly; no Shift key or per-region confirmation is needed."),
        shiny::radioButtons("draw_mode", "Mouse tool",
          choices = c("Lasso select" = "lasso", "Box select" = "select",
                      "Pan" = "pan", "Zoom" = "zoom"), selected = "lasso", inline = TRUE),
        shiny::p("Check regions to remove. Orange cells are the final deletion set; overlaps count once."),
        shiny::checkboxGroupInput("checked_regions", "Saved regions", choices = character()),
        shiny::actionButton("select_all", "Select all"),
        shiny::actionButton("select_none", "Select none"), shiny::hr(),
        shiny::actionButton("undo_region", "Undo last region"),
        shiny::actionButton("clear_regions", "Clear regions"), shiny::hr(),
        shiny::strong(shiny::textOutput("selection_count")),
        shiny::actionButton("delete_checked", "Remove cells in checked regions",
                            class = "btn-danger", disabled = TRUE),
        shiny::actionButton("cancel", "Cancel and exit"),
        shiny::p("Removal is saved as a new object; the original object and file are preserved."),
        shiny::p(paste0("Selector version: ", MULTI_REGION_SELECTOR_VERSION)),
        shiny::div(style = "color:#b42318;", shiny::textOutput("status"))
      ),
      shiny::mainPanel(width = 9,
        shiny::p("Multiple numbered outlines stay visible. Checked regions use solid purple; unchecked regions use dotted gray."),
        plotly::plotlyOutput("umap_regions", height = "720px")
      )
    )
  )
  server <- function(input, output, session) {
    regions <- shiny::reactiveVal(list())
    # Authoritative server state: never derive accumulated selections from the
    # asynchronous echo of updateCheckboxGroupInput().
    checked_ids <- shiny::reactiveVal(character())
    serial <- 0L
    seen_events <- character()
    finished <- FALSE
    status <- shiny::reactiveVal("")
    selected <- shiny::reactive({
      multi_region_cells(regions(), checked_ids(), cells)
    })
    sync_regions <- function(value, checked) {
      regions(value)
      checked_ids(intersect(as.character(checked), names(value)))
      region_ids <- as.character(names(value))
      labels <- vapply(region_ids, function(id) {
        paste0(id, " (", length(value[[id]]$barcodes), " cells)")
      }, character(1))
      shiny::updateCheckboxGroupInput(session, "checked_regions",
        choices = stats::setNames(region_ids, labels), selected = checked_ids())
    }
    shiny::observeEvent(input$region_checks, {
      change <- input$region_checks
      visible <- as.character(unlist(change$visible, use.names = FALSE))
      checked <- as.character(unlist(change$checked, use.names = FALSE))
      # Keep new regions that were not yet visible when the user clicked.
      ids <- names(regions())
      checked_ids(union(setdiff(checked_ids(), visible), intersect(checked, ids)))
    }, ignoreInit = TRUE, priority = 100)
    output$umap_regions <- plotly::renderPlotly({
      # A single trace with barcode customdata avoids trace-index/cell-index confusion.
      p <- plotly::plot_ly(frame, x = ~x, y = ~y, customdata = ~barcode,
        text = ~paste0(barcode, "<br>", group), hoverinfo = "text",
        type = "scattergl", mode = "markers", source = "multi_region",
        marker = list(size = 4, color = base_colors, opacity = 0.8),
        selected = list(marker = list(opacity = 1)),
        unselected = list(marker = list(opacity = 0.8)), showlegend = FALSE)
      p <- plotly::layout(p, dragmode = "lasso", uirevision = "fixed-original-embedding",
        xaxis = list(title = colnames(embedding)[1], zeroline = FALSE,
                     range = padded_range(frame$x), autorange = FALSE),
        yaxis = list(title = colnames(embedding)[2], zeroline = FALSE,
                     range = padded_range(frame$y), autorange = FALSE),
        margin = list(l = 65, r = 20, b = 55, t = 20))
      p <- plotly::config(p, displaylogo = FALSE,
                          displayModeBar = TRUE,
                          modeBarButtonsToAdd = c("select2d", "lasso2d"))
      htmlwidgets::onRender(p, "function(el, x) {
        if (el.multiRegionSelectedHandler) el.removeListener('plotly_selected', el.multiRegionSelectedHandler);
        var pending = [], serial = 0;
        var epoch = Date.now().toString(36) + '-' + Math.random().toString(36).slice(2);
        el.multiRegionAcknowledge = function(ids) {
          ids = Array.isArray(ids) ? ids : [ids];
          pending = pending.filter(function(item) { return ids.indexOf(item.event_id) < 0; });
        };
        el.multiRegionSelectedHandler = function(e) {
          if (!e || !e.points || !e.points.length) return;
          var barcodes = e.points.map(function(p) { return p.customdata; });
          var vertices = null;
          if (e.lassoPoints) vertices = e.lassoPoints;
          else if (e.range) {
            var a = e.range.x, b = e.range.y;
            vertices = {x:[a[0],a[1],a[1],a[0]], y:[b[0],b[0],b[1],b[1]]};
          }
          // Send every unacknowledged stroke. Shiny can coalesce fast updates,
          // but the latest batch still contains every region and IDs deduplicate retries.
          pending.push({event_id:epoch + '-' + (++serial), barcodes:barcodes, vertices:vertices});
          Shiny.setInputValue('region_events', JSON.stringify(pending), {priority:'event'});
          // Stored shapes are separate from Plotly's replaceable active selection.
          // Clear the latter immediately so a new drag can start inside an old area.
          Plotly.relayout(el, {selections:[]});
        };
        el.on('plotly_selected', el.multiRegionSelectedHandler);
        Shiny.setInputValue('plot_ready', epoch, {priority:'event'});
      }")
    })
    add_region <- function(event) {
      picked <- unique(as.character(unlist(event$barcodes, use.names = FALSE)))
      if (!length(picked)) return()
      if (anyNA(picked) || !all(picked %in% cells)) {
        status("Selection contained unknown cells and was not added. Please draw again.")
        return()
      }
      vertices <- event$vertices
      if (!is.null(vertices)) {
        vx <- as.numeric(unlist(vertices$x)); vy <- as.numeric(unlist(vertices$y))
        if (length(vx) != length(vy) || length(vx) < 3L ||
            any(!is.finite(c(vx, vy)))) vertices <- NULL else
          vertices <- data.frame(x = vx, y = vy)
      }
      current <- regions()
      serial <<- serial + 1L
      id <- sprintf("R%03d", serial)
      current[[id]] <- list(barcodes = cells[cells %in% picked], vertices = vertices,
                            created = as.character(Sys.time()))
      sync_regions(current, union(checked_ids(), id))
      status("")
    }
    shiny::observeEvent(input$region_events, {
      # Decode explicitly: Shiny may simplify an array of event objects into a
      # data.frame, which loses the intended list-of-strokes structure.
      events <- tryCatch(jsonlite::fromJSON(input$region_events, simplifyVector = FALSE),
                         error = function(e) NULL)
      if (!is.list(events)) {
        status("Selection data could not be read. Please draw again.")
        return()
      }
      ack <- character()
      for (event in events) {
        if (!is.list(event)) next
        id <- as.character(event$event_id)
        if (length(id) != 1L || is.na(id) || !nzchar(id)) next
        if (!id %in% seen_events) {
          add_region(event)
          seen_events <<- c(seen_events, id)
        }
        ack <- c(ack, id)
      }
      session$sendCustomMessage("multi_region_ack", list(ids = as.list(ack)))
    }, ignoreInit = TRUE)
    shiny::observeEvent(list(input$draw_mode, input$plot_ready), {
      shiny::req(input$plot_ready)
      mode <- input$draw_mode
      if (!mode %in% c("lasso", "select", "pan", "zoom")) return()
      plotly::plotlyProxyInvoke(plotly::plotlyProxy("umap_regions", session),
                                "relayout", list(dragmode = mode, selections = list()))
    })
    shiny::observeEvent(input$select_all, {
      sync_regions(regions(), names(regions()))
    })
    shiny::observeEvent(input$select_none, {
      sync_regions(regions(), character())
    })
    shiny::observeEvent(input$undo_region, {
      current <- regions()
      if (!length(current)) return()
      current <- head(current, -1L)
      sync_regions(current, intersect(checked_ids(), names(current)))
    })
    shiny::observeEvent(input$clear_regions, { sync_regions(list(), character()) })
    shiny::observe({
      shiny::req(input$plot_ready)
      current <- regions()
      checked <- checked_ids()
      chosen <- selected()
      shapes <- lapply(names(current), function(id) {
        vertices <- current[[id]]$vertices
        if (is.null(vertices)) return(NULL)
        path <- paste0("M ", paste(vertices$x, vertices$y, sep = ",", collapse = " L "), " Z")
        list(type = "path", path = path, xref = "x", yref = "y", layer = "above",
             line = list(color = if (id %in% checked) "#7C3AED" else "#9CA3AF",
                         width = 2, dash = if (id %in% checked) "solid" else "dot"),
             fillcolor = "rgba(0,0,0,0)")
      })
      shapes <- Filter(Negate(is.null), shapes)
      annotations <- Filter(Negate(is.null), lapply(names(current), function(id) {
        vertices <- current[[id]]$vertices
        if (is.null(vertices)) return(NULL)
        list(x = min(vertices$x), y = max(vertices$y), xref = "x", yref = "y",
             text = id, showarrow = FALSE, xanchor = "left", yanchor = "bottom",
             font = list(color = if (id %in% checked) "#7C3AED" else "#6B7280"),
             bgcolor = "rgba(255,255,255,0.85)")
      }))
      colors <- base_colors
      colors[cells %in% chosen] <- "#F26D4F"
      proxy <- plotly::plotlyProxy("umap_regions", session)
      plotly::plotlyProxyInvoke(proxy, "restyle",
        list("marker.color" = list(colors), selectedpoints = list(NULL)), list(0L))
      # Remove Plotly's editable transient selection so the next drag can start
      # inside a previous region; our stored outlines remain in shapes.
      plotly::plotlyProxyInvoke(proxy, "relayout",
                                list(shapes = shapes, annotations = annotations, selections = list()))
    })
    output$selection_count <- shiny::renderText({
      paste0("Saved ", length(regions()), " regions; pending removal: ", length(selected()),
             " / ", length(cells), " cells; retained: ", length(cells) - length(selected()), ".")
    })
    shiny::observe({
      disabled <- !length(selected()) || length(selected()) == length(cells) || finished
      session$sendCustomMessage("toggle_delete_checked", list(disabled = disabled))
    })
    output$status <- shiny::renderText(status())
    shiny::observeEvent(input$delete_checked, {
      if (finished) return()
      status("Submitting removal request; please wait.")
      session$sendCustomMessage("toggle_delete_checked", list(disabled = TRUE))
      result <- tryCatch(multi_region_commit(object, regions(),
        checked_ids(), out_dir, reduction),
        error = function(e) {
          status(conditionMessage(e))
          session$sendCustomMessage("toggle_delete_checked", list(disabled = FALSE))
          NULL
        })
      if (!is.null(result)) {
        finished <<- TRUE
        status("Removal request received; returning to the R session.")
        on_done(result)
      }
    })
    shiny::observeEvent(input$cancel, {
      if (!finished) { finished <<- TRUE; on_done(NULL) }
    })
    session$onSessionEnded(function() {
      if (!finished) { finished <<- TRUE; on_done(NULL) }
    })
  }
  shiny::shinyApp(ui, server)
}

multi_region_cell_selector <- function(object, out_dir, reduction = "umap",
                                        group.by = NULL, launch.browser = interactive(),
                                        host = "127.0.0.1", port = NULL) {
  app <- build_multi_region_cell_app(object, out_dir, reduction, group.by)
  shiny::runApp(app, launch.browser = launch.browser, host = host, port = port)
}
