#!/usr/bin/env Rscript

# Build one mouse-level save_list entry from stored LINEMAP transition inputs.
# Set MOUSE_ID to the mouse ID to rebuild.

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "mit_linemap_helpers.R"))

mouse_id <- Sys.getenv("MOUSE_ID", unset = "3522")
kptracer_dir <- Sys.getenv("KPTRACER_DATA_DIR", unset = default_kptracer_dir)
out_dir <- Sys.getenv("LINEMAP_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
transition_file <- Sys.getenv(
  "LINEMAP_TRANSITION_INPUTS",
  unset = file.path(out_dir, "linemap_transition_inputs_fixed_all.rds")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading KPTracer data from: ", kptracer_dir)
tables <- load_kptracer_tables(kptracer_dir)

message("Loading transition inputs from: ", transition_file)
transition_inputs <- readRDS(transition_file)

build_and_save_one_mouse <- function(mouse_id, tables, transition_inputs, out_dir) {
  message("Building mouse entry: ", mouse_id)
  entry <- build_mouse_save_entry(mouse_id, tables, transition_inputs)
  save_list <- setNames(list(entry), mouse_id)

  mouse_dir <- file.path(out_dir, paste0("mouse_", mouse_id))
  dir.create(mouse_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(save_list, file.path(mouse_dir, paste0("save_list_", mouse_id, ".rds")))
  saveRDS(entry$sample, file.path(mouse_dir, paste0("sample_", mouse_id, ".rds")))
  saveRDS(entry$sequence, file.path(mouse_dir, paste0("sequence_", mouse_id, ".rds")))
  saveRDS(entry$D, file.path(mouse_dir, paste0("distance_with_C0_", mouse_id, ".rds")))
  write.csv(entry$table, file.path(mouse_dir, paste0("cell_table_", mouse_id, ".csv")), row.names = FALSE)

  message("Wrote mouse outputs to: ", mouse_dir)
  message("n_cells: ", entry$meta$n_cells)
  message("n_intBC: ", entry$meta$n_intBC)

  data.frame(
    mouse_id = mouse_id,
    n_cells = entry$meta$n_cells,
    n_intBC = entry$meta$n_intBC,
    out_dir = mouse_dir,
    stringsAsFactors = FALSE
  )
}

if (toupper(mouse_id) == "ALL") {
  stop("This script rebuilds one mouse at a time. Set MOUSE_ID to a single mouse ID, for example MOUSE_ID=3522.")
}

build_and_save_one_mouse(mouse_id, tables, transition_inputs, out_dir)
