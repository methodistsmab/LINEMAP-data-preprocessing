#!/usr/bin/env Rscript

# Combine one-mouse save_list RDS files into one analysis-ready save_list object.
# This does not require raw KPTracer data.

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(script_file))

out_dir <- Sys.getenv("LINEMAP_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
pattern <- Sys.getenv("LINEMAP_MOUSE_PATTERN", unset = "save_list_.*[.]rds$")
combined_file <- Sys.getenv("LINEMAP_COMBINED_SAVE_LIST", unset = file.path(out_dir, "save_list_all_mice_linemap.rds"))

mouse_files <- list.files(out_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
mouse_files <- mouse_files[!basename(mouse_files) %in% basename(combined_file)]
if (length(mouse_files) == 0L) {
  stop("No mouse save_list files found under: ", out_dir)
}

save_list <- list()
for (file in mouse_files) {
  x <- readRDS(file)
  if (!is.list(x) || is.null(names(x))) {
    warning("Skipping unnamed list file: ", file)
    next
  }
  save_list[names(x)] <- x
}

save_list <- save_list[order(names(save_list))]
saveRDS(save_list, combined_file)

message("Combined ", length(save_list), " mouse entries.")
message("Wrote: ", combined_file)
