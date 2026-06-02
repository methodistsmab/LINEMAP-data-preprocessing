#!/usr/bin/env Rscript

# Build LINEMAP transition inputs for the KPTracer/MIT mouse data.
# The output stores the three sgRNA transition matrices and their distance
# dictionaries for the downstream mouse-level reconstruction script.

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(script_file))
source(file.path(script_dir, "mit_linemap_helpers.R"))

kptracer_dir <- Sys.getenv("KPTRACER_DATA_DIR", unset = default_kptracer_dir)
out_dir <- Sys.getenv("LINEMAP_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
division <- as.numeric(Sys.getenv("LINEMAP_DIVISION", unset = "20"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading KPTracer data from: ", kptracer_dir)
tables <- load_kptracer_tables(kptracer_dir)
data.mutation <- tables$data.mutation

make_P <- function(values, readout_name) {
  values <- as.character(values)
  values <- values[!is.na(values)]
  states <- unique(c("NC", sort(setdiff(unique(values), "NC"))))
  P <- estimate_transition_matrix(
    sequence = matrix(values, ncol = 1),
    division = division,
    initial_state = "NC",
    true_states = states,
    missing_values = c("", "-", "MISSING"),
    include_missing_as_state = FALSE
  )
  rownames(P) <- colnames(P) <- states
  message(readout_name, ": ", length(states), " states")
  P
}

message("Calculating transition matrix for each guide RNA...")

P1 <- make_P(data.mutation$r1, "r1")
P2 <- make_P(data.mutation$r2, "r2")
P3 <- make_P(data.mutation$r3, "r3")

message("Building missing-aware LINEMAP distance dictionaries...")
tabs_base <- build_pair_loglik_tables(
  P_or_list = list(r1 = P1, r2 = P2, r3 = P3),
  division = division,
  idx0 = 1,
  missing_label = "MISSING",
  group = c("r1", "r2", "r3")
)

transition_inputs <- list(
  P = list(r1 = P1, r2 = P2, r3 = P3),
  dictionary = tabs_base
)

out_file <- file.path(out_dir, "linemap_transition_inputs_fixed_all.rds")
saveRDS(transition_inputs, out_file)
message("Wrote: ", out_file)
