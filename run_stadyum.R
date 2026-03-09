suppressPackageStartupMessages({
  library(STADyUM)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# ---- Inputs (CLI args) ----
args <- commandArgs(trailingOnly = TRUE)
plus_bw   <- if (length(args) >= 1) args[1] else "SRR28248970.plus.bw"
minus_bw  <- if (length(args) >= 2) args[2] else "SRR28248970.minus.bw"
gtf       <- if (length(args) >= 3) args[3] else "data/gencode.v47.annotation.gtf"
sample_id <- if (length(args) >= 4) args[4] else "SRR28248970"
outdir    <- if (length(args) >= 5) args[5] else file.path("outputs", sample_id)
steric    <- if (length(args) >= 6) as.logical(args[6]) else FALSE

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("plus_bw:   ", plus_bw)
message("minus_bw:  ", minus_bw)
message("gtf:       ", gtf)
message("sample:    ", sample_id)
message("outdir:    ", outdir)
message("stericH:   ", steric)

# ---- Strand-aware greedy keep-longest (simple) overlap remover ----
remove_overlaps_greedy <- function(gr) {
  gr <- gr[order(width(gr), decreasing = TRUE)]
  hits <- findOverlaps(gr, gr, ignore.strand = FALSE)
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  if (length(hits) == 0) return(gr)
  drop_idx <- unique(subjectHits(hits)[subjectHits(hits) > queryHits(hits)])
  gr[-drop_idx]
}

# ---- Import GTF ----
gtf_gr <- import(gtf)

# ---- Pick ONE transcript per gene: the longest transcript ----
tx <- gtf_gr[gtf_gr$type == "transcript"]
tx <- GenomeInfoDb::keepStandardChromosomes(tx, pruning.mode = "coarse")
tx <- tx[!is.na(strand(tx)) & strand(tx) %in% c("+", "-")]

gid <- mcols(tx)$gene_id
if (is.null(gid)) stop("No gene_id found in transcript entries of the GTF.")

tx_split <- split(tx, gid)
tx_rep_list <- lapply(tx_split, function(x) x[which.max(width(x))])
tx_rep <- unlist(GRangesList(tx_rep_list), use.names = FALSE)
mcols(tx_rep)$gene_id <- names(tx_split)   # keep the gene_id attached
names(tx_rep) <- mcols(tx_rep)$gene_id

# ---- Pause regions: TSS +/- 250 ----
tss <- resize(tx_rep, width = 1, fix = "start")
pause_regions <- promoters(tss, upstream = 250, downstream = 250)
names(pause_regions) <- names(tx_rep)

# ---- Gene body: filter long first, then trim ----
gene_body <- tx_rep[width(tx_rep) > 2000]
start(gene_body) <- ifelse(strand(gene_body) == "+", start(gene_body) + 500, start(gene_body))
end(gene_body)   <- ifelse(strand(gene_body) == "-", end(gene_body) - 500, end(gene_body))
gene_body <- gene_body[width(gene_body) > 1000]

# ---- Ensure pause/gene_body share same genes ----
common0 <- intersect(names(pause_regions), names(gene_body))
pause_regions <- pause_regions[common0]
gene_body     <- gene_body[common0]

# ---- Remove overlaps (strand-aware) ----
pause_regions <- remove_overlaps_greedy(pause_regions)
gene_body     <- remove_overlaps_greedy(gene_body)

common <- intersect(names(pause_regions), names(gene_body))
pause_regions <- pause_regions[common]
gene_body     <- gene_body[common]

message("Representative transcripts: ", length(tx_rep))
message("After length filters:      ", length(common0))
message("After overlap removal:     ", length(common))

# ---- Fit model (EM) ----
# Note: steric hindrance mode may require omegaScale for your dataset.
# If steric=TRUE fails, rerun with steric=FALSE.
fit <- STADyUM::estimateTranscriptionRates(
  plus_bw,
  minus_bw,
  pause_regions,
  gene_body,
  name = sample_id,
  stericHindrance = steric
)

# ---- Save outputs ----
# ---- Save outputs ----
rates_obj <- STADyUM::rates(fit)

# Convert to base data.frame (may include list-columns)
rates_tbl <- as.data.frame(rates_obj)

csv_file <- file.path(outdir, paste0(sample_id, "_STADyUM_rates.csv"))
rds_fit  <- file.path(outdir, paste0(sample_id, "_STADyUM_fit.rds"))
rds_tbl  <- file.path(outdir, paste0(sample_id, "_STADyUM_rates.rds"))

# Drop list-columns for CSV (CSV cannot store lists)
is_listcol <- vapply(rates_tbl, is.list, logical(1))
if (any(is_listcol)) {
  message("Dropping list-columns for CSV: ",
          paste(names(rates_tbl)[is_listcol], collapse = ", "))
  rates_tbl_csv <- rates_tbl[, !is_listcol, drop = FALSE]
} else {
  rates_tbl_csv <- rates_tbl
}

write.csv(rates_tbl_csv, csv_file, row.names = FALSE)

# Keep full objects in RDS (safe for list-columns)
saveRDS(fit, rds_fit)
saveRDS(rates_obj, rds_tbl)

message("Wrote: ", csv_file, " (", nrow(rates_tbl_csv), " genes)")
message("Wrote: ", rds_fit)
message("Wrote: ", rds_tbl)

# ---- “Vignette-style” plots (written to files) ----
STADyUM::plotChiDistrib(
  fit,
  file = file.path(outdir, paste0(sample_id, "_ChiDistribution.png"))
)

STADyUM::plotExpectedVsActualPauseSiteCounts(
  fit,
  file = file.path(outdir, paste0(sample_id, "_ExpectedVsActualPauseCounts.png"))
)

STADyUM::plotMeanPauseDistrib(
  fit,
  file = file.path(outdir, paste0(sample_id, "_MeanPauseDistrib.png"))
)

message("DONE plots -> ", outdir)