library(STADyUM)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)

# ---- Inputs (edit these paths to match your setup) ----
args <- commandArgs(trailingOnly=TRUE)
plus_bw  <- if (length(args) >= 1) args[1] else "SRR28248970.plus.bw"
minus_bw <- if (length(args) >= 2) args[2] else "SRR28248970.minus.bw"
gtf      <- if (length(args) >= 3) args[3] else "data/gencode.v47.annotation.gtf"
sample_id <- if (length(args) >= 4) args[4] else "SRR28248970"

message("plus_bw:  ", plus_bw)
message("minus_bw: ", minus_bw)
message("gtf:      ", gtf)
message("sample:   ", sample_id)

# Remove overlapping ranges using a greedy keep-longest strategy.
# Strand-aware: only ranges on the same strand are considered overlapping,
# because STADyUM uses strand-specific BigWig files.
remove_overlaps_greedy <- function(gr) {
  # Sort by width descending so we preferentially keep longer genes.
  gr <- gr[order(width(gr), decreasing=TRUE)]
  # Strand-aware overlap detection (default ignore.strand=FALSE).
  hits <- findOverlaps(gr, gr)
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  if (length(hits) == 0) return(gr)
  # Build set of indices to drop: for each overlapping pair, drop the shorter
  # one (higher index, since we sorted by decreasing width).
  drop_idx <- unique(subjectHits(hits)[subjectHits(hits) > queryHits(hits)])
  gr[-drop_idx]
}

gtf_gr <- import(gtf)

# ---- choose ONE transcript per gene: the longest transcript ----
tx <- gtf_gr[gtf_gr$type == "transcript"]
tx <- GenomeInfoDb::keepStandardChromosomes(tx, pruning.mode="coarse")
tx <- tx[!is.na(strand(tx)) & strand(tx) %in% c("+","-")]

gid <- mcols(tx)$gene_id
if (is.null(gid)) stop("No gene_id found in transcript entries of the GTF.")

# For each gene_id, keep transcript with max width.
tx_split <- split(tx, gid)
tx_rep <- suppressWarnings(do.call(c, lapply(tx_split, function(x) x[which.max(width(x))])))

names(tx_rep) <- mcols(tx_rep)$gene_id  # name by gene_id for matching

# ---- pause regions (TSS +/- 250) from representative transcript ----
tss <- resize(tx_rep, width=1, fix="start")
pause_regions <- promoters(tss, upstream=250, downstream=250)
names(pause_regions) <- names(tx_rep)

# ---- gene body: filter long transcripts first, then trim ----
gene_body <- tx_rep[width(tx_rep) > 2000]  # prevent negative widths after trimming
start(gene_body) <- ifelse(strand(gene_body) == "+", start(gene_body) + 500, start(gene_body))
end(gene_body)   <- ifelse(strand(gene_body) == "-", end(gene_body) - 500, end(gene_body))
gene_body <- gene_body[width(gene_body) > 1000]

# ---- ensure pause/gene_body have same genes ----
common0 <- intersect(names(pause_regions), names(gene_body))
pause_regions <- pause_regions[common0]
gene_body     <- gene_body[common0]

# ---- remove overlapping ranges (strand-aware, greedy) ----
pause_regions <- remove_overlaps_greedy(pause_regions)
gene_body     <- remove_overlaps_greedy(gene_body)

common <- intersect(names(pause_regions), names(gene_body))
pause_regions <- pause_regions[common]
gene_body     <- gene_body[common]

message("Representative transcripts: ", length(tx_rep))
message("After length filters: ", length(common0))
message("After overlap removal: ", length(common))

output_csv <- paste0(sample_id, "_STADyUM_rates.csv")

fit <- STADyUM::estimateTranscriptionRates(
  plus_bw,
  minus_bw,
  pause_regions,
  gene_body,
  sample_id
)

res <- as.data.frame(STADyUM::rates(fit))
write.csv(res, output_csv, row.names=FALSE)
message("DONE: wrote ", output_csv, " (", nrow(res), " genes)")
