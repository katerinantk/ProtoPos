library(STADyUM)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)

plus_bw  <- "/home/katerinan/protopos/SRR28248970.plus.bw"
minus_bw <- "/home/katerinan/protopos/SRR28248970.minus.bw"
gtf      <- "/home/katerinan/protopos/data/gencode.v47.annotation.gtf"

# Drop ANY overlapping ranges (minimal pain, makes STADyUM happy)
drop_any_overlaps <- function(gr) {
  repeat {
    hits <- findOverlaps(gr, gr, ignore.strand=TRUE)
    hits <- hits[queryHits(hits) != subjectHits(hits)]
    if (length(hits) == 0) break
    gr <- gr[-unique(queryHits(hits))]
  }
  gr
}

gtf_gr <- import(gtf)

# ---- choose ONE transcript per gene: the longest transcript ----
tx <- gtf_gr[gtf_gr$type == "transcript"]
tx <- GenomeInfoDb::keepStandardChromosomes(tx, pruning.mode="coarse")
tx <- tx[!is.na(strand(tx)) & strand(tx) %in% c("+","-")]

gid <- mcols(tx)$gene_id
if (is.null(gid)) stop("No gene_id found in transcript entries of the GTF.")

# For each gene_id, keep transcript with max width
tx_split <- split(tx, gid)
tx_rep <- suppressWarnings(do.call(c, lapply(tx_split, function(x) x[which.max(width(x))])))

names(tx_rep) <- mcols(tx_rep)$gene_id  # name by gene_id for matching

# ---- pause regions (TSS ± 250) from representative transcript ----
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

# ---- remove remaining overlaps (much fewer now) ----
pause_regions <- drop_any_overlaps(pause_regions)
gene_body     <- drop_any_overlaps(gene_body)

common <- intersect(names(pause_regions), names(gene_body))
pause_regions <- pause_regions[common]
gene_body     <- gene_body[common]

message("Representative transcripts: ", length(tx_rep))
message("After length filters: ", length(common0))
message("After overlap removal: ", length(common))

fit <- STADyUM::estimateTranscriptionRates(
  plus_bw,
  minus_bw,
  pause_regions,
  gene_body,
  "SRR28248970"
)

res <- as.data.frame(STADyUM::rates(fit))
write.csv(res, "SRR28248970_STADyUM_rates.csv", row.names=FALSE)
message("DONE: wrote SRR28248970_STADyUM_rates.csv (", nrow(res), " genes)")
