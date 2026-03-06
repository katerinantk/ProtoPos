Run the ProtoPos PRO-seq processing pipeline.

User request: $ARGUMENTS

The full pipeline (proseqanalysis.sh) runs:
1. fasterq-dump (SRA download)
2. fastp (adapter trimming, quality filtering)
3. bowtie2 alignment to hg38 (splice-unaware, appropriate for nascent RNA)
4. umi_tools dedup (UMI-based deduplication)
5. Protopos.py extraction (1 bp reads at Pol II active site)

Activate conda first: source ~/miniconda3/bin/activate protopos

Key biological correctness rules:
- PRO-seq: 5' end of read = Pol II active site
- GRO-seq: 3' end of read = Pol II active site
- mNET-seq: 3' end of read = Pol II active site
- Coordinate system: BAM is 0-based, GTF is 1-based — explicit conversion required
- Always use bowtie2 (splice-unaware), never STAR, for nascent RNA

After running:
1. Report output BAM files and their locations
2. Run basic QC: samtools flagstat on output BAMs
3. Verify strand-specific signal around known TSSs
