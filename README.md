# ProtoPos
This repository contains a pipeline for processing PRO-seq nascent transcription data. The workflow downloads raw FASTQ files, performs quality control and trimming, aligns reads with STAR, removes rRNA/tRNA contamination, extracts strand-specific 3′ ends and generates genome coverage tracks. It also contains a python pipeline that computes the average PRO-seq signal around the TSS. 
