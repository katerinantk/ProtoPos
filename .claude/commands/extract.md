Extract Pol II active site positions from aligned BAMs using Protopos.py.

User request: $ARGUMENTS

Activate conda first: source ~/miniconda3/bin/activate protopos

Default invocation:
```bash
source ~/miniconda3/bin/activate protopos && python /mnt/d/Lab/Botos_Transcription_Simulation_Project/ProtoPos/Protopos.py
```

Key parameters:
- protocol: "pro_seq" (5' end), "gro_seq" (3' end), or "mnet_seq" (3' end)
- Input: Deduplicated BAMs in alignments/
- Output: *.3prime.sorted.bam or *.5prime.sorted.bam with 1 bp reads

After running:
1. Report output files and their sizes
2. Verify the correct end was extracted (5' for PRO-seq, 3' for GRO/mNET-seq)
3. Run samtools idxstats to check chromosome distribution
