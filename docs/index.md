# About IsoQuant

IsoQuant is a tool for the genome-based analysis of long RNA reads, such as PacBio or
Oxford Nanopores. IsoQuant allows to reconstruct and quantify transcript models with
high precision and decent recall. If the reference annotation is given, IsoQuant also
assigns reads to the annotated isoforms based on their intron and exon structure.
IsoQuant further performs annotated gene, isoform, exon and intron quantification.
If reads are grouped (e.g. according to cell type), counts are reported according to the provided grouping.

IsoQuant consists of two stages, which generate its own output:

1. Reference-based analysis. Runs only if reference annotation is provided. Performs read-to-isoform assignment,
splice site correction and abundance quantification for reference genes/transcripts.
2. Transcript discovery. Reconstructs transcript models and performs abundance quantification for discovered isoforms.

Latest IsoQuant version can be downloaded from [https://github.com/ablab/IsoQuant/releases/latest](https://github.com/ablab/IsoQuant/releases/latest).

### IsoQuant pipeline
![Pipeline](isoquant_pipeline.png)

