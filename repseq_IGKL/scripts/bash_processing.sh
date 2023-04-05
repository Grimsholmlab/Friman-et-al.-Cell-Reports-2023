## Process fastq files including UMI extraction and full length IGKL repseq assembly (documentation: https://migec.readthedocs.io/en/latest/index.html and https://mixcr.com/)

migec CheckoutBatch -cute barcodes_file.txt MIG_UMI_OUTPUT/

migec Histogram MIG_UMI_OUTPUT/ histogram/

migec AssembleBatch --force-collision-filter --force-overseq 5 MIG_UMI_OUTPUT/ histogram/ assemble/

mixcr align -p kAligner2 -s hsa -OreadsLayout=Collinear -OvParameters.geneFeatureToAlign=VTranscript --report alignmentReport.txt files_R1.t5.cf.fastq files_R2.t5.cf.fastq alignments.vdjca

mixcr assemble -OassemblingFeatures=VDJRegion -OseparateByC=true -ObadQualityThreshold=15 -OqualityAggregationType=Average --report assembleReport.txt alignments.vdjca clones.clns

mixcr exportClones -c IGL -o -t clones.clns clones.txt

## Process output from IMGT/HighV-QUEST web tool using ChangeO (documentation: https://changeo.readthedocs.io/en/stable/)

MakeDb.py imgt -i files.txz -s files.fasta --extended

# Cluster sequences into clonal groups, determine --dist using distToNearest() and findThreshold() (documentation: https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/)
DefineClones.py -d files_db-pass.tsv --act set --model ham --norm len --dist

CreateGermlines.py -d files-pass_clone-pass.tsv -g dmask --cloned -r IMGT_Human_IGHV.fasta IMGT_Human_IGHD.fasta IMGT_Human_IGHJ.fasta

# Build phylogenetic trees using IgPhyML (documentation: https://igphyml.readthedocs.io/en/stable/index.html)
BuildTrees.py -d files_clone-pass_germ-pass.tsv --outname ex --log ex.log --ncdr3 --collapse --igphyml --clean all --nproc 1