# gorg-classifier-v2
Recruit, classify, and annotate reads against the Global Ocean Reference Genome dataset of SAGs

# Usage
Assuming you have:
* A local nextflow.config file (with a profile named 'charlie' in this example)
* Reads to classify (in this case, within the directory demo_reads/ )
  
```
nextflow run gorg_classifier_v2.1.nf \
--seqs './demo_reads/100ksub{1,2}.fq.gz' \
--db 'gorg_dark_v1' \
-profile 'charlie'
```

# db options
Currently gorg_dark_v1 is the only DB in V2, but additional versions will soon be released.
