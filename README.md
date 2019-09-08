## VCF projection from chromosome to transcriptome

The program uses ncls, a fast interval tree search program 
that considerably speeds up the ''stabbing'' problem of 
finding out the variant that overlaps the transcriptome interval. 

```

GTF_FILE = 'Mus_musculus.GRCm38.91.gtf'
VCF_FILE = 'mgp_v4_indels_NOD_PWK.vcf'


from vcfproj import vcfproj
vcf_gtf = vcfproj.projection(GTF_FILE, VCF_FILE)

```

`vcf_gtf` is a dataframe, that would contain the following columns 
```
[

'chrom_x', - chromosome of variation 
'start', - start of transcript
'end', - end of transcript
'gene', - gene name
'txome', - transcript name
'vcf_index', - ignore
'relative_pos', - relative position of variation on transcript
'transcript_length', - length of transcript
'chrom_y', - same as chrom_x
'pos', - position of variation
'id', - ID
'ref', - REF
'alt', - ALT
'qual', - QUAL
'filter', - FILTER
'info', - INFO
'format', - FORMAT
'samples' - SAMPLES 

]
```
