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
