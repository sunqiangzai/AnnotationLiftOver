# AnnotationLiftOver

The gene sequence and gene structure variation between different accessions/line are very important for natural variation research.
In general the genomic variants between different lines/accession and reference genome is accessible.
And the gene structure (always in GFF/GTF format) of reference is well established.
While the gene structure of variant accessions/lines is also very interesting.\
This pipeline tries to lift the reference gene structure to variant accessions/lines.
The [alternative alignment problem](https://www.ncbi.nlm.nih.gov/pubmed/25701572) could lead to false positive splice sites disturb or ORF-shift predication.
Here we solved this problem by a dynamic programming algorithm.\
And I believe this pipeline could help to quantify the gene expression for non-reference accession/line and detect the difference expression level across different accession/line.
## install
go to the source folder:\
`make`\
this command will generate a file named `AnnotationLiftOver`

## run it
Program `AnnotationLiftOver`\
Usage:  `AnnotationLiftOver` <command> [options]\
Commands:\
`getGenomeSequence`              create fasta file using reference genome sequence and variant calling records\
`coordinateLiftOver`             get the corresponding coordinate of reference coordinate at another accession\
`gffCoordinateLiftOver`          transform the reference GFF/GTF coordinate to the coordinate of another accession,
                               don't care about the completeness of start/stop codon, splice sites\
`getSequences`                   get the protein/CDS/gene sequence with genome sequence file and GFF/GTF file\
`annotationAndExonerateAndNovo`  transform the reference GFF/GTF coordinate to the coordinate of another accession,
                               complementing with trying to keep the ORF/splice sites and complete as possible by
                               genome sequence alignment. And then complementing with aligning CDS sequence of
                               reference to the genome sequence of target line and then complementing with
                               aligning protein sequence of reference to the genome sequence of target line
                               Finally the missing gene annotation would be complemented with
                               other annotation inputted
                               
### getGenomeSequence
Usage:    * getGenomeSequence -r reference -v variants -o output\
Options\
   -h        produce help message\
   -r        reference genome in fasta format\
   -v        variant calling result in vcf/sdi format\
   -prefix   prefix for vcf records\
   -o        the output pseudo genome sequence in fasta format\

** -prefix is the prefix of chromosome name for vcf/sdi variant records. Like the chromosome in TAIR10 reference genome is Chr1, Chr2, Chr3, Chr4 and Chr5. While the chromosomes in vcf files from the 1001 genomes project were indicated with 1, 2, 3, 4 and 5.
So `-prefix Chr` should be set to make the software work properly. If this parameter is not set correctly, the software would act as no variant records in the input vcf/sdi file.   

### Contact
Bug report? Any question? Any suggestion? Any requirement?\
Please feel free to send E-mail to songbaoxing168@163.com.

### acknowledgement
I thank:\
 Hequan Sun (MPIPZ) for the discussion of algorithm design for gene structure alignment\
 Lukas Baumgarten (MPIPZ) for bug reporting

