# AnnotationLiftOver
A pipeline to liftover the gene structure of reference genome to other lines/accessions\

## install
go to the source folder:\
`make`\
this command will generate a file named `AnnotationLiftOver`\

## run it
Program `AnnotationLiftOver`\
Usage:  `AnnotationLiftOver` <command> [options]\
Commands:\
`getGenomeSequence`              create fasta file using reference genome sequence and variants calling records\
`coordinateLiftOver`             get the corresponding coordinate or reference coordinate at another accession\
`gffCoordinateLiftOver`          transform the reference GFF/GTF coordinate to the coordinate of another accession,
                               don't care about the completeness of start/stop codon, splice sites\
`getSequences`                   get the protein/CDS/gene sequence with genome sequence file and GFF/GTF file\
`annotationAndExonerateAndNovo`  transform the reference GFF/GTF coordinate to the coordinate of another accession,
                               complementing with trying to keep the ORF/splice sites and complete as possible by
                               genome sequence alignment. And then complementing with aligning CDS sequence of
                               reference to the genome sequence of target line and then complementing with
                               aligning protein sequence of reference to the genome sequence of target line
                               Finally the missing gene annotation would be complemented with
                               other annotation inputted\
                               
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
   