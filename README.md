# InDelFixer
An insertion and deletion fixing aligner for 454, Illumina and PacBio.
This java command line application aligns Next-Generation Sequencing (NGS) and third-generation reads
to a set of reference sequences, by a fast k-mer matching and removes indels, causing
frame shifts. In addition, only a specific region can be considered. 

The output is in SAM format.
- - -

#### PREREQUISITES TO RUN:
 - JDK 7 (http://jdk7.java.net/)
 - Get latest version: https://github.com/armintoepfer/InDelFixer/downloads

## RUN:
#### 454/Roche:
`java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta`
 
#### Fasta / PacBio ccs:
`java -jar InDelFixer.jar -i libCase102.fasta -g referenceGenomes.fasta -pacbio`
 
#### Illumina paired end:
`java -jar InDelFixer.jar -i libCase102_R1.fastq -ir libCase102_R2.fastq -g referenceGenomes.fasta -illumina`

### Affine GAP costs
Gap costs for the used Smith-Waterman can be set with
```
-gop 3 (gap open)
-gex 1 (gap extend)

or
-illumina (46 open / 10 extend)
-pacbio (10 open / 10 extend)
-roche (10 open/ 8 extend)
```

### Line breaks
In the case that a single fastq entry is longer than four lines, which is caused by line breaks in the sequence and quality string, use `-flat`.

### Remove InDels:
Remove insertions and deletions with `-adjust`

### Extract region:
In addition, only a specific region can be extracted with `-r begin-end`, for example a certain gene:
  `java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta -r 342-944`
* * *
### Help:
Further help can be showed by running without additional parameters:
    `java -jar InDelFixer.jar`

### BAM output:
In order to convert the `reads.sam` into the BAM format, please install samtools and run:

    samtools view -bS reads.sam > out.bam; 
    samtools sort out.bam reads; 
    samtools index reads.bam; 
    rm out.bam;

### COMPILE (only for dev):
Install Maven 3 (http://maven.apache.org/)

    cd InDelFixer
    mvn clean package
    java -jar InDelFixer/target/InDelFixer.jar
* * *
# CONTACT:
    Armin Töpfer
    armin.toepfer (at) gmail.com
    http://www.bsse.ethz.ch/cbg/people/toepfera

# LICENSE:
 GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0
