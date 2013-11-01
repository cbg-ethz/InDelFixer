# InDelFixer
A sensitive aligner for 454, Illumina and PacBio data, employing a full Smith-Waterman alignment against a reference.  

This java command line application aligns Next-Generation Sequencing (NGS) and third-generation reads
to a set of reference sequences, by a prior fast k-mer matching and removes indels, causing
frame shifts. In addition, only a specific region can be considered. 

An iterative refinement of the alignment can be performed, by alignment against the consensus sequence with wobbles.

The output is in SAM format.

## Features
 - Fully multithreaded
 - Performes a full Smith-Waterman alignment
 - Multiple sets of affine gap costs can be used to find optimal alignment for each read
 - Paired-end reads are properly paired with SAM Flags
 - Accepts multiple reference genomes with wobbles
 - Iterative alignment against a consensus with wobbles to increase alignment quality

#### PREREQUISITES TO RUN:
 - JDK 7 (http://jdk7.java.net/)
 - Get latest version: https://github.com/armintoepfer/InDelFixer/releases

## RUN:
#### 454/Roche:
`java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta`

But I encourage to convert the sff to fastq with `sff2fastq input.sff -o input.fastq`
<b>sff2fastq</b> can be installed with:
```
git clone git://github.com/indraniel/sff2fastq.git;
cd sff2fastq;
make;
```
 
#### Fasta / PacBio ccs:
`java -jar InDelFixer.jar -i libCase102.fasta -g referenceGenomes.fasta`

For PacBio input, please use `-noHashing` since the PacBio error rate is too high for a reliable kmer-matching.
 
#### Illumina paired end:
`java -jar InDelFixer.jar -i libCase102_R1.fastq -ir libCase102_R2.fastq -g referenceGenomes.fasta`

### High quality alignment
With parameter `-sensitive`, multiple affine gap costs are tested for each read and the best alignment is kept.

### Affine GAP costs
Gap costs for the used Smith-Waterman can be set with
```
-gop 3 (gap open)
-gex 1 (gap extend)
```
Predefined: 10 open & 3 extend. Tested with with PacBio, Illumina and 454 data on HIV, HCV and HBV data.

### Iterative refinement
The alignment can be improved by aligning against the consensus sequence. The parameter `-refine INT` takes a positive number as input and activates the iterative refinement. Only works if the alignment is against _one_ reference genome.

#### Remove conserved deletions:
During the iterative alignment, conserved deletions can be removed with `-rmDel`.

### Remove frame-shift causing deletions
With parameter `-fix`, frame-shift causing deletions are replaced with the consensus sequence.

### Line breaks
In the case that a single fastq entry is longer than four lines, which is caused by line breaks in the sequence and quality string, use `-flat`.

### Extract region:
In addition, only a specific region can be extracted with `-r begin-end`, for example a certain gene:
  `java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta -r 342-944`
  
### FILTER 
```
  -l      INT    : Minimal read-length prior alignment (default 0)
  -la     INT    : Minimal read-length after alignment (default 0)
  -ins    DOUBLE : The maximum percentage of insertions allowed [range 0.0 - 1.0] (default 1.0)
  -del    DOUBLE : The maximum percentage of deletions allowed [range 0.0 - 1.0] (default 1.0)
  -sub    DOUBLE : The maximum percentage of substitutions allowed [range 0.0 - 1.0] (default 1.0)
  -maxDel INT    : The maximum number of consecutive deletions allowed (default no filtering)
```

## Help:
Further help can be shown by running without additional parameters:
    `java -jar InDelFixer.jar`

### BAM output:
In order to convert the `reads.sam` into the BAM format, please install samtools and run:

    samtools view -bS reads.sam > out.bam; 
    samtools sort out.bam reads; 
    samtools index reads.bam; 
    rm out.bam;

### COMPILE (only for dev):
Install [Maven 3](http://maven.apache.org/)

    cd InDelFixer
    mvn -DartifactId=samtools -DgroupId=net.sf -Dversion=1.9.6 -Dpackaging=jar -Dfile=src/main/resources/jars/sam-1.96.jar -DgeneratePom=false install:install-file
    mvn clean package
    java -jar target/InDelFixer.jar

## CONTACT:
    Armin Töpfer
    armin.toepfer (at) gmail.com
    http://www.bsse.ethz.ch/cbg/people/armintoepfer
    
## Contributions
 [Armin Töpfer](http://www.bsse.ethz.ch/cbg/people/armintoepfer)  
 [David Seifert](http://www.bsse.ethz.ch/cbg/people/dseifert)  
 [Alexander Artyomenko](http://alan.cs.gsu.edu/NGS/?q=content/alexander-artyomenko)

## LICENSE:
 GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0
