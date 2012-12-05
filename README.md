# InDelFixer
A insertion and deletion fixing aligner for 454, Illumina and PacBio.
This java command line application aligns Next-Generation Sequencing (NGS) and third-generation reads
to a set of reference sequences, by a fast k-mer matching and removes indels, causing
frame shifts. In addition, only a specific region can be considered.
- - -

#### PREREQUISITES TO RUN:
 - JDK 7 (http://jdk7.java.net/)
 - Get latest version: https://github.com/armintoepfer/InDelFixer/downloads

## RUN:
#### 454/Roche:
`java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta`
 
#### Fasta / PacBio ccs:
`java -jar InDelFixer.jar -i libCase102.fasta -g referenceGenomes.fasta`
 
#### Illumina paired end:
  `java -jar InDelFixer.jar -i libCase102_R1.fastq -ir libCase102_R2.fastq -g referenceGenomes.fasta`

### Extract region:
 In addition only a specific region can be extracted, for example a certain gene:
  `java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta -r 342-944`
### Help:
 Further help can be showed by running without additional parameters:
    `java -jar InDelFixer.jar`

## PREREQUISITES COMPILE (only for dev):
 - Maven 3 (http://maven.apache.org/)

## INSTALL (only for dev):
    cd InDelFixer
    mvn clean package
    java -jar InDelFixer/target/InDelFixer.jar

# CONTACT:
    Armin Töpfer
    armin.toepfer (at) gmail.com
    http://www.bsse.ethz.ch/cbg/people/toepfera

# LICENSE:
 GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0