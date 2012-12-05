# CONTENT:
This java command line application aligns Next-Generation Sequencing (NGS) reads
to a set of reference sequences, by a fast k-mer matching and removes indels, causing
frame shifts. In addition, only a specific region can be considered.

## PREREQUISITES TO RUN:
 - JDK 7 (http://jdk7.java.net/)
 - Get latest version: https://github.com/armintoepfer/InDelFixer/downloads

## RUN:
 (454/Roche)
    java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta
 or
 (Plain sequences, for example PacBio)
    java -jar InDelFixer.jar -i libCase102.fasta -g referenceGenomes.fasta
 or
 (Illumina paired end)
    java -jar InDelFixer.jar -i libCase102_R1.fastq -ir libCase102_R2.fastq -g referenceGenomes.fasta

## USE / SHOW HELP:
    java -jar InDelFixer.jar

## Extract reads of particular regions
    java -jar InDelFixer.jar -i libCase102.sff -g referenceGenomes.fasta -r 342-944

## PREREQUISITES COMPILE (only for dev):
 - Maven 3 (http://maven.apache.org/)

## INSTALL (only for dev):
 cd InDelFixer
 mvn clean package
 java -jar InDelFixer/target/InDelFixer.jar --help

# CONTACT:
 Armin Töpfer
 armin.toepfer (at) gmail.com
 http://www.bsse.ethz.ch/cbg/people/toepfera

# LICENSE:
 GNU GPLv3 http://www.gnu.org/licenses/gpl-3.0