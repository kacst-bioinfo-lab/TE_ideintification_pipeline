# Repeat library construction


To construct reliable and comprehensive repeat libraries is a challenging task due to the variation in repeat structure and the difficulty of assembling repeats in genome sequences. As many elements vary considerably in genetic structure and sequence, the only means of achieving reliable results when identifying and annotating TEs is to practice complementary approaches. A flowchart describing our overall approach to TE identification is given in the Figure below. The specific methods for each type are detailed below. We employ _**de novo**_ signature-based detection programs that rely upon prior knowledge concerning the sharing between different TEs of standard architectural features necessary for the process of transposition. 

Examples of classification according to similarity to known TEs include records in databases like **Repbase** and protein profiles retrieved from the **Pfam** database. Unfortunately, only well-described TEs that have a robust structural signature can be discovered by these methods. Some TEs do not have such characteristics and thus cannot be distinguished by this approach. In contrast to homology-based methods, signature-based methods are less biased by similarity to the set of known elements.

![Pipeline](https://github.com/kacst-bioinfo-lab/TE_ideintification_pipeline/blob/c96653e55dde2ea0faba57d802f41a262f81f816/figure_1%20.jpg)

# Class 1
### LTR  (long terminal repeat) retrotransposons

LTR retrotransposons represent the largest genomic mass among all repeats. As a result, it is very important to collect this type of element with high confidence. We annotated LTR-RTs using [the PhyLTR pipeline](https://github.com/mcsimenc/PhyLTR/).

Briefly, the elements are collected using  [LTRharvest](http://www.genometools.org/tools/gt_ltrharvest.html)  and filtered by  [LTRdigest](http://www.genometools.org/tools/gt_ltrdigest.html)  LTRdigest searches for homologs in the putative LTR-RTs using HMMER3 and a set of TE-related pHMMs we provided from Pfam and GyDB. and other custom programs. Both LTRharvest and LTRdigest are in a package called  [GenomeTools](https://www.computer.org/csdl/trans/tb/2013/03/ttb2013030645-abs.html).


**creating suffix array**
Commands:
 ```sh
gt suffixerator -db Input.fna -indexname Input.fna.index -dna -suf -lcp -des -ssp
```
* -dna:          input is DNA sequence
* -suf:          output suffix array (suftab) to file
* -lcp:          output lcp table (lcptab) to file
* -des:          output sequence descriptions to file
* -ssp:          output sequence separator positions to file

### [LTRharvest](http://www.genometools.org/tools/gt_ltrharvest.html) 
Searches the input sequence for direct repeats (LTRs) that are separated by a given distance (default 1 kb) and outside of which are apparent target site duplications (TSDs). Candidates distinguished by  [LTRharvest](http://www.genometools.org/tools/gt_ltrharvest.html) .
 ```sh
gt ltrharvest -index Input.fna.index -gff3 Input.fna.ltrharvest.out.gff -seqids yes -v yes -mintsd 4 -maxtsd 20 -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -minlenltr 100 -maxlenltr 1000 -mindistltr 100 -maxdistltr 15000 -vic 60
```
* -minlenltr: specify minimum length for each LTR.
* -maxlenltr: specify maximum length for each LTR.
* -mindistltr: specify minimum distance of LTR startpositions.
* -maxdistltr: specify maximum distance of LTR startpositions.
* -mintsd: specify minimum length for each TSD.
* -maxtsd: specify maximum length for each TSD.
* -vic: specify the number of nucleotides (to the left and to the right that will be searched for TSDs and/or motifs around 5' and 3' boundary of predicted LTR retrotransposons.
* -mat: specify matchscore for extension-alignment.
* -mis: specify mismatchscore for extension-alignment.
* -ins: specify insertionscore for extension-alignment.
* -del: specify deletionscore for extension-alignment.
* -v: verbose mode.
* -seqids: use sequence descriptions instead of sequence numbers in GFF3 output.
* -gff3: specify GFF3 outputfilename.

### [LTRdigest](http://www.genometools.org/tools/gt_ltrdigest.html) 
LTRdigest searches for homologs in the putative LTR-RTs using HMMER3 and a set of TE-related pHMMs we provided from  [Pfam](http://pfam.xfam.org/) and  [GyDB](https://gydb.org/index.php/Main_Page).
 ```sh
gt -j 20 ltrdigest -matchdescstart -outfileprefix Input.fna.LTRdigest -hmms Dfam.hmms -seqfile Input.fna
```
* -matchdescstart: exactly match the sequence descriptions from the input files for the desired sequence IDs (in GFF3) from the beginning to the first whitespace.
* -outfileprefix: prefix for output files (e.g. 'foo' will create files called 'foo_*.csv' and 'foo_*.fas').
* -hmms: profile HMM models for domain detection in HMMER3 format.
* -seqfile: set the sequence file from which to take the sequences Began extracting LTR_retrotransposon sequences from LTRharvest GFF.

**Began extracting LTR_retrotransposon sequences from LTRharvest GFF**

[bedtools](https://bedtools.readthedocs.io/en/latest/)
```sh
Usage: bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>

bedtools getfasta -fi Input.fna -s -bed Input.LTRharvest_LTR_retrotransposons.gff > Input.LTRharvest_LTR_retrotransposons.fasta 
```
* -fi  Input FASTA file
* -fo  Output file (opt., default is STDOUT
* -bed  BED/GFF/VCF file of ranges to extract from -fi
* -s  Force strandedness. If the feature occupies the antisense, strand, the sequence will be reverse complemented.
## Classification
#### Dfam Classification
 Using [HMMER](http://hmmer.org/) ([Dfma)](https://dfam.org/home) Data base.

 ```sh
 Usage: nhmmer [options] <query hmmfile|alignfile|seqfile> <target seqfile>
 
nhmmer --tblout DfamClassification/Input.fna.nhmmer_DfamHits.table --incE 1e-05 -E 10 --cpu 20 Dfam_ERV_LTR.hmm LTRharvest_LTR_retrotransposons.fasta 
```
* --tblout <f> : save parseable table of hits to file
* --incE <x> : consider sequences <= this E-value threshold as significant  [0.01]  (x>0)
* -E <x> : report sequences <= this E-value threshold in output  [10.0]  (x>0)
 ##### Extracting best hits from nhmmer on Dfam results

 ```sh
nhmmer_table2columns.py < Input.fna.nhmmer_DfamHits.table > Input.LTR_retrotransposon_DfamBestHits.tab
 ```
  [nhmmer_table2columns.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/nhmmer_table2columns.py)
##### Adding best hits from nhmmer on Dfam results to LTRdigest GFF
 ```sh
gffAddAttr.py -gff Input.fna.LTRdigest.withORFs_gt_300bp.gff -attr dfamClassification -map Input.fna.LTR_retrotransposon_DfamBestHits.tab -mapKey ID -restrictType LTR_retrotransposon -replaceIfNone > Input.fna.LTRdigest.withDfam.gff 
```
[gffAddAttr.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/gffAddAttr.py)

#### RepBase Classification 
Using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (Basic Local Alignment Search Tool) [RepBase](https://www.girinst.org/repbase/) data base.
```sh
tblastx -db Repbase_ERV_LTR.fasta -query Input.LTRharvest_LTR_retrotransposons.fasta -evalue 1e-05 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" -num_threads 20 -max_hsps 25 Output.fna.tblastx_Repbase.tab 
```
* -db: BLAST database name
* -evalue: Expectation value (E) threshold for saving hits
* -outfmt: 7 can be additionally configured to produce a custom format specified by space delimited format specifiers, or by a token specified by the delim keyword.
* -max_hsps: Set maximum number of HSPs per subject sequence to save for each query



 ##### Extracting best hits from tblastx on Repbase results, using [best_blast_hit.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/best_blast_hit.py).
```sh
best_blast_hit.py < Output.fna.tblastx_Repbase.tab > Output.fna.LTR_retrotransposon_RepbaseBestHits.tab
```


##### adding best hits from tblastx on Repbase results to LTRdigest GFF:
```sh
gffAddAttr.py -gff Input.LTRdigest.withDfam.gff -attr repbaseClassification -map Input.LTR_retrotransposon_RepbaseBestHits.tab -mapKey ID -restrictType LTR_retrotransposon -replaceIfNone > Input.LTRdigest.withDfam.withRepbase.gff 
```
[gffAddAttr.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/gffAddAttr.py)



> **Note:** The **Full pipeline depsited in zenodo** Further details regarding Class1 TEs Identification, classification, divergent, scripts. please install the script. [![DOI](https://doi.org/10.5281/zenodo.4256534)].

## Non-LTR (long terminal repeat) retrotransposons
Here, we began with the recognized genomic coordinates of LTR-RTs identified in the previous step. These candidates were masked with maskfasta from [BedTools](https://bedtools.readthedocs.io/en/latest/content/installation.html) to avoid conflicts or duplicate hits. 
```sh
Usage: bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>

bedtools maskfasta -fi InputGenome.fa -bed LTRdigestClassifiedNoFP.gff -fo InputGenome_Msked.fa
 ```
Next, open reading frame sequences were extracted from the masked genome by applying the getorf tool from [EMBOSS v6.4.0.0.](http://emboss.open-bio.org/html/use/ch02s07.html) The minimum ORF size was set to 500 bp in anticipation of detecting the apyrimidinic endonuclease (APE) gene (which is 600--800 bp in 97\% of inspected non-LTR elements)
```sh
getorf -sequence InputGenome_Msked.fa -outseq InputGenome_Msked.fa -minsize 500 -find 2 
 ```
* -minsize: Minimum nucleotide size of ORF
* -find: 2 (Nucleic sequences between STOP codons)

the exploration of the genomic sequences with [MGEScan-non-LTR](https://github.com/COL-IU/mgescan) , which identifies and classifies non-LTR TEs in genomic sequences using probabilistic models based on the structure of the 12 established non-LTR TE clades. More precisely, we used MGEScan-non-LTR and hmmsearch from [HMMER](http://hmmer.org/)  with two separate hidden Markov model (HMM) profiles, one for the reverse transcriptase (RT) gene and one for the endonuclease (APE) gene, both of which are well conserved among non-LTR TEs.


##### Create the Python virtual environment for MGEScan
```sh
virtualenv $MGESCAN_VENV
source $MGESCAN_VENV/bin/activate
```
##### Installation:
```sh
git clone https://github.com/COL-IU/mgescan.git
cd mgescan
python setup.py install
```
##### Usage:
```sh
 mgescan nonltr <genome_dir> [--output=<data_dir>] [--mpi=<num>]
```
# Class 2
All eukaryotic DNA transposons reported so far belong to a single category of elements which use the so-called ``cut-and-paste'' transposition mechanism, except Helitrons, which transpose by rolling-circle replication. Here, we employed methodologies for the detection of DNA transposons in the studied genomes based on the initial identification of TIR, and non-autonomous elements such as miniature inverted-repeat elements (MITEs) and helitron. 

MITEs are DNA-based elements that have TIRs but lack a transposase gene, and their well-defined structural features make them suitable for discovery by computational approaches. a valuable tool for detecting MITEs in eukaryotic genomes, [MiteFinder](https://github.com/jhu99/miteFinder). this tool is capable of detecting both perfect and imperfect inverted repeats through a string matching approach. It computes a new function to cluster MITE sequences into different MITE families in several steps. 
 1. First, it builds a k-mer index and seeks inverted repeats. Then, all sequence candidates are distinguished by the presence of a TIR pair of default length and a TSD pair. 
 2. the scaffolds are divided into multiple sequence fragments that overlap by 800 bp, which is the maximum length of MITEs, to guarantee that all inverted repeats are identified. 
 3. pairs of TIRs having lengths in the range of 50-800 bp are retained, and the remainder used as seeds for MITE candidates in the next step. 
 4. identified sequences are compared with MITEs in the Repbase database using blastn. Those with high similarity are considered valid positives, and those with low similarity as false positives. For each MITE cluster, the sequence with the highest blast score was selected via [best_blast_hit.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/best_blast_hit.py) as the representative family sequence. The tool was executed with default parameters, except for the use of a confidence-score threshold of 0.5 to exclude low-confidence candidates.

### [MiteFinder](https://github.com/jhu99/miteFinder)
```sh
miteFinder -input Input.fna -output Output_Mite.fa -pattern_scoring pattern_scoring.txt -threshold 0.5
```
* -threshold  Threshold of removing mite candidates with low-confidence score.
* -pattern_scoring  The path of a scoring file.

**Clustering** ([CD-Hit](http://weizhong-lab.ucsd.edu/cd-hit/))
```sh
CD-Hit -i Output_Mite.fa -o Clustered-Output_Mite.fa -c 0.8
```
**Classification** ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download))
```sh
blastn -task blastn -query SpeciesMITEs.fa -subject Clustered-Output_Mite.fa -outfmt 6 > out.csv
```
**Extract Best hits** [best_blast_hit.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/best_blast_hit.py)
```sh
best_blast_hit.py < out.csv > bestHits
```
### [HelitronScanner](https://sourceforge.net/projects/helitronscanner/)
HelitronScanner relies on sequence matches between trained local combinational variables (LCVs) and genome . Specifically, it scores 5$'$ and 3$'$ termini based on a training set of published Helitrons, and then merges the coordinates and scores for putative Helitron-like sequences. Multiple stages are examined, using the HelitronScanner tool, First, scan the **Head/Tail** with a default threshold is 1, which allows any matches to LCVs. The location scores are further filtered in the **pairends** stage, considering a threshold of 6 for Helitron’s both ends.
```sh
java -jar HelitronScanner.jar scanHead -lcv_filepath -genome -buffer_size -output HeadOut
java -jar HelitronScanner.jar scanTail -lcv_filepath -genome -buffer_size -output TailOut
java -jar HelitronScanner.jar pairends -hs HeadOut -ts TailOut -o Pair.fa -ht 6 -tt 6
java -jar HelitronScanner.jar draw -pscore Pair.fa -genome genome.fa -output FinalHelitron.fa --pure
 ```
 **Clustering** ([CD-Hit](http://weizhong-lab.ucsd.edu/cd-hit/))
 ```sh
CD-Hit -i FinalHelitron.fa -o Clustered_Helitron.fa -c 0.8
```
**Classification** ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download))
```sh
blastn -task blastn -query SpeciesHelitron.fa -subject Clustered_Helitron.fa -outfmt 6 > out.csv
```
 **Extract Best hits** [best_blast_hit.py](https://github.com/mcsimenc/PhyLTR/blob/master/scripts/best_blast_hit.py)
```sh
best_blast_hit.py < out.csv > bestHits
```
## Combine Library
libraries generated using the programs mentioned above merged, pluse the TE repeat sequences present in study species can be extracte from Dfam Consensus-20170127 and RepBase-20181026 using the script queryRepeatDatabase shipped with [RepeatMasker](http://www.repeatmasker.org/RepeatMasker.html). The results of both steps above were combined and processed for masking and annotation using [RepeatMasker](http://www.repeatmasker.org/RepeatMasker.html). 


## Annotation, and mask the genome using [RepeatMasker](http://www.repeatmasker.org/RepeatMasker.html)
```sh
RepeatMasker -nolow -norna -a -s -gff -cutoff 225 -no_is -pa 40 -lib Combined_Lib.fa Input_Genome.fa -dir Output
 ```
* -cutoff: Sets cutoff score for masking repeats when using -lib
* -no_is: Skips bacterial insertion element check
* -lib: Allows use of a custom library
* -norna: Does not mask small RNA (pseudo) genes 
* -nolow: Does not mask low_complexity DNA or simple repeats
* -a:Writes alignments in .align output file
* -s: Slow search; 0-5% more sensitive, 2-3 times slower than default
* -gff:Creates an additional Gene Feature Finding format output

 Additionally, the copy number of each TE and determined genome coverage from the RepeatMasker output files (.out), which correspond to the number of insertions identified in the masked genomes. The remaining unmasked portion of the genome is scanned using [RepeatModeler](http://www.repeatmasker.org/RepeatModeler.html) with default settings to detect any unclassified TEs such as TIRs that were missed by structure-based TE identification. 

## Detecting repetitive sequences missed by previouse Structure-based tools using [RepeatModeler](http://www.repeatmasker.org/RepeatModeler.html)
The remaining sequence is then used as the input for RepeatModeler:

```sh
BuildDatabase -name MaskedGenome Output.masked.fna
RepeatModeler -database MaskedGenome
 ```

 
Among the sequences generated by RepeatModeler, some are associated with identities and others are not. The script [repeatmodeler_parse.pl](http://www.hrt.msu.edu/uploads/535/78637/CRL_Scripts1.0.tar.gz) separates elements without identities into repeatmodeler_unknowns and the others in repeatmodeler_identities.
* Finally the candidate duplication removal assisted using  [seqket](https://bioinf.shenwei.me/seqkit/) rmdup. 
* Filtered library used as an input for reannotation using [RepeatMasker](http://www.repeatmasker.org/RepeatMasker.html) To estimate the final Transposable elements in the genome.
## Citation
Ibrahim, M.A., Al-Shomrani, B.M., Simenc, M., Alharbi, S.N., Alqahtani, F.H., Al-Fageeh, M.B. and Manee, M.M., 2021. Comparative analysis of transposable elements provides insights into genome evolution in the genus Camelus. BMC Genomics, 22(1), pp.1-16.
