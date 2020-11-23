
# **MÉTABARCODING ANALYSIS**

## DATA

The paired end reads (R1 and R2) are placed in a directory identified by the name of the corresponding amplicon: "AC" for example.
These directories ("AC", "FO", etc.) are placed in the "reads" directory in "rawdata".

	*******************************
		Sampling nomenclature
	*******************************
	Sampling name example : DJAG_04_2
	with :
	1. Molecule type		ADN	D
	2. Daylight		JOUR:J	NUIT:N
	3. Oxygen		OXIQUE:O	ANOXIQUE:A
	4. Fraction		G:10-50	P:INF10
	5. Date		04	06	09	11
	6. Replicate		RIEN OU 1	2

	*******************************
		Run name nomenclature
	*******************************
	Run name example : 191029_MELISSE_D78BK
	with :
	1. Date : 191029
	2. Sequencer : MELISSE
	3. Flowcell id : D78BK
	Other sequencer names : MELISSE MIMOSA PLATINE

	*******************************
		Paired files nomenclature
	*******************************
	Name file examples :
	CIN_CBOSTA_1_1_D78BK.12BA207_clean.fastq.gz
	CIN_CBOSTA_1_2_D78BK.12BA207_clean.fastq.gz
	with :
	1. Project code : CIN
	2. Material name : CB
	3. Sequencing type : S = Solexa
	4. Bank type : TA : TA - Targeted DNAseq
	5. Lane number : 1
	6. Read number : 1
	7. Flowcell id : D78BK
	8. Index : 12BA207

	*******************************
		Sequence header nomenclature
	*******************************
	Sequence name example :
	@M1:A6U0C:1:1101:17317:1262/2
	@H3:C29DNACXX:2:1101:1349:2239/1
	with :
	1. Sequencer id : @M1 ; M = Miseq; H = Hiseq
	2. Flowcell id : A6U0C
	3. Lane number : 1
	4. Tile number : 1101
	5. x coordinate on the tile : 17317
	6. y coordinate on the tile : 1262

## Pre-process

A table containing the name ("AC", "FO", etc.), the condition (DNOG, DJAP, etc.) and the region (V4 or V9) corresponding to each amplicon must be produced.

* Genoscope provides the following table: "data-inf.txt".
* The "reads" directory and the "data-inf.txt" table are placed in the "rawdata" directory located in "Microstore-metabarcoding".

### 1. Pre-process and PANAM2 installation

Define current directory: `cd Microstore-metabarcoding/script/`

* **Requires R version 3.6.3**

A conda environment is created to install R 3.6.3: `conda create -y -n REnv -c conda-forge r-base=3.6.3 ; conda activate REnv`

A script is used to install the dependencies necessary for the analyzes: `Rscript 0_Install_dependencies.R`

Run script responsible for formatting reads and installing PANAM2: `bash 1_Pre-process.sh`

    WARNING : one line in panam2.pl must be modified ! :
    replace line 169 : $panam_ini=$path_results."/panam.ini";
    by : $panam_ini=$path_results.$ini; #GB200114

### 2. PANAM2 initialization

Different initialization files (".ini") are created (from the "test2-panam2.ini" file found in PANAM2 and corresponding to the data sets used):
    
    * V4-panam2-095.ini (√)
        Parameters :
        * DOMAIN    eukaryota
        * PATH_RESULTS    V4-result-095
        * DEMUL_FOLDER    V4
        * FORWARD_PRIMER_SEQUENCE    GTGYCAGCMGCCGCGGTA
        * REVERSE_PRIMER_SEQUENCE    TTGGYRAATGCTTTCGC
        * CLUSTERING_CUTOFF    0.95
        * LOWER_SEQUENCE_LENGTH_CUTOF    200
        * MAX_SEQUENCE_LENGTH    500
        * MIN_OVERLAP_LENGTH    50
        * MISMATCH_OVERLAP    0
        * SCORE_FIX    0.90
        * MERGE_SEQ    perfect_match

    * V9-panam2-097.ini (√)
        Parameters :
        * DOMAIN    eukaryota
        * PATH_RESULTS    V9-result-097
        * DEMUL_FOLDER    V9
        * FORWARD_PRIMER_SEQUENCE    TTGTACACACCGCCC
        * REVERSE_PRIMER_SEQUENCE    CCTTCYGCAGGTTCACCTAC
        * CLUSTERING_CUTOFF    0.97
        * LOWER_SEQUENCE_LENGTH_CUTOF    100
        * MAX_SEQUENCE_LENGTH    300
        * MIN_OVERLAP_LENGTH    50
        * MISMATCH_OVERLAP    0
        * SCORE_FIX    0.90
        * MERGE_SEQ    perfect_match

These initialization files are placed in the "PANAM2" directory

### 3. PANAM2 ANALYSIS

Run PANAM2 :

* V4 with clustering threshold = 0.95: `nohup perl panam2.pl -ini V4-panam2-095.ini > V4-095.out`

* V9 with clustering threshold = 0.97: `nohup perl panam2.pl -ini V9-panam2-097.ini > V9-097.out`

### 4. R ANALYSIS

#### A. Duplicate analysis

First, return in script directory : `cd ../../script`

A first analysis is carried out in order to test the difference between amplicons of the same duplicate: `Rscript 8_Analyse_Duplicat_AFC.R input output`

for exemple: `Rscript 8_Analyse_Duplicat_AFC.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch`

Enter Info (region / 0.0005filter / rarefy).

**PS : to execute script silently** : `nohup Rscript 8_Analyse_Duplicat_AFC.R input output region(V4 or V9) 0.0005 filter(yes or no) rarefy(yes or no)`

for example : `nohup Rscript 8_Analyse_Duplicat_AFC.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch V4 yes yes > V4-095_Analyse_Duplicat.out`

#### B. Composition Analysis

The composition and abundance analysis is performed with a second script: `Rscript 9_Analyse_Composition.R input output`

for example :  `Rscript 9_Analyse_Composition.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch`

Enter Info (region / mode / group[optional] / 0.0005filter / taxonomy / unify /rarefy).

**PS : to execute script silently** : `nohup Rscript 9_Analyse_Composition.R input output region(V4 or V9) 0.0005filter(yes or no) taxonomy(NN, LCA, or Best_HIT) unify(yes or no) rarefy(yes or no) mode(Superphylum or Phylum) phylum(Fungi, Alveolata, etc)[optional]`

for example : `nohup Rscript 9_Analyse_Composition.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch V4 yes LCA yes yes Phylum Alveolata > V4-095_Analyse_Composition.out`

#### C. Rarefaction curves and diversity indices

Rarefaction curves and diversity indices can be calculated using a third script : `Rscript 10_Calcul_Rarecurve.R input output`

for example :  `Rscript 10_Calcul_Rarecurve.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch`

Enter Info (region / 0.0005filter).

**PS : to execute script silently** : `nohup Rscript 10_Calcul_Rarecurve.R input output region(V4 or V9) 0.0005filter(yes or no)`

for example : `nohup Rscript 10_Calcul_Rarecurve.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch V4 yes > V4-095_Composition_Rarefy.out`
