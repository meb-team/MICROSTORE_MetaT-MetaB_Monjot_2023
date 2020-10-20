
# **ANALYSE DE MÉTABARCODING**

## DATA

Les reads paired end (R1 et R2) sont placés dans un répertoire identifié par le nom de l'amplicon correspondant : "AC" par exemple.
Ces répertoires ("AC","FO", etc) sont placés dans le répertoire "reads".

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

## Lancement rapide : Pre-process

Une table contenant le nom ("AC", "FO", etc), la condition (DNOG, DJAP, etc) et la région (V4 ou V9) correspondant à chaque amplicon doit être réalisée.

* Le Genoscope fournit la table suivante : "data-inf.txt".
* Le répertoire "reads" et la table "data-inf.txt" sont placés dans le répertoire "rawdata" situé dans "Microstore-metabarcoding".

### 1. Pre-process et installation de PANAM

Définir le répertoire : `cd Microstore-metabarcoding/script/`

Lancer le script chargé de la mise en forme des reads et de l'installation de PANAM : `bash 0_Pre-process.sh`

    WARNING : Une ligne dans panam2.pl doit être modifiée ! :
    remplacer la ligne 169 : $panam_ini=$path_results."/panam.ini";
    par : $panam_ini=$path_results.$ini; #GB200114

### 2. Initialisation de PANAM

Différents fichiers d'initialisation ".ini" sont créés (à partir du fichier "test2-panam2.ini" présent dans PANAM2 et correspondant au jeux de données utilisés) :
    
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

    * V4-panam2-097.ini
        Parameters :
        * DOMAIN    eukaryota
        * PATH_RESULTS    V4-result-097
        * DEMUL_FOLDER    V4
        * FORWARD_PRIMER_SEQUENCE    GTGYCAGCMGCCGCGGTA
        * REVERSE_PRIMER_SEQUENCE    TTGGYRAATGCTTTCGC
        * CLUSTERING_CUTOFF    0.97
        * LOWER_SEQUENCE_LENGTH_CUTOF    200
        * MAX_SEQUENCE_LENGTH    500
        * MIN_OVERLAP_LENGTH    50
        * MISMATCH_OVERLAP    0
        * SCORE_FIX    0.90
        * MERGE_SEQ    perfect_match

    * V9-panam2-095.ini
        Parameters :
        * DOMAIN    eukaryota
        * PATH_RESULTS    V9-result-095
        * DEMUL_FOLDER    V9
        * FORWARD_PRIMER_SEQUENCE    TTGTACACACCGCCC
        * REVERSE_PRIMER_SEQUENCE    CCTTCYGCAGGTTCACCTAC
        * CLUSTERING_CUTOFF    0.95
        * LOWER_SEQUENCE_LENGTH_CUTOF    100
        * MAX_SEQUENCE_LENGTH    300
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

Ces fichiers d'initialisation ".ini" sont placés dans le répertoire "PANAM2"

### 3. ANALYSES PANAM 

Pour lancer PANAM2 : `cd ../dataPANAM/PANAM2`

* V4 avec seuil de clusterisation = 0.95 : `nohup perl panam2.pl -ini V4-panam2-095.ini > V4-095.out`

* V4 avec seuil de clusterisation = 0.97 : `nohup perl panam2.pl -ini V4-panam2-097.ini > V4-097.out`

* V4 avec seuil de clusterisation = 0.95 : `nohup perl panam2.pl -ini V4-panam2-095-pear.ini > V4-095-pear.out`

* V9 avec seuil de clusterisation = 0.95 : `nohup perl panam2.pl -ini V9-panam2-095.ini > V9-095.out`

* V9 avec seuil de clusterisation = 0.97 : `nohup perl panam2.pl -ini V9-panam2-097.ini > V9-097.out`

* V9 avec seuil de clusterisation = 0.97 : `nohup perl panam2.pl -ini V9-panam2-097-pear.ini > V9-097-pear.out`


### 4. ANALYSES R

* **Nécessite R version 3.6.3**

Un environnement conda est alors créé afin d'installer R 3.6.3 : `cd ../../script ; conda create -y -n REnv -c conda-forge r-base=3.6.3 ; conda activate REnv`

Un script permet d'installer les dépendences nécessaires aux analyses : `Rscript 7_Install_dependencies.R`

#### A. Analyse des duplicats

Une première analyse est réalisée afin de tester la différence entre amplicons d'un même duplicat : `Rscript 8_Analyse_Duplicat_AFC.R input output`

for exemple : `Rscript 8_Analyse_Duplicat_AFC.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch`

Renseigner Infos (Région / Rarecurve).

**PS : pour exécuter le script en nohup** : `nohup Rscript 8_Analyse_Duplicat_AFC.R input output region(V4 or V9) rarecurve(yes or no)`

for exemple : `nohup Rscript 8_Analyse_Duplicat_AFC.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch V4 yes > V4-095_Analyse_Duplicat.out`

#### B. Analyse de composition

L'analyse de composition et d'abondance est réalisée avec un second script : `Rscript 9_Composition_Rarefy.R input output`

par exemple :  `Rscript 9_Composition_Rarefy.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch`

Renseigner Infos  (Région / Rarecurve).

**PS : pour exécuter le script en nohup** : `nohup Rscript 9_Composition_Rarefy.R input output region(V4 or V9) rarecurve(yes or no)`

for exemple : `nohup Rscript 9_Composition_Rarefy.R ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch V4 yes > V4-095_Composition_Rarefy.out`
