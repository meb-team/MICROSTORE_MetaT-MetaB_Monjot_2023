
# ANALYSE DE MÉTABARCODING

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


## 1. TRI DES AMPLICONS V4 ET V9

Une table contenant le nom ("AC", "FO", etc), la condition (DNOG, DJAP, etc) et la région (V4 ou V9) correspondant à chaque amplicon doit être réalisée.

* Le Genoscope fournit la table suivante : "data-inf.txt".
* Le répertoire "reads" et la table "data-inf.txt" sont placés dans le répertoire "Microstore-metabarcoding".

Définir le répertoire : `cd Microstore-metabarcoding/`

Formater la table "data-inf.txt" au format voulu : `cd rawdata ; cat data-inf.txt | awk -F"\t" '{ print $2"_"$3"_"$5 }' | awk -F"_" '{ print $2"_"$3"_"$NF }' | cut -c1-10 | awk -F"_" '{ print $1"_"$2"_"$3"_" }' > sortV4-V9`

"sortV4-V9" définit tous les échantillons de la manière suivante : "Nom_Condition_Région_" (_exemple : "AC_DJAG_V9__" )

PS : Après certaines vérifications, l'amplicon "FR" noté par le génoscope comme amplicon de la région V4 est en réalité un amplicon de la région V9 : 
* Longueur des reads R1 : **151pb**
* Longueur des reads R2 : **151pb**
* Les reads ne "mergent" pas avec une longueur minimum de **200 pb** et un overlap de **50 pb**.

Cette information est alors modifiée manuellement dans la table "sortV4-V9" :  `sed -i 's/FR_DJAG_V4/FR_DJAG_V9/g' sortV4-V9`

Finalement, le tri (V4 vs V9) est réalisé avec le script "1_sortV4-9.sh" : `cd ../script ; bash 1_sort_V4-9.sh`

## 2. CONCATÉNER LES AMPLICONS RÉALISÉS À TRAVERS DEUX FLOW CELLS

Certains amplicons sont représentés par plusieurs fichiers de reads R1 et plusieurs fichiers de reads R2 car ils ont été séquencés avec plusieurs flow cells. Ces fichiers doivent être concaténés afin de n'avoir qu'un seul fichier de reads R1 et un seul fichier de reads R2 : `bash 2_concate.sh`

## 3. PRÉPARATION DU RÉPERTOIRE "dataPANAM"

Afin de débuter les analyses, les reads V4 et V9 sont copiés dans un répertoire "dataPANAM" :  `bash 3_sort.sh`

## 4. RE-IDENTIFICATION DES READS POUR PANAM

PANAM utilise la nomenclature suivante : "x_y_R1_z.fastq.gz" et "x_y_R2_z.fastq(.gz)". Il est alors nécessaire de renommer les reads.

_Par exemple :_
* _les reads :	"CIN_DMOSTA_2_1_HV2Y7BCX2.12BA289_clean.fastq"_
* _sont renommés en : "DMOSTA_2_R1_clean.fastq"_

Les reads sont renommés grâce au script "4_rawdata-rename.sh" : `bash 4_rawdata-rename.sh`

PS : les reads seront zippés ultérieurement pour être utilisés dans PANAM.

## 5. VÉRIFICATION AVEC VSEARCH (ÉTAPE OPTIONNELLE)

* **Nécessite vsearch**

Avant de lancer PANAM, la présence des primers utilisés ainsi que le "merging" des reads sont vérifiés. Selon la région (V4 ou V9), les paramètres suivants sont différents :
* Longueur minimale de l'overlap : **-fastq_minovlen**
* Longueur maximale des reads : **-fastq_maxmergelen**
* Longueur minimale des reads : **-fastq_minlen**

### Reads V4 

Paramètres d'assemblage :
* Longueur maximale V4 : **500pb**
* Longueur minimale V4 : **200pb**
* Longueur minimale de l'overlap : **50pb**
    
Les vérifications sont réalisées avec le script "5_testV4-findPrimer-assemblage.sh" : `nohup bash 5_testV4-findPrimer-assemblage.sh`

Les résultats sont visibles dans le fichier "V4-assembly.out" présent dans le répertoire "V4-testPrimer" dans "dataPANAM" suite à la commande suivante : `mv nohup.out ../dataPANAM/V4-testPrimer/V4-assembly.out`

Ce script assemble les reads R1 et R2 et supprime les chimères. Les séquences "nonchimera-" sont placées dans le répertoire "V4-testPrimer" dans "dataPANAM". Le pourcentage de séquences contenant les primers est obtenu avec le script "6_sort-onlyPrimerV4.sh" : `bash 6_sort-onlyPrimerV4.sh`

Pour cette étape, les primers suivants sont recherchés :

**F : 18SV4 515F**

    F :   5" GTG[CT]CAGC[AC]GCCGCGGTA "3
        *	3" CAC[GA]GTCG[TG]CGGCGCCAT "5
        *	5" TACCGCGGC[GT]GCTG[AG]CAC "3
    PS : `grep "^GTG[CT]CAGC[AC]GCCGCGGTA.*GCGAAAGCATT[CT][AG]CCAA$" nonchimera-*OSTA_2_R1_clean.fastq.fasta | wc -l`

**R : 18SV4 941R**

    R :   5" TTGG[CT][AG]AATGCTTTCGC "3
        *	3" AACC[GA][TC]TTACGAAAGCG "5
        *	5" GCGAAAGCATT[CT][AG]CCAA "3
    PS : `grep "^TTGG[CT][AG]AATGCTTTCGC.*TACCGCGGC[GT]GCTG[AG]CAC$" nonchimera-*OSTA_2_R1_clean.fastq.fasta | wc -l`

Les statistiques sont visibles dans le fichier "V4-primer.out" présent dans le répertoire "V4-testPrimer".

### V9 Reads

Paramètres d'assemblage :
* Longueur maximale V9 : **300pb**
* Longueur minimale V9 : **100pb**
* Longueur minimale de l'overlap : **50pb**

Les vérifications sont réalisées avec le script "5_testV9-findPrimer-assemblage.sh" : `nohup bash 5_testV9-findPrimer-assemblage.sh`

Les résultats sont visibles dans le fichier "V9-assembly.out" présent dans le répertoire "V9-testPrimer" dans "dataPANAM" suite à la commande suivante : `mv nohup.out ../dataPANAM/V9-testPrimer/V9-assembly.out`

Ce script assemble les reads R1 et R2 et supprime les chimères. Les séquences "nonchimera-" sont placées dans le répertoire "V9-testPrimer" dans "dataPANAM". Le pourcentage de séquences contenant les primers est obtenu avec le script "6_sort-onlyPrimerV9.sh" : `bash 6_sort-onlyPrimerV9.sh`

Pour cette étape, les primers suivants sont recherchés :

**F : 18SV9 1389F**

    F :   5" TTGTACACACCGCCC "3
        *	3" AACATGTGTGGCGGG "5
        *	5" GGGCGGTGTGTACAA "3
    PS : `grep "^TTGTACACACCGCCC.*GTAGGTGAACCTGC[AG]GAAGG$" nonchimera-*OSTA_2_R1_clean.fastq.fasta | wc -l`

**R : 18SV9 1510R**

    R :   5" CCTTC[CT]GCAGGTTCACCTAC "3
        *	3" GGAAG[GA]CGTCCAAGTGGATG "5
        *	5" GTAGGTGAACCTGC[AG]GAAGG "3
    PS : `grep "^CCTTC[CT]GCAGGTTCACCTAC.*GGGCGGTGTGTACAA$" nonchimera-*OSTA_2_R1_clean.fastq.fasta | wc -l`

Les statistiques sont visibles dans le fichier "V9-primer.out" présent dans le répertoire "V9-testPrimer".

## 6. INSTALLATION ET INITIALISATION DE PANAM 

Installation de PANAM dans le répertoire "dataPANAM"  : 

`cd ../dataPANAM ; git clone https://github.com/panammeb/PANAM2.git`

`cd PANAM2 ; perl setup.pl`

`chmod 777 bd_ssrna`

    WARNING : Une ligne dans panam.pl doit être modifiée ! :
    remplacer la ligne 169 : $panam_ini=$path_results."/panam.ini";
    par : $panam_ini=$path_results.$ini; #GB200114

Les reads V4 et V9 sont déplacés dans le répertoire PANAM2 : `cd .. ; mv V4 PANAM2/V4 ; gzip PANAM2/V4/*.fastq ; mv V9 PANAM2/V9 ; gzip PANAM2/V9/*.fastq`

Différents fichiers d'initialisation ".ini" sont créés (à partir du fichier "test2-panam2.ini" présent dans PANAM2) :
	
    * V4-panam2-095.ini (√)
		Parameters :
		* DOMAIN	eukaryota
		* PATH_RESULTS	V4-result-095
		* DEMUL_FOLDER    V4
		* FORWARD_PRIMER_SEQUENCE	GTGYCAGCMGCCGCGGTA
		* REVERSE_PRIMER_SEQUENCE	TTGGYRAATGCTTTCGC
		* CLUSTERING_CUTOFF	0.95
		* LOWER_SEQUENCE_LENGTH_CUTOF	200
		* MAX_SEQUENCE_LENGTH	500
		* MIN_OVERLAP_LENGTH	50
		* MISMATCH_OVERLAP	0
		* SCORE_FIX	0.90
		* MERGE_SEQ	perfect_match

	* V4-panam2-097.ini
		Parameters :
		* DOMAIN	eukaryota
		* PATH_RESULTS	V4-result-097
		* DEMUL_FOLDER    V4
		* FORWARD_PRIMER_SEQUENCE	GTGYCAGCMGCCGCGGTA
		* REVERSE_PRIMER_SEQUENCE	TTGGYRAATGCTTTCGC
		* CLUSTERING_CUTOFF	0.97
		* LOWER_SEQUENCE_LENGTH_CUTOF	200
		* MAX_SEQUENCE_LENGTH	500
		* MIN_OVERLAP_LENGTH	50
		* MISMATCH_OVERLAP	0
		* SCORE_FIX	0.90
		* MERGE_SEQ	perfect_match

	* V9-panam2-095.ini
		Parameters :
		* DOMAIN	eukaryota
		* PATH_RESULTS	V9-result-095
		* DEMUL_FOLDER    V9
		* FORWARD_PRIMER_SEQUENCE	TTGTACACACCGCCC
		* REVERSE_PRIMER_SEQUENCE	CCTTCYGCAGGTTCACCTAC
		* CLUSTERING_CUTOFF	0.95
		* LOWER_SEQUENCE_LENGTH_CUTOF	100
		* MAX_SEQUENCE_LENGTH	300
		* MIN_OVERLAP_LENGTH	50
		* MISMATCH_OVERLAP	0
		* SCORE_FIX	0.90
		* MERGE_SEQ	perfect_match

	* V9-panam2-097.ini (√)
		Parameters :
		* DOMAIN	eukaryota
		* PATH_RESULTS	V9-result-097
		* DEMUL_FOLDER    V9
		* FORWARD_PRIMER_SEQUENCE	TTGTACACACCGCCC
		* REVERSE_PRIMER_SEQUENCE	CCTTCYGCAGGTTCACCTAC
		* CLUSTERING_CUTOFF	0.97
		* LOWER_SEQUENCE_LENGTH_CUTOF	100
		* MAX_SEQUENCE_LENGTH	300
		* MIN_OVERLAP_LENGTH	50
		* MISMATCH_OVERLAP	0
		* SCORE_FIX	0.90
		* MERGE_SEQ	perfect_match

Ces fichiers d'initialisation ".ini" sont placés dans le répertoire "PANAM2"

## 7. ANALYSES PANAM 

Pour lancer PANAM2 : `cd PANAM2`

* V4 avec seuil de clusterisation = 0.95 : `nohup perl panam2.pl -ini V4-panam2-095.ini > V4-095.out`

* V4 avec seuil de clusterisation = 0.97 : `nohup perl panam2.pl -ini V4-panam2-097.ini > V4-097.out`

* V9 avec seuil de clusterisation = 0.95 : `nohup perl panam2.pl -ini V9-panam2-095.ini > V9-095.out`

* V9 avec seuil de clusterisation = 0.97 : `nohup perl panam2.pl -ini V9-panam2-097.ini > V9-097.out`

## 8. ANALYSES R (Utilisées ici sur les données V4-095 (seuil de clusterisation fixé à 0.95%))

* **Nécessite R version 3.6.3**

Un environnement conda est alors créé afin d'installer R 3.6.3 : `cd ../../script ; conda create -y -n REnv -c conda-forge r-base=3.6.3 ; conda activate REnv`

Un script permet d'installer les dépendences nécessaires aux analyses : `Rscript Install_dependencies.R`

Une première analyse est réalisée afin de tester la différence entre amplicons d'un même duplicat : `Rscript Auto-analyse_AFC_OTUonly.R`

L'analyse de composition et d'abondance est réalisée avec un second script : `Rscript Script-Composition-OTU-Rarefy.R`

Les résultats d'analyses se trouvent dans le répertoire "Analyse-Composition-Rarefy"
