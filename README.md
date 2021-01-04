# **METABARCODING ANALYSIS**

**This is a workflow to process the metabarcoding sequencing data from the Microstore project (2018).**

## DATA

The compressed paired end reads (R1.fastq.gz and R2.fastq.gz) are placed in a directory identified by the name of the corresponding amplicon: "AC" for example.
These directories ("AC", "FO", etc.) are placed in "reads" directory.

    *******************************
        Sampling nomenclature
    *******************************
    Sampling name example : DJAG_04_2
    with :
    1. Molecule type        ADN    D
    2. Daylight        JOUR:J    NUIT:N
    3. Oxygen        OXIQUE:O    ANOXIQUE:A
    4. Fraction        G:10-50    P:INF10
    5. Date        04    06    09    11
    6. Replicate        RIEN OU 1    2

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

## PANAM2 ANALYSIS

A table containing the name ("AC", "FO", etc.), the condition (DNOG, DJAP, etc.), the date (04, 06, 09 or 11), the region (V4 or V9) and replicate ID (1 or 2) corresponding to each amplicon must be produced.

* Here, Genoscope provides the following table: "data-inf.txt".

The "reads" directory and the "data-inf.txt" table are placed in the "rawdata" directory located in "Microstore-metabarcoding".

### 1. Pre-processing and PANAM2 installation

Define current directory: `cd Microstore-metabarcoding/`

Run script responsible for formatting reads and installing PANAM2: `bash PANAM_preprocess.sh`

PANAM2 is installed in dataPANAM directory : `cd /dataPANAM/PANAM2/`

    WARNING : one line in panam2.pl must be modified ! :
    replace line 169 : $panam_ini=$path_results."/panam.ini";
    by : $panam_ini=$path_results.$ini; #GB200114

### 2. PANAM2 initialization

Different initialization files (".ini") are created (from the "test2-panam2.ini" file found in "PANAM2" and corresponding to the data sets used):
    
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

### 3. PANAM2 processing

Run PANAM2 :

* V4 with clustering threshold = 0.95: `nohup perl panam2.pl -ini V4-panam2-095.ini > V4-095.out`

* V9 with clustering threshold = 0.97: `nohup perl panam2.pl -ini V9-panam2-097.ini > V9-097.out`

## R ANALYSIS

### 1. R installation (+ dependencies)

Return in initial directory : `cd ../../`

* **Requires R version 3.6.3**

Miniconda is downloaded: `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

and installed: `bash Miniconda3-latest-Linux-x86_64.sh`
        
To finilize the miniconda installation, close the currently terminal and open a new one.
        
Define current directory: `cd Microstore-metabarcoding/`
        
Install R and dependencies: `bash R_setup.sh`

### 2. R initialization

Initialization file (".ini") is created (from the "file.ini" file found in "Microstore-metabarcoding"):

For example:
        
        ## Input OTU table [1]
            INPUT    dataPANAM/PANAM2/V4-result-unified-095-215/OTU_distribution_tax.txt
        ## Output file (enter "auto" for automatic result file) [2]
            OUTPUT    auto
        ## Region [3]
            REGION    V4
        ## Taxonomic level (Superphylum, Phylum or Class) [8]
            MODE    Phylum
        ## Division for Phylum or Class level [9]
            DIVISION    Fungi
        ## Filter (if yes enter filter mode : "Bokulich", 'Doubleton' or "OnlyOne" ; if no enter "no") [4]
            FILTER    OnlyOne
        ## Taxonomy mode (NN, LCA or BH) [5]
            TAXONOMY    NN
        ## Unify duplicat (yes or no) [6]
            UNIFY    no
        ## Rarefy global data (yes or no) [7]
            RAREFY    yes
        ## Rarecurve Calcul (yes or no) ? It make a while [10]
            RARECURVE    no
        ## Note ? (optional)
            IDENTIFIER    095-Vsearch215

### 3. R processing

Run R workflow: `bash R_process.sh file.ini`

Results are generated in "result" file.
