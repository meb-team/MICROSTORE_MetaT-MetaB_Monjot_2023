# **METABARCODING ANALYSIS WITH DADA2**

**This is a workflow to process the metabarcoding sequencing data from the Microstore project (2018) with DADA2.**

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

## DADA2 ANALYSIS

A table containing the name ("AC", "FO", etc.), the condition (DNOG, DJAP, etc.), the date (04, 06, 09 or 11), the region (V4 or V9) and replicate ID (1 or 2) corresponding to each amplicon must be produced.

* Here, Genoscope provides the following table: "data-inf.txt".

The "reads" directory and the "data-inf.txt" table are placed in the "rawdata" directory located in "Metabarcoding_DADA2".

### 1. Pre-processing and reads re-orientation

First, define current directory: `cd Metabarcoding_DADA2/`

Second, to start DADA2 analysis, the V4 and V9 reads must be pooled according to their replicates ID. 

To complete this, run script 1 : 

    * script: DADA2_1_preprocess.sh

For example : `bash DADA2_1_preprocess.sh`

At the end of this stage, four result files will be generated : "V4", "V9", "V4-unified" and "V9-unified" ; respectively reads V4 and V9 not pooled and the reads V4 and V9 pooled. These files will be created in reads directory located in dataDADA2 directory (also create at this stage).

Third, reads must be re-orient for ASV processing. In fact, Illumina adaptators are randomly linked to the DNA fragments during de sequencing process and the reads R1 and R2 are represent by approximately 50% of forward reads and  50% of reverse reads.

To complete this, run script 2 :
    
    * script: DADA2_2_reorient.sh
    * argument 1: input file for reorientation (V4, V9, V4-unified or V9-unified)
    * argument 2: number of threads to process data
    * argument 3: forward primer
    * argument 4: reverse primer

For example: `bash DADA2_2_reorient_py.sh V4-unified 16 "GTG[CT]CAGC[AC]GCCGCGGTA" "TTGG[CT][AG]AATGCTTTCGC"`

### 2. DADA2 initialization

### 1. R installation (+ dependencies)

Return in initial directory : `cd ../`

Miniconda is downloaded: `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

and installed: `bash Miniconda3-latest-Linux-x86_64.sh`

To finilize the miniconda installation, close the currently terminal and open a new one.

Define current directory: `cd Metabarcoding_DADA2/`

Install R and dependencies: `bash R_3_setup.sh`

### 2. R initialization

Initialization file (".ini") is created (from the "V4-DADA2.ini" file found in "Metabarcoding_DADA2"):

For example:
        
        ## Raw data directory [1]
        INPUT    V4-unified-correct-paired

        ## Database [2]
        DATABASE    pr2_version_4.14.0_SSU_dada2.fasta.gz

        ## Result directory [3]
        RESULT    V4-unified-correct-paired-out

        ## Minimum read length [4]
        MINLEN    200

        ## Maximum read length [5]
        MAXLEN    500

        ## Maximum N nucleotide in reads [6]
        MAXN    0

        ## Minimum of the overlap region between forward and reverse reads [7]
        MINOVERLAP    50

        ## Mismatch value accepted on the overlap region between forward and reverse reads [8]
        MAXMISMATCH    0

        ## Number of threads to process data [9]
        NTHREADS    16

        ## Primer Forward [10]
        FWD    GTGYCAGCMGCCGCGGTA

        ## Primer Reverse [11]
        REV    TTGGYRAATGCTTTCGC

        ## Region [12]
        REGION    V4

        ## Filter (if yes enter filter mode : "Bokulich", "Singleton", "Doubleton" or "OnlyOne" ; if no enter "no") [13]
        FILTER    Singleton

        ## Rarefy global data (yes or no) [14]
        RAREFY    yes

        ## Rarecurve Calcul (yes or no) ? It make a while [15]
        RARECURVE    yes

        ## Process Composition script (yes or no) [16]
        COMPOSITION    yes

        ## Output for the Composition script [17]
        OUTCOMP    V4-unified-correct-paired-out-compo

### 3. R processing

Run DADA2 workflow: `bash R_4_DADA2_process.sh V4-DADA2.ini`

Run Taxonomic analysis script: `bash R_5_Taxonomic_analysis.sh V4-DADA2.ini`

The previous step annotate taxonomic references using marine trait table described in Ramond et al. 2019. The resulting table is generated in "result" file which is located in dataDADA2 directory.

Prepare your trait table using the previously generated table and place it on rawdata directory in tsv format or use our trait table already present in rawdata directory.

Run Functionnal analysis script:
    
    * script: R_6_Functional_analysis.sh
    * argument 1: name of the trait table in tsv format
    * argument 2: name of the taxonomic analysis results directory
    * argument 3: Region (V4 or V9)

For example: `bash R_6_Functional_analysis.sh Trait_Table_Final.tsv V4-unified-correct-paired-out-compo V4`

### 4. Implementation of metatranscriptomic data

First, inputs and metadata file must be placed in rawdata directory.

* Genoscope provides the following metadata table: "data-inf-metaT.txt".

* We use unigene table and taxonomy table provided by Damien Courtine: respectively "main_table.mapping.unique.read_per_kb.noHuman.noConta.noMetazoa.annot.tsv" and "table_taxonomy.perUnigene.allUnigenes.tsv" ; KO id definition table generated from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg and the database used in the trait-based approcach (i.e. PR2).

Initialization file (".ini") is then created (from the "MetaT.ini" file found in "Metabarcoding_DADA2"):

For example:
    
        ## Enter unigene table (.tsv) (located in rawdata directory) [1]
        INPUT    main_table.mapping.unique.read_per_kb.noHuman.noConta.noMetazoa.annot.tsv

        ## Result path (the same of the composition script output) [2]
        OUTPUT    V4-unified-correct-paired-out-compo

        ## Unigene Taxonomy path (located in rawdata directory) [3]
        TAX    table_taxonomy.perUnigene.allUnigenes.tsv

        ## Database path used in the trait study (located in database directory) [4]
        DATABASE    pr2_version_4.14.0_SSU_dada2.fasta.gz


Run metatranscriptomic analysis script: `bash R_7_Metatranscriptomic_analysis.sh MetaT.ini`

Results are generated in "result" file wich is located in dataDADA2 directory.
