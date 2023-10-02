# **Functional diversity of microbial eukaryotes in a meromictic lake: coupling between metatranscriptomic and a trait-based approach**

**This is a workflow to reproduce analysis conduced in Monjot et al., 2023**

First, clone github repository: `git clone https://github.com/meb-team/MICROSTORE_MetaT-MetaB_Monjot_2023.git`

Second, define current directory: `cd MICROSTORE_MetaT-MetaB_Monjot_2023`

                                        *******************************
                                            Directory organization
                                        *******************************
    Microstore_Analysis_Monjot_et_al._2023 
    |-> rawdata (sub-directory for reads, traits table, metadata)
        |-> metadata_metaB (metadata of metabarcoding data)
        |-> metadata_metaT (metadata of metatranscriptomic data)
        |-> BTTwellannoted.csv (Marine traits table by Ramond et al., 2019)
        |-> Table_Supp_1.tsv (own traits table (Monjot et al., 2023))
        |-> ko_to_hierarchy.txt (KO id definition table from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg)
    |-> script (sub-directory for minor scripts)
        |-> 0_download_metaB_data.sh (script for downloading metabarcoding raw reads)
        |-> 1_Pre-process.sh (pre-processing script)
        |-> 2_Orient_reads_parallel.sh (script to launch parallelization of reads reorientation script)
        |-> 2_Orient_reads.py (reads reorientation script)
        |-> 3_Install_dependencies.R (dependencies installation script)
        |-> 4_Dada2.R (DADA2 pipeline script)
        |-> 5_Analyse_Composition_ASV_DADA2.R (script to analyze taxonomic diversity)
        |-> 6_Preview_Tax_ref.R (script to extract taxonomic ref in order to launch trait-based analysis)
        |-> 7_Calcul_Rarecurve.R (script to process rarefaction curve and diversity index)
        |-> 8A_Kmean_Clusterization.R (script to clusterize taxonomic references using custom traits table)
        |-> 8B_Analyse_Function_Table.R (script to start trait-based analysis using previously defined cluster names)
        |-> 9_Duplicat_Metatrans_analize.R (script to start comparison between duplicates of conditions)
        |-> 10_KO_Metatrans_DESeq2.R (DESeq2 pipeline script)
        |-> environment_REnv_Monjot_2023A.yml (conda environment)
        |-> REnv_Monjot_2023A_packages (local repository to R packages installation)
    |-> Preprocess_setup.sh (script to launch conda environment setup)
    |-> Downloading_metaB_rawdata.sh (script to launch  metabarcoding raw reads download)
    |-> Downloading_metaT_rawdata.sh (script to launch  metatranscriptomic catalog genes download)
    |-> DADA2_1_preprocess.sh (script to launch the preprocessing stage)
    |-> DADA2_2_reorient_py.sh (script to launch raw reads reorientation)
    |-> R_3_setup.sh (script to launch R dependencies installation scripts)
    |-> R_4_DADA2_process.sh (script to launch DADA2 pipeline)
    |-> R_5_Taxonomic_analysis.sh (script to launch Taxonomic diversity analysis)
    |-> R_6A_Kmean_Clusterization.sh (script to launch kmean clusterization)
    |-> R_6B_Functional_analysis.sh (script to launch Trait-based analysis)
    |-> R_7_Metatranscriptomic_analysis.sh (script to launch DESeq2 pipeline and metatranscriptomic analysis)
    |-> Retrieve_Figures.sh (script to retrieve published figures from result directory)
    |-> V4-DADA2.ini (initialization file for metabarcoding analysis)
    |-> MetaT.ini (initialization file for metatranscriptomic analysis)
    

## Download metabarcoding data

Raw reads is downloaded from ENA archive under PRJEB61527 accession number: `bash Downloading_metaB_rawdata.sh`

The compressed paired end reads (*R1.fastq.gz* and *R2.fastq.gz*) are automatically placed in a *reads/* sub-directory in *rawdata/* directory.

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
    

A table containing the name ("DM", "DN", etc.), the condition (DNOG, DJAP, etc.), the date (04, 06, 09 or 11), the region (V4) and replicate ID (1 or 2) corresponding to each amplicon is needed.

* We provides the following table: *metadata_metaB* which is located in the *rawdata/* directory.

## Preprocessing and dependencies installation

### 1. Conda environment installation

Install miniconda following the standard procedure (https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html)

Then, install conda environment with the following script: `bash Preprocess_setup.sh`

    * This takes 1 min on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

This installs the following tools:

    * KronaTools v2.8
    * cutadapt v4.1
    * r-base v4.2.2
    * imagemagick v7.1.0_52
    * python v3.9.15
    * perl v5.32.1
    * cmake v3.22.1
    * parallel v20230722
    * numerous libraries (details in script/environment_REnv_Monjot_2023A.yml)

### 2. Pre-processing and reads re-orientation
    
To start DADA2 analysis, the reads must be pooled according to their replicates ID. 

* script: DADA2_1_preprocess.sh

Run pre-processing script : `bash DADA2_1_preprocess.sh`

    * This takes just a few seconds on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

At the end of this stage, some result files are generated : *V4* and *V4-unified* ; respectively reads V4 not pooled and the reads V4 pooled. These files are created in *reads/* sub-directory located in *dataDADA2/* directory (also created at this stage).

Reads must be re-orient for ASV processing. In fact, Illumina adaptators are randomly linked to the DNA fragments during sequencing process and the reads R1 and R2 are represented by approximately 50% of forward reads and 50% of reverse reads.

To complete this, run following script :
    
* script: DADA2_2_reorient.sh
* argument 1: input file for reorientation (V4 or V4-unified)
* argument 2: number of threads to process data
* argument 3: forward primer
* argument 4: reverse primer

=> `bash DADA2_2_reorient_py.sh V4-unified 16 "GTG[CT]CAGC[AC]GCCGCGGTA" "TTGG[CT][AG]AATGCTTTCGC"`

    * This takes 1 min on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

### 3. R installation (+ dependencies) if necessary

Install R dependencies with the following script: `bash R_3_setup.sh`

    * This takes 77 min on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

This operation installs the following R packages:

    * GUniFrac              * reshape2                  * vegan
    * ggplot2               * SARTools                  * gtools
    * tidyr                 * ComplexHeatmap            * SciViews
    * dplyr                 * DESeq2                    * paletteer
    * cowplot               * dada2                     * cluster
    * ggrepel               * elementalist              * stringr
    * ggsci                 * ggh4x                     * NbClust
    * scales                * R.utils                   * BiocManager
    * varhandle             * scatterplot3d             * tibble
    * treemap               * car                       * svglite
    * VennDiagram           * rgl                       * treemapify
    * FactoMineR            * data.table                * hrbrthemes   
    * RColorBrewer          * devtools                  * psych
    * factoextra            * ggExtra                   * ggpubr
    * reshape2              * gplots

## Metabarcoding analysis

### 1. R initialization

Initialization file *.ini* is created from the *V4-DADA2.ini*:

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

### 2. R processing

#### A. DADA2 pipeline

Run DADA2 workflow:

* script: R_4_DADA2_process.sh
* argument 1: V4-DADA2.ini

=> `bash R_4_DADA2_process.sh V4-DADA2.ini`

    * This takes 84 min on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

#### B. Taxonomic diversity analysis

Run Taxonomic analysis script:

* script: R_5_Taxonomic_analysis.sh
* argument 1: V4-DADA2.ini

=> `bash R_5_Taxonomic_analysis.sh V4-DADA2.ini`

    * This takes 35 min on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

#### C. Trait-based analysis

Previous steps annotate taxonomic references using marine traits table provided by Ramond et al., 2019. The resulting table is generated in *result/* sub-directory which is located in *dataDADA2/* directory.

Prepare your traits table using the previously generated table and place it on *rawdata/* directory in *.tsv* format or use our trait table already present in *rawdata/* directory: *Table_Supp_1.tsv*.

Run kmean clusterization to define clusters:

* script: R_6A_Kmean_Clusterization.sh
* argument 1: name of the trait table in tsv format
* argument 2: name of the taxonomic analysis results directory
* argument 3: Region (V4 or V9)
    
=> `bash R_6A_Kmean_Clusterization.sh Table_Supp_1.tsv V4-unified-correct-paired-out-compo V4`

    * This takes just a few seconds on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

This operation produces a *cluster_name.tsv* table providing the number of clusters as well as a *Def_group_by_factor.svg* figure in *rawdata/* directory. You have to complete the *.tsv* file with names of clusters using the *Def_group_by_factor.svg* figure.

In Monjot et al. 2023, authors characterized clusters using modalities of each morpho-physio-phenological trait. The functional groups have been named consequently: 1&5) Parasites (PARA): characterized by their feeding strategy, symbiosis type (the majority has been described as having a parasitic lifestyle) and organic covers or naked; 2) Saprotrophs (SAP): characterized by saprotrophic feeding strategy, attached lifestyle and mainly absence of biotic interaction reported in the literature; 3) Heavy-cover- and Swimmer-photoautotrophs (HCOV and SWAT): characterized by plastids presence and osmotrophic ingestion mode with either mineral cover (e.g. siliceous) and non-swimming abilities or organic cover and swimming abilities; 4&7) Strict-heterotrophs (HET): characterized by phagotrophic feeding strategy and the absence of plastids; 6) Mixotrophs (MIXO): considered as mixoplankton (photo-osmo-phago-mixotrophs according to Mitra et al. (2023) and characterized by the presence of chloroplast, motility and their feeding strategy (Phagotrophic in majority); 8) Floater- and colonial-photoautotrophs (FLAT): characterized by non-swimming abilities, plastids presence and osmotrophic ingestion mode; 9) Endophyte (END): characterized by their attached life-style, feeding strategy, biotic interaction with plants and mostly organic covers.

For example:

    nbGroup	name
    Cluster 1	PARA
    Cluster 2	SAP
    Cluster 3	HCOV and SWAT
    Cluster 4	HET
    Cluster 5	PARA
    Cluster 6	MIXO
    Cluster 7	HET
    Cluster 8	FLAT
    Cluster 9	END

Then, run functional analysis script:
    
* script: R_6B_Functional_analysis.sh
* argument 1: name of the trait table in tsv format
* argument 2: name of the taxonomic analysis results directory
* argument 3: Region (V4 or V9)

=> `bash R_6B_Functional_analysis.sh Table_Supp_1.tsv V4-unified-correct-paired-out-compo V4`

    * This takes 3 min on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

## Metatranscriptomic analysis

Unigenes catalog was downloaded from ZENODO archive (available under https://doi.org/10.5281/zenodo.8376850 DOI link): `bash Downloading_metaT_rawdata.sh`

* This script downloads unigenes table and taxonomy table: respectively *main_table.mapping.unique.read_per_kb.noHuman.noConta.noMetazoa.annot.tsv* and *table_taxonomy.perUnigene.allUnigenes.tsv*.

Metadata file as well as KO id definition table must be placed in *rawdata/* directory.

* We provide the following metadata table: *metadata_metaT* as well as KO id definition table generated from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg .

Initialization file *.ini* is then create from the *MetaT.ini*:

For example:
    
    ## Enter unigene table (.tsv) (located in rawdata directory) [1]
    INPUT    main_table.mapping.unique.read_per_kb.noHuman.noConta.noMetazoa.annot.tsv

    ## Result path (the same of the composition script output) [2]
    OUTPUT    V4-unified-correct-paired-out-compo

    ## Unigene Taxonomy path (located in rawdata directory) [3]
    TAX    table_taxonomy.perUnigene.allUnigenes.tsv

    ## Database path used in the trait study (located in database directory) [4]
    DATABASE    pr2_version_4.14.0_SSU_dada2.fasta.gz

Then, run metatranscriptomic analysis script: 

* script: R_7_Metatranscriptomic_analysis.sh
* argument 1: MetaT.ini

=> `bash R_7_Metatranscriptomic_analysis.sh MetaT.ini`

    * This takes 294 min (≈5 hours) on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

## Retrieve article figures

To retrieve article figures, run following script:

* script: Retrieve_Figures.sh
* argument 1: V4-unified-correct-paired-out-compo

=> `bash Retrieve_Figures.sh V4-unified-correct-paired-out-compo`

    * This takes just a few seconds on Intel(R) Xeon(R) CPU E5-2670 with 512 Go of RAM

The resulting figures can be found in *Monjot_etal_2023/* directory.
