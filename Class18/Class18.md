Class 18
================

## Investigating Cancer Genomics Datasets

# The GenomicDataCommons R Package

The workflow will be:

Install packages if not already installed Load libraries Identify and
download somatic variants for a representative TCGA dataset, in this
case pancreatic adenocarcinoma. Use maftools to provide rich summaries
of the data.

``` r
library(GenomicDataCommons)
```

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'GenomicDataCommons'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(TCGAbiolinks)
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

    ## Registered S3 method overwritten by 'R.oo':
    ##   method        from       
    ##   throw.default R.methodsS3

``` r
library(maftools)
```

Let’s check the status of ‘GenomicDataCommons’ or GDC.

``` r
status()
```

    ## $commit
    ## [1] "e588f035feefee17f562b3a1bc2816c49a2b2b19"
    ## 
    ## $data_release
    ## [1] "Data Release 16.0 - March 26, 2019"
    ## 
    ## $status
    ## [1] "OK"
    ## 
    ## $tag
    ## [1] "1.20.0"
    ## 
    ## $version
    ## [1] 1

# Querying the GDC from R

We will typically start our interaction with the GDC by searching the
resource to find data that we are interested in investigating further.
In GDC speak this is called “Querying GDC metadata”. Metadata here
refers to the extra descriptive information associated with the actual
patient data (i.e. ‘cases’) in the GDC.

The are four main sets of metadata that we can query, namely projects(),
cases(), files(), and annotations(). We will start with projects()

``` r
projects <- getGDCprojects()
```

If you use the View(projects) function call you can see all the project
names (such as Neuroblastoma, Pancreatic Adenocarcinoma, etc.) along
with their project IDs (such as TARGET-NBL, TCGA-PAAD, etc.) and
associated information.

``` r
head(projects)
```

    ##   dbgap_accession_number
    ## 1                   <NA>
    ## 2              phs000466
    ## 3                   <NA>
    ## 4                   <NA>
    ## 5              phs001444
    ## 6              phs000471
    ##                                                          disease_type
    ## 1 Cystic, Mucinous and Serous Neoplasms, Adenomas and Adenocarcinomas
    ## 2                                    Clear Cell Sarcoma of the Kidney
    ## 3                                               Mesothelial Neoplasms
    ## 4                                        Adenomas and Adenocarcinomas
    ## 5                     Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
    ## 6                                               High-Risk Wilms Tumor
    ##   releasable released state
    ## 1      FALSE     TRUE  open
    ## 2      FALSE     TRUE  open
    ## 3      FALSE     TRUE  open
    ## 4      FALSE     TRUE  open
    ## 5      FALSE     TRUE  open
    ## 6      FALSE     TRUE  open
    ##                                                                                     primary_site
    ## 1 Rectosigmoid junction, Unknown, Rectum, Colon, Connective, subcutaneous and other soft tissues
    ## 2                                                                                         Kidney
    ## 3                                              Heart, mediastinum, and pleura, Bronchus and lung
    ## 4   Other and unspecified parts of biliary tract, Gallbladder, Liver and intrahepatic bile ducts
    ## 5                                                                                    Lymph Nodes
    ## 6                                                                                         Kidney
    ##     project_id           id
    ## 1    TCGA-READ    TCGA-READ
    ## 2  TARGET-CCSK  TARGET-CCSK
    ## 3    TCGA-MESO    TCGA-MESO
    ## 4    TCGA-CHOL    TCGA-CHOL
    ## 5 NCICCR-DLBCL NCICCR-DLBCL
    ## 6    TARGET-WT    TARGET-WT
    ##                                                  name tumor
    ## 1                               Rectum Adenocarcinoma  READ
    ## 2                    Clear Cell Sarcoma of the Kidney  CCSK
    ## 3                                        Mesothelioma  MESO
    ## 4                                  Cholangiocarcinoma  CHOL
    ## 5 Genomic Variation in Diffuse Large B Cell Lymphomas DLBCL
    ## 6                               High-Risk Wilms Tumor    WT

Moving onto cases() we can use an example from the package associated
publication to answer our first from question above (i.e. find the
number of cases/patients across different projects within the GDC):

``` r
cases_by_project <- cases() %>%
  facet("project.project_id") %>%
  aggregations()
head(cases_by_project$project.project_id)
```

    ##          key doc_count
    ## 1      FM-AD     18004
    ## 2 TARGET-NBL      1127
    ## 3  TCGA-BRCA      1098
    ## 4 TARGET-AML       988
    ## 5  TARGET-WT       652
    ## 6   TCGA-GBM       617

Note that the facet() and aggregations() functions here are from the
GenomicDataCommons package and act to group all cases by the project id
and then count them up.

Write the R code to make a barplot of the cases per project. Lets plot
this data with a log scale for the y axis (log=“y”), rotated axis labels
(las=2) and color the bar coresponding to the TCGA-PAAD project.

``` r
x <- cases_by_project$project.project_id

# Make a custom color vector for our plot
colvec <- rep("lightblue", nrow(x))
colvec[x$key == "TCGA-PAAD"] <- "red"

# Plot with 'log' for y axis and rotate labels with 'las'
par(mar = c(9, 4, 2, 2))  
barplot(x$doc_count, names.arg=x$key, log="y", col=colvec, las=2)
```

![](Class18_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Lets explore some other functions from the related TCGAbiolinks package.

We can use the getSampleFilesSummary() function to determine for a given
project how many cases and what type of data we have available for each
case:

``` r
samp <- getSampleFilesSummary("TCGA-PAAD")
```

    ## Accessing information for project: TCGA-PAAD

    ## Using 'state_comment' as value column. Use 'value.var' to override

    ## Aggregation function missing: defaulting to length

``` r
head(samp)
```

    ##            .id Biospecimen_Biospecimen Supplement
    ## 1 TCGA-2J-AAB1                                 14
    ## 2 TCGA-2J-AAB4                                 14
    ## 3 TCGA-2J-AAB6                                 14
    ## 4 TCGA-2J-AAB8                                 14
    ## 5 TCGA-2J-AAB9                                 14
    ## 6 TCGA-2J-AABA                                 14
    ##   Biospecimen_Slide Image_Diagnostic Slide
    ## 1                                        1
    ## 2                                        1
    ## 3                                        1
    ## 4                                        1
    ## 5                                        1
    ## 6                                        1
    ##   Biospecimen_Slide Image_Tissue Slide Clinical_Clinical Supplement
    ## 1                                    1                            8
    ## 2                                    1                            8
    ## 3                                    1                            8
    ## 4                                    1                            8
    ## 5                                    1                            8
    ## 6                                    1                            8
    ##   Copy Number Variation_Copy Number Segment_Genotyping Array_Affymetrix SNP 6.0
    ## 1                                                                             2
    ## 2                                                                             2
    ## 3                                                                             2
    ## 4                                                                             2
    ## 5                                                                             2
    ## 6                                                                             2
    ##   Copy Number Variation_Gene Level Copy Number Scores_Genotyping Array_Affymetrix SNP 6.0
    ## 1                                                                                       1
    ## 2                                                                                       1
    ## 3                                                                                       1
    ## 4                                                                                       1
    ## 5                                                                                       1
    ## 6                                                                                       1
    ##   Copy Number Variation_Masked Copy Number Segment_Genotyping Array_Affymetrix SNP 6.0
    ## 1                                                                                    2
    ## 2                                                                                    2
    ## 3                                                                                    2
    ## 4                                                                                    2
    ## 5                                                                                    2
    ## 6                                                                                    2
    ##   DNA Methylation_Methylation Beta Value_Methylation Array_Illumina Human Methylation 450
    ## 1                                                                                       1
    ## 2                                                                                       1
    ## 3                                                                                       1
    ## 4                                                                                       1
    ## 5                                                                                       1
    ## 6                                                                                       1
    ##   Sequencing Reads_Aligned Reads_miRNA-Seq_Illumina
    ## 1                                                 1
    ## 2                                                 1
    ## 3                                                 1
    ## 4                                                 1
    ## 5                                                 1
    ## 6                                                 1
    ##   Sequencing Reads_Aligned Reads_RNA-Seq_Illumina
    ## 1                                               1
    ## 2                                               1
    ## 3                                               1
    ## 4                                               1
    ## 5                                               1
    ## 6                                               1
    ##   Sequencing Reads_Aligned Reads_WXS_Illumina
    ## 1                                           2
    ## 2                                           2
    ## 3                                           2
    ## 4                                           2
    ## 5                                           2
    ## 6                                           2
    ##   Simple Nucleotide Variation_Aggregated Somatic Mutation_WXS
    ## 1                                                           4
    ## 2                                                           4
    ## 3                                                           4
    ## 4                                                           4
    ## 5                                                           4
    ## 6                                                           4
    ##   Simple Nucleotide Variation_Annotated Somatic Mutation_WXS
    ## 1                                                          4
    ## 2                                                          4
    ## 3                                                          4
    ## 4                                                          4
    ## 5                                                          4
    ## 6                                                          4
    ##   Simple Nucleotide Variation_Masked Somatic Mutation_WXS
    ## 1                                                       4
    ## 2                                                       4
    ## 3                                                       4
    ## 4                                                       4
    ## 5                                                       4
    ## 6                                                       4
    ##   Simple Nucleotide Variation_Raw Simple Somatic Mutation_WXS
    ## 1                                                           4
    ## 2                                                           4
    ## 3                                                           4
    ## 4                                                           4
    ## 5                                                           4
    ## 6                                                           4
    ##   Transcriptome Profiling_Gene Expression Quantification_RNA-Seq
    ## 1                                                              3
    ## 2                                                              3
    ## 3                                                              3
    ## 4                                                              3
    ## 5                                                              3
    ## 6                                                              3
    ##   Transcriptome Profiling_Isoform Expression Quantification_miRNA-Seq
    ## 1                                                                   1
    ## 2                                                                   1
    ## 3                                                                   1
    ## 4                                                                   1
    ## 5                                                                   1
    ## 6                                                                   1
    ##   Transcriptome Profiling_miRNA Expression Quantification_miRNA-Seq
    ## 1                                                                 1
    ## 2                                                                 1
    ## 3                                                                 1
    ## 4                                                                 1
    ## 5                                                                 1
    ## 6                                                                 1
    ##     project
    ## 1 TCGA-PAAD
    ## 2 TCGA-PAAD
    ## 3 TCGA-PAAD
    ## 4 TCGA-PAAD
    ## 5 TCGA-PAAD
    ## 6 TCGA-PAAD

Now we can use GDCquery() function to focus in on a particular data type
that we are interested in. For example, to answer our second question
from above - namely ‘find all gene expression data files for all
pancreatic cancer patients’:

``` r
query <- GDCquery(project = "TCGA-PAAD",
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification")
```

    ## --------------------------------------

    ## o GDCquery: Searching in GDC database

    ## --------------------------------------

    ## Genome of reference: hg38

    ## --------------------------------------------

    ## oo Accessing GDC. This might take a while...

    ## --------------------------------------------

    ## ooo Project: TCGA-PAAD

    ## --------------------

    ## oo Filtering results

    ## --------------------

    ## ooo By data.type

    ## ----------------

    ## oo Checking data

    ## ----------------

    ## ooo Check if there are duplicated cases

    ## Warning: There are more than one file for the same case. Please verify query results. You can use the command View(getResults(query)) in rstudio

    ## ooo Check if there results for the query

    ## -------------------

    ## o Preparing output

    ## -------------------

``` r
ans <- getResults(query)
```

How to get a list of the ‘data.category’ argument for each project:

``` r
data.category <- TCGAbiolinks:::getProjectSummary("TCGA-PAAD")
data.category$data_categories
```

    ##   case_count file_count               data_category
    ## 1        178        912     Transcriptome Profiling
    ## 2        185        737       Copy Number Variation
    ## 3        183       1480 Simple Nucleotide Variation
    ## 4        184        195             DNA Methylation
    ## 5        185        211                    Clinical
    ## 6        185        740            Sequencing Reads
    ## 7        185       1032                 Biospecimen

``` r
head(ans)
```

    ##   data_release                      data_type
    ## 1  12.0 - 16.0 Gene Expression Quantification
    ## 2  12.0 - 16.0 Gene Expression Quantification
    ## 3  12.0 - 16.0 Gene Expression Quantification
    ## 4  12.0 - 16.0 Gene Expression Quantification
    ## 5  12.0 - 16.0 Gene Expression Quantification
    ## 6  12.0 - 16.0 Gene Expression Quantification
    ##                   updated_datetime
    ## 1 2018-11-30T10:59:48.252246+00:00
    ## 2 2018-11-30T10:59:48.252246+00:00
    ## 3 2018-11-30T10:59:48.252246+00:00
    ## 4 2018-11-30T10:59:48.252246+00:00
    ## 5 2018-11-30T10:59:48.252246+00:00
    ## 6 2018-11-30T10:59:48.252246+00:00
    ##                                              file_name
    ## 1  27f41c98-3658-48c9-a257-d9ecad3276c8.FPKM-UQ.txt.gz
    ## 2  2132e4b3-f882-4d05-9cb6-be9a24c510c0.FPKM-UQ.txt.gz
    ## 3 25b906db-fd54-4b85-b67e-e421826bd794.htseq.counts.gz
    ## 4     1b54f883-d25d-4f6c-9398-7d7dcfb75653.FPKM.txt.gz
    ## 5     7051e52f-069d-48d5-966e-064a01bf2725.FPKM.txt.gz
    ## 6     6c73911d-8d0a-4a5a-9251-4ded7ea70fef.FPKM.txt.gz
    ##                                  submitter_id
    ## 1 27f41c98-3658-48c9-a257-d9ecad3276c8_uqfpkm
    ## 2 2132e4b3-f882-4d05-9cb6-be9a24c510c0_uqfpkm
    ## 3  25b906db-fd54-4b85-b67e-e421826bd794_count
    ## 4   1b54f883-d25d-4f6c-9398-7d7dcfb75653_fpkm
    ## 5   7051e52f-069d-48d5-966e-064a01bf2725_fpkm
    ## 6   6c73911d-8d0a-4a5a-9251-4ded7ea70fef_fpkm
    ##                                file_id file_size
    ## 1 1e0be26b-da47-4e56-926c-9afc259a2bd0    528102
    ## 2 2803abb9-bf21-4b75-a3c9-687c87b75701    496856
    ## 3 4d98a852-c30a-437b-a214-ee81015d674b    250164
    ## 4 af9e4765-462e-4fa3-819c-f890456647ba    494993
    ## 5 62809aca-c5ac-47f1-bfed-b82d811a797d    504209
    ## 6 e8ae7104-605a-4bd3-b58d-b85da7f4fcc3    525367
    ##                          cases                                   id
    ## 1 TCGA-FB-AAPS-01A-12R-A39D-07 1e0be26b-da47-4e56-926c-9afc259a2bd0
    ## 2 TCGA-FB-AAQ3-01A-11R-A41B-07 2803abb9-bf21-4b75-a3c9-687c87b75701
    ## 3 TCGA-HZ-8005-01A-11R-2204-07 4d98a852-c30a-437b-a214-ee81015d674b
    ## 4 TCGA-FB-AAPQ-01A-11R-A41B-07 af9e4765-462e-4fa3-819c-f890456647ba
    ## 5 TCGA-HV-A7OL-01A-11R-A33R-07 62809aca-c5ac-47f1-bfed-b82d811a797d
    ## 6 TCGA-3A-A9IB-01A-21R-A39D-07 e8ae7104-605a-4bd3-b58d-b85da7f4fcc3
    ##                   created_datetime                           md5sum
    ## 1 2016-05-26T21:20:54.657488-05:00 23b40da735dd52e9867927e253b7b135
    ## 2 2016-05-29T10:32:57.442916-05:00 94e1c8501c2b0120c5ca238dcaffb500
    ## 3 2016-05-26T21:21:27.312488-05:00 38acc858f2e4390e10c82b9423ca4f3c
    ## 4 2016-05-30T18:57:37.259818-05:00 2736865e0da315e5e5c72efd2f013b4a
    ## 5 2016-05-29T10:07:29.441428-05:00 4387ffbcfc574513b3910fb8d1ac24d6
    ## 6 2016-05-29T10:15:47.187935-05:00 86c8f9940f4dd39c899df7c17f341d29
    ##   data_format access    state version           data_category
    ## 1         TXT   open released       1 Transcriptome Profiling
    ## 2         TXT   open released       1 Transcriptome Profiling
    ## 3         TXT   open released       1 Transcriptome Profiling
    ## 4         TXT   open released       1 Transcriptome Profiling
    ## 5         TXT   open released       1 Transcriptome Profiling
    ## 6         TXT   open released       1 Transcriptome Profiling
    ##              type experimental_strategy   project
    ## 1 gene_expression               RNA-Seq TCGA-PAAD
    ## 2 gene_expression               RNA-Seq TCGA-PAAD
    ## 3 gene_expression               RNA-Seq TCGA-PAAD
    ## 4 gene_expression               RNA-Seq TCGA-PAAD
    ## 5 gene_expression               RNA-Seq TCGA-PAAD
    ## 6 gene_expression               RNA-Seq TCGA-PAAD
    ##                            analysis_id        analysis_updated_datetime
    ## 1 593c4dd1-0fde-4a26-9d37-0d8dfb2c3693 2018-09-10T15:08:41.786316-05:00
    ## 2 5c56e496-5fa6-4880-a297-0339d7519a7c 2018-09-10T15:08:41.786316-05:00
    ## 3 52be80a4-7268-4caa-aaff-d449a13a3a41 2018-09-10T15:08:41.786316-05:00
    ## 4 a51cf7ed-5f4a-486f-9b28-b62ece59c185 2018-09-10T15:08:41.786316-05:00
    ## 5 c1473452-e376-413a-ad15-9b9e3efe94cd 2018-09-10T15:08:41.786316-05:00
    ## 6 1abb4868-1793-46c4-bccd-90c0934cc1b2 2018-09-10T15:08:41.786316-05:00
    ##          analysis_created_datetime
    ## 1 2016-05-26T21:20:54.657488-05:00
    ## 2 2016-05-29T10:32:57.442916-05:00
    ## 3 2016-05-26T21:21:27.312488-05:00
    ## 4 2016-05-30T18:57:37.259818-05:00
    ## 5 2016-05-29T10:07:29.441428-05:00
    ## 6 2016-05-29T10:15:47.187935-05:00
    ##                         analysis_submitter_id analysis_state
    ## 1 27f41c98-3658-48c9-a257-d9ecad3276c8_uqfpkm       released
    ## 2 2132e4b3-f882-4d05-9cb6-be9a24c510c0_uqfpkm       released
    ## 3  25b906db-fd54-4b85-b67e-e421826bd794_count       released
    ## 4   1b54f883-d25d-4f6c-9398-7d7dcfb75653_fpkm       released
    ## 5   7051e52f-069d-48d5-966e-064a01bf2725_fpkm       released
    ## 6   6c73911d-8d0a-4a5a-9251-4ded7ea70fef_fpkm       released
    ##                 analysis_workflow_link analysis_workflow_type
    ## 1 https://github.com/NCI-GDC/htseq-cwl        HTSeq - FPKM-UQ
    ## 2 https://github.com/NCI-GDC/htseq-cwl        HTSeq - FPKM-UQ
    ## 3 https://github.com/NCI-GDC/htseq-cwl         HTSeq - Counts
    ## 4 https://github.com/NCI-GDC/htseq-cwl           HTSeq - FPKM
    ## 5 https://github.com/NCI-GDC/htseq-cwl           HTSeq - FPKM
    ## 6 https://github.com/NCI-GDC/htseq-cwl           HTSeq - FPKM
    ##   analysis_workflow_version   tissue.definition
    ## 1                        v1 Primary solid Tumor
    ## 2                        v1 Primary solid Tumor
    ## 3                        v1 Primary solid Tumor
    ## 4                        v1 Primary solid Tumor
    ## 5                        v1 Primary solid Tumor
    ## 6                        v1 Primary solid Tumor

``` r
nrow(ans)
```

    ## [1] 546

## Designing a personalized cancer vaccine

# Protein sequences from healthy and tumor tissue

``` r
library(bio3d)
```

Read in the fasta sequence files

``` r
seq <- read.fasta("lecture18_sequences.fa")
seq
```

    ##              1        .         .         .         .         .         60 
    ## P53_wt       MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
    ## P53_mutant   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP
    ##              **************************************** ******************* 
    ##              1        .         .         .         .         .         60 
    ## 
    ##             61        .         .         .         .         .         120 
    ## P53_wt       DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ## P53_mutant   DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ##              **** ******************************************************* 
    ##             61        .         .         .         .         .         120 
    ## 
    ##            121        .         .         .         .         .         180 
    ## P53_wt       SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ## P53_mutant   SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ##              ************************************************************ 
    ##            121        .         .         .         .         .         180 
    ## 
    ##            181        .         .         .         .         .         240 
    ## P53_wt       RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ## P53_mutant   RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ##              ******************************** *************************** 
    ##            181        .         .         .         .         .         240 
    ## 
    ##            241        .         .         .         .         .         300 
    ## P53_wt       SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
    ## P53_mutant   SCMGGMNRRPILTIITLEV-----------------------------------------
    ##              ******************                                           
    ##            241        .         .         .         .         .         300 
    ## 
    ##            301        .         .         .         .         .         360 
    ## P53_wt       PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
    ## P53_mutant   ------------------------------------------------------------
    ##                                                                           
    ##            301        .         .         .         .         .         360 
    ## 
    ##            361        .         .         .  393 
    ## P53_wt       GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    ## P53_mutant   ---------------------------------
    ##                                                
    ##            361        .         .         .  393 
    ## 
    ## Call:
    ##   read.fasta(file = "lecture18_sequences.fa")
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   2 sequence rows; 393 position columns (259 non-gap, 134 gap) 
    ## 
    ## + attr: id, ali, call

``` r
mismatches <- which(seq$ali[1,] != seq$ali[2,])
```

Or… use a function in
    bio3d

``` r
conserv(seq, method = "identity")
```

    ##   [1] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ##  [18] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ##  [35] 1.0 1.0 1.0 1.0 1.0 1.0 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ##  [52] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.5 1.0 1.0 1.0
    ##  [69] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ##  [86] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [103] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [120] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [137] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [154] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [171] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [188] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [205] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [222] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [239] 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    ## [256] 1.0 1.0 1.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [273] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [290] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [307] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [324] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [341] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [358] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [375] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    ## [392] 0.5 0.5

There are a lot of mismatches at the end. We don’t want the gaps.

``` r
gaps <- gap.inspect(seq)
gaps <- gaps$t.inds
```

We can subtract the gaps from the mismatches

``` r
tumor.spec <- mismatches[!mismatches %in% gaps]
```

Let’s see the mutation
sites

``` r
ids <- paste(seq$ali[1,tumor.spec], tumor.spec, seq$ali[2,tumor.spec], sep = "")
ids
```

    ## [1] "D41L"  "R65W"  "R213V" "D259V"

Find 9-mers with mismatches.

``` r
start.ind <- tumor.spec - 8
end.ind <- tumor.spec + 8
  
tumor.seq <- matrix("-", nrow=length(ids), ncol=17)
rownames(tumor.seq) <- ids

for( i in 1:length(tumor.spec) ) {
  tumor.seq[i,] <- seq$ali[2,start.ind[i]:end.ind[i]]
}

tumor.seq
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ## D41L  "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "L"  "D"   "L"   "M"   "L"  
    ## R65W  "D"  "P"  "G"  "P"  "D"  "E"  "A"  "P"  "W"  "M"   "P"   "E"   "A"  
    ## R213V "Y"  "L"  "D"  "D"  "R"  "N"  "T"  "F"  "V"  "H"   "S"   "V"   "V"  
    ## D259V "I"  "L"  "T"  "I"  "I"  "T"  "L"  "E"  "V"  "-"   "-"   "-"   "-"  
    ##       [,14] [,15] [,16] [,17]
    ## D41L  "S"   "P"   "D"   "D"  
    ## R65W  "A"   "P"   "P"   "V"  
    ## R213V "V"   "P"   "Y"   "E"  
    ## D259V "-"   "-"   "-"   "-"

Blank out the empty sections

``` r
tumor.seq[tumor.seq == "-"] <- ""
```

Write this to a fasta file for processing

``` r
write.fasta(ids = ids, seqs = tumor.seq, file = "subsequences.fa")
```
