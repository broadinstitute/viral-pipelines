Submitting viral sequences to NCBI
==================================

Step 1: Register your BioProject
--------------------------------
*If you want to add samples to an existing BioProject, skip to Step 2.*

1. Go to: https://submit.ncbi.nlm.nih.gov and login (new users - create new login).
1. Go to the Submissions tab and select BioProject - click on New Submission.
1. Follow the onscreen instructions and then click submit - you will receive a BioProject ID (``PRJNA###``) via email almost immediately.


Step 2: Register your BioSamples
--------------------------------

1. Go to: https://submit.ncbi.nlm.nih.gov and login.
1. Go to the Submissions tab and select BioSample - click on New Submission.
1. Follow instructions, selecting "batch submission type" where applicable.
1. The metadata template to use is likely: "Pathogen affecting public health".
1. Follow template instructions (careful about date formatting) and submit as .txt file.
1. You will receive BioSamples IDs (``SAMN####``) via email (often 1-2 days later).

Step 3: Set up an NCBI author template
--------------------------------------
*If different author lists are used for different sets of samples, create a new .sbt file for each list*

1. Go to: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ 
1. Fill out the form including all authors and submitter information (if unpublished, the reference title can be just a general description of the project).
1. At the end of the form, include the BioProject number from Step 1 but NOT the BioSample number'
1. Click "create template" which will download an .sbt file to your computer'
1. Save file as "authors.sbt" or similar. If you have multiple author files, give each file a different name and prep your submissions as separate batches, one for each authors.sbt file.

Step 4: Set up the BioSample map file
-------------------------------------

1. Set up an Excel spreadsheet in exactly the format below:
 ---------  -------------
 sample     BioSample
 sample1-1  SAMNxxxxxxxxx
 sample2-1  SAMNxxxxxxxxx
 ---------  -------------
1. The BioSample is the BioSample number (i.e., ``SAMNxxxxxxxx``) given to you by NCBI.
1. The sample name should match the FASTA header (not necessarily the file name).
  a. Make sure your FASTA headers include segment numbers (i.e., IRF001-1) -- viral-ngs will fail otherwise! 
  a. If submitting a segmented virus (i.e., Lassa virus), each line should be a different segment, see example below (assumes sample2 is a 2-segmented virus)
  a. For samples with multiple segments, the BioSample number should be the same for all segments
     ---------  ------------
     sample     BioSample
     sample1-1  SAMN04488486
     sample2-1  SAMN04488657
     sample2-2  SAMN04488657
     sample3-1  SAMN04489002
     ---------  ------------
1. Save the file as as a tab delimited text file (e.g. "biosample-map.txt").
1. If preparing the file on a Mac computer in Microsoft Excel (which saves tab files in a 20th-century era OS9 format), ensure that tabs and newlines are entered correctly by opening the file (via the command line) in an editor such as Nano and unchecking the [Mac-format] option (in Nano: edit the file, save the file, then click OPTION-M). You can also opt to create this file directly in a text editor, ensuring there is exactly one tab character between columns (i.e., sample<tab>BioSample in the first row). Command line converters such as ``mac2unix`` also work.

Step 5: Set up the metadata file (aka Source Modifier Table)
------------------------------------------------------------
1. Set up an Excel spreadsheet in exactly the format below
  a. This example shows sample2 as a 2-segmented virus.
  a. All data should be on the same line (there are 9 columns). Here they are shown as separate tables simply for space reasons.
  a. The "Sequence_ID" should match the "sample" field in the BioSample map (see Step 4). Note that this should match the FASTA header.
  a. Shown are the some of the fields we typically use in NCBI submissions, but fields can be added or removed to suit your sample needs. Other fields we often include are: "isolation_source" (i.e., serum), "collected_by" (i.e., Redeemer's University), and "genotype". Here is the full list of fields accepted by NCBI: https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html.
  a. The database cross-reference (db_xref) field number can be obtained by navigating to https://www.ncbi.nlm.nih.gov/taxonomy, searching for the organism of interest, and copying the "Taxonomy ID" number from the webpage.
  ------------ ---------------- -------- --------------------------------------------- --------------------- ---------- ------------ ------------ -----------------------------------------------------------------------------------
  Sequence_ID  collection_date  country  isolate                                       organism              lab_host   host         db_xref      note
  sample1-1    10-Mar-2014      Nigeria  Ebola virus/H.sapiens-tc/GIN/2014/Makona-C05  Zaire ebolavirus      Vero cells Homo sapiens taxon:186538 Harvest date: 01-Jan-2016; passaged 2x in cell culture (parent stock: SAMN01110234)
  sample2-1    12-Mar-2014      Nigeria  Lassa virus Macenta                           Lassa mammarenavirus  Vero cells Homo sapiens taxon:11620  
  sample2-2    12-Mar-2014      Nigeria  Lassa virus Macenta                           Lassa mammarenavirus  Vero cells Homo sapiens taxon:11620  
  sample3-1    16-Mar-2014      Nigeria  Ebola virus/H.sapiens-tc/GIN/2014/Makona-1121 Zaire ebolavirus      Vero cells Homo sapiens taxon:186538 This sample was collected by Dr. Blood from a very sick patient.
  ------------ ---------------- -------- --------------------------------------------- --------------------- ---------- ------------ ------------ -----------------------------------------------------------------------------------
1. The data in this table is what actually shows up on NCBI with the genome. In many cases, it is a subset of the metadata you submitted when you registered the BioSamples.
1. Save this table as sample_meta.txt. If you make the file in Excel, double check the date formatting is preserved when you save -- it should be dd-mmm-yyyy format.
1. If preparing the file on a Mac computer in Microsoft Excel (which saves tab files in a 20th-century era OS9 format), ensure that tabs and newlines are entered correctly by opening the file (via the command line) in an editor such as Nano and unchecking the [Mac-format] option (in Nano: edit the file, save the file, then click OPTION-M). You can also opt to create this file directly in a text editor, ensuring there is exactly one tab character between columns (i.e., sample<tab>BioSample in the first row). Command line converters such as ``mac2unix`` also work.

Step 6: Prepare requisite input files for your submission batches
-----------------------------------------------------------------

1. TO DO -- more description here (authors.sbt file, your biosample-map.txt file, and your sample_meta.txt)
1. The reference genome you provide should be annotated in the way you want your genomes annotated on NCBI. If one doesn't exist, see the addendum below about creating your own feature list.
1. Note that you will have to run the pipeline separately for each virus you are submitting AND separately for each author list.
Copy your 


Step 7: Run the genbank submission pipeline
-------------------------------------------

1. TO DO -- more description here -- <genbank>
1. For each sample, you will see a .sqn, .gbf, .val, and .tbl file. You should also see an errorsummary.val file that you can use to check for annotation errors (or you can check the .val file for each sample individually). Ideally, your samples should be error-free before you submit them to NCBI. For an explanation of the cryptic error messages, see: https://www.ncbi.nlm.nih.gov/genbank/genome_validation/.
1. TO DO -- moltype?
1. Check your .gbf files for a preview of what your genbank entries will look like. Once you are happy with your files, zip up all of the .sqn files (for all of the samples you are submitting, regardless of author list or organism) and email the .zip file to gb-sub@ncbi.nlm.nih.gov.
1. It often takes 2-8 weeks to receive a response and accession numbers for your samples. Do follow up if you haven’t heard anything for a few weeks!

