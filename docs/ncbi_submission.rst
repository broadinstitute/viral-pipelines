Submitting viral sequences to NCBI
==================================

Register your BioProject
------------------------
*If you want to add samples to an existing BioProject, skip to Step 2.*

1. Go to: https://submit.ncbi.nlm.nih.gov and login (new users - create new login).
#. Go to the Submissions tab and select BioProject - click on New Submission.
#. Follow the onscreen instructions and then click submit - you will receive a BioProject ID (``PRJNA###``) via email almost immediately.


Register your BioSamples
------------------------

1. Go to: https://submit.ncbi.nlm.nih.gov and login.
#. Go to the Submissions tab and select BioSample - click on New Submission.
#. Follow instructions, selecting "batch submission type" where applicable.
#. The metadata template to use is likely: "Pathogen affecting public health" (Pathogen.cl.1.0.xlsx).
#. Follow template instructions to fill in the sheet. Pay particular attention to the Excel comments that are attached to each column header: they describe the intended content for these columns, the valid formatting, and controlled vocabulary.

   a. For example, "organism" should always match the long name that is given by the NCBI Taxonomy database for that species.
   b. Date fields seem to have multiple acceptable formats, but we prefer ISO8601 (YYYY-MM-DD) just to reduce ambiguity.
   c. You will likely need to duplicate your sample_name to the host_subject_id column (or something like it)--if you do not, then any samples that happen to have the same attribute values will trigger an error when trying to register new BioSamples because they look like duplicates. Assuming that your sample_names are one-to-one corresponding to a human patient, host_subject_id is probably the most appropriate place to duplicate the value in order to make all entries unique.

#. Export to text and submit as .txt file. You will receive BioSamples IDs (``SAMN####``) via email (often 1-2 days later).
#. If you wish to amend/correct any metadata in your submissions, you can always do so at a future time -- however, you will need BioSample IDs before any of the following steps, so it's best to register as soon as you have collection_date and sample_name for everything. This can be a super-set of anything you submit to NCBI in the future (Genbank or SRA), so we typically register BioSamples for every viral sample we *attempt* to sequence, regardless of whether we successfully sequenced it or not.


Set up an NCBI author template
------------------------------
*If different author lists are used for different sets of samples, create a new .sbt file for each list*

1. Go to: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ 
#. Fill out the form including all authors and submitter information (if unpublished, the reference title can be just a general description of the project).
#. At the end of the form, include the BioProject number from Step 1 but NOT the BioSample number'
#. Click "create template" which will download an .sbt file to your computer'
#. Save file as "authors.sbt" or similar. If you have multiple author files, give each file a different name and prep your submissions as separate batches, one for each authors.sbt file.


Set up the BioSample map file
-----------------------------

1. Set up an Excel spreadsheet in exactly the format below:

 =========  =============
 sample     BioSample
 sample1-1  SAMNxxxxxxxxx
 sample2-1  SAMNxxxxxxxxx
 =========  =============

2. The BioSample is the BioSample number (i.e., ``SAMNxxxxxxxx``) given to you by NCBI.
3. The sample name should match the FASTA header (not necessarily the file name).

  a. Make sure your FASTA headers include segment numbers (i.e., IRF001-1) -- viral-ngs will fail otherwise! 
  b. If submitting a segmented virus (i.e., Lassa virus), each line should be a different segment, see example below (assumes sample2 is a 2-segmented virus)
  c. For samples with multiple segments, the BioSample number should be the same for all segments

     =========  =============
     sample     BioSample
     sample1-1  SAMN04488486
     sample2-1  SAMN04488657
     sample2-2  SAMN04488657
     sample3-1  SAMN04489002
     =========  =============

4. Save the file as as a tab delimited text file (e.g. "biosample-map.txt"). This file can describe *more* samples than you plan to run in a submission batch (the extras will be ignored).
5. If preparing the file on a Mac computer in Microsoft Excel (which saves tab files in a 20th-century era OS9 format), ensure that tabs and newlines are entered correctly by opening the file (via the command line) in an editor such as Nano and unchecking the [Mac-format] option (in Nano: edit the file, save the file, then click OPTION-M). You can also opt to create this file directly in a text editor, ensuring there is exactly one tab character between columns (i.e., sample<tab>BioSample in the first row). Command line converters such as ``mac2unix`` also work.


Set up the metadata file (aka Source Modifier Table)
----------------------------------------------------
1. Set up an Excel spreadsheet in exactly the format below

   a. This example shows sample2 as a 2-segmented virus.
   b. All data should be on the same line (there are 9 columns). Here they are shown as separate tables simply for space reasons.
   c. The "Sequence_ID" should match the "sample" field in the BioSample map (see Step 4). Note that this should match the FASTA header.
   d. Shown are the some of the fields we typically use in NCBI submissions, but fields can be added or removed to suit your sample needs. Other fields we often include are: "isolation_source" (i.e., serum), "collected_by" (i.e., Redeemer's University), and "genotype". Here are more details and examples provided by NCBI: https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html. A longer list of accepted column headers is provided here: https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html.
   e. The database cross-reference (db_xref) field number can be obtained by navigating to https://www.ncbi.nlm.nih.gov/taxonomy, searching for the organism of interest, and copying the "Taxonomy ID" number from the webpage.

    ===========  ===============  =======  =============================================  =====================  ==========  ============  ============  ====================================================================================
    Sequence_ID  collection_date  country  isolate                                        organism               lab_host    host          db_xref       note
    sample1-1    10-Mar-2014      Nigeria  Ebola virus/H.sapiens-tc/GIN/2014/Makona-C05   Zaire ebolavirus       Vero cells  Homo sapiens  taxon:186538  Harvest date: 01-Jan-2016; passaged 2x in cell culture (parent stock: SAMN01110234)
    sample2-1    12-Mar-2014      Nigeria  Lassa virus Macenta                            Lassa mammarenavirus   Vero cells  Homo sapiens  taxon:11620   
    sample2-2    12-Mar-2014      Nigeria  Lassa virus Macenta                            Lassa mammarenavirus   Vero cells  Homo sapiens  taxon:11620   
    sample3-1    16-Mar-2014      Nigeria  Ebola virus/H.sapiens-tc/GIN/2014/Makona-1121  Zaire ebolavirus       Vero cells  Homo sapiens  taxon:186538  This sample was collected by Dr. Blood from a very sick patient.
    ===========  ===============  =======  =============================================  =====================  ==========  ============  ============  ====================================================================================

2. The data in this table is what actually shows up on NCBI with the genome. In many cases, it is a subset of the metadata you submitted when you registered the BioSamples.
3. Save this table as sample_meta.txt. If you make the file in Excel, double check the date formatting is preserved when you save -- it should be dd-mmm-yyyy format. This file can describe *more* samples than you plan to run in a submission batch (the extras will be ignored).
4. If preparing the file on a Mac computer in Microsoft Excel (which saves tab files in a 20th-century era OS9 format), ensure that tabs and newlines are entered correctly by opening the file (via the command line) in an editor such as Nano and unchecking the [Mac-format] option (in Nano: edit the file, save the file, then click OPTION-M). You can also opt to create this file directly in a text editor, ensuring there is exactly one tab character between columns (i.e., sample<tab>BioSample in the first row). Command line converters such as ``mac2unix`` also work.


Prepare requisite input files for your submission batches
---------------------------------------------------------

1. Stage the above files you've prepared and other requisite inputs into the environment you plan to execute the :doc:`genbank` WDL workflow. If that is Terra, push these files into the appropriate GCS bucket, if DNAnexus, drop your files there. If you plan to execute locally (e.g. with ``miniwdl run``), move the files to an appropriate directory on your machine. The files you will need are the following:

   a. The files you prepared above: the submission template (authors.sbt), the biosample map (biosample-map.txt), and the source modifier table (sample_meta.txt)
   #. All of the assemblies you want to submit. These should be in fasta files, one per genome. Multi-segment/multi-chromosome genomes (such as Lassa virus, Influenza A, etc) should contain all segments within one fasta file.
   #. Your reference genome, as a fasta file. Multi-segment/multi-chromosome genomes should contain all segments within one fasta file. The fasta sequence headers should be Genbank accession numbers.
   #. Your reference gene annotations, as a series of TBL files, one per segment/chromosome. These must correspond to the accessions in you reference genome.
   #. A genome coverage table as a two-column tabular text file (optional, but helpful).
   #. The organism name (which should match what NCBI taxonomy calls the species you are submitting for). This is a string input to the workflow, not a file.
   #. The sequencing technology used. This is a string input, not a file.

#. The reference genome you provide should be annotated in the way you want your genomes annotated on NCBI. If one doesn't exist, see the addendum below about creating your own feature list.
#. Note that you will have to run the pipeline separately for each virus you are submitting AND separately for each author list.


Run the genbank submission pipeline
-----------------------------------

1. Run the :doc:`genbank` WDL workflow. Most of the metadata files described above (BioSample map, source modifier table, genome coverage table) are allowed to be a super-set of the samples you are submitting--the extra metadata will be ignored by the workflow. The samples that are included in this batch are the ones you provide to the ``assemblies_fasta`` input field. Any missing samples in the metadata inputs should not cause failures, but will produce less descriptive submission files.
#. The :doc:`genbank` workflow performs the following steps: it aligns your assemblies against a Genbank reference sequence, transfers gene annotation from that Genbank reference into your assemblies' coordinate spaces, and then takes your genomes, the transferred annotations, and all of the sample metadata prepared above, and produces a zipped bundle that you send to NCBI. There are two zip bundles: ``sequins_only.zip`` is the file to email to NCBI. ``all_files.zip`` contains a full set of files for your inspection prior to submission.
#. In the ``all_files.zip`` output, for each sample, you will see a ``.sqn``, ``.gbf``, ``.val``, and ``.tbl`` file. You should also see an ``errorsummary.val`` file that you can use to check for annotation errors (or you can check the ``.val`` file for each sample individually). Ideally, your samples should be error-free before you submit them to NCBI unless you're confident enough in the genomic evidence for unusual coding effects and frameshifts. For an explanation of the cryptic error messages, see: https://www.ncbi.nlm.nih.gov/genbank/genome_validation/.
#. We currently use a bioconda wrapper of NCBI's `tbl2asn` tool called `tbl2asn-forever`. This works around some deficiencies in NCBI's tool but has the side effect of setting the submission date to Jan 1, 2019 for all submission, regardless of today's date. Unfortunately, until NCBI releases a fixed tool, you will need to search-replace the date in the SQN files in a text editor prior to submission.
#. Check your ``.gbf`` files for a preview of what your genbank entries will look like. Once you are happy with your files email the ``sequins_only.zip`` file to gb-sub@ncbi.nlm.nih.gov.
#. It often takes 2-8 weeks to receive a response and accession numbers for your samples. Do follow up if you haven’t heard anything for a few weeks!
