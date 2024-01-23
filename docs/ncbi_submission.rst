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
   b. Date fields seem to have multiple acceptable formats, but we prefer ISO8601 (YYYY-MM-DD) just to reduce ambiguity. Note that this format will trigger a warning when uploading, if you don't have HH:MM time values as well (it will suggest an edit for you).
   c. You will likely need to duplicate your sample_name to the host_subject_id column (or something like it)--if you do not, then any samples that happen to have the same attribute values will trigger an error when trying to register new BioSamples because they look like duplicates. Assuming that your sample_names are one-to-one corresponding to a human patient, host_subject_id is probably the most appropriate place to duplicate the value in order to make all entries unique.
   d. Populate the isolate column using the naming convention you want to apply to this organism (most viral species have a specific, structured naming convention you should follow). Our workflow will re-use this value for the Genbank record name.

#. Export to text and submit as .txt file. You will receive BioSamples IDs (``SAMN####``) via email (often 1-2 days later).
#. After NCBI accepts your submission and registers your samples, retrieve the text-formatted "attributes table" associated with this submission from the portal at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on "Download attributes file with BioSample accessions". You will need this file later. Do not use the file that was attached to the NCBI response email--it does not contain the full record and is formatted differently.
#. If you wish to amend/correct any metadata in your submissions, you can always do so at a future time -- however, you will need BioSample IDs before any of the following steps, so it's best to register as soon as you have collection_date and sample_name for everything. This can be a super-set of anything you submit to NCBI in the future (Genbank or SRA), so we typically register BioSamples for every viral sample we *attempt* to sequence, regardless of whether we successfully sequenced it or not.


Set up an NCBI author template
------------------------------
*If different author lists are used for different sets of samples, create a new .sbt file for each list*

1. Go to: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ 
#. Fill out the form including all authors and submitter information (if unpublished, the reference title can be just a general description of the project).
#. At the end of the form, include the BioProject number from Step 1 but NOT the BioSample number'
#. Click "create template" which will download an .sbt file to your computer'
#. Save file as "authors.sbt" or similar. If you have multiple author files, give each file a different name and prep your submissions as separate batches, one for each authors.sbt file.


Prepare requisite input files for your submission batches
---------------------------------------------------------

1. Stage the above files you've prepared and other requisite inputs into the environment you plan to execute the :doc:`genbank` WDL workflow. If that is Terra, push these files into the appropriate GCS bucket, if DNAnexus, drop your files there. If you plan to execute locally (e.g. with ``miniwdl run``), move the files to an appropriate directory on your machine. The files you will need are the following:

   a. The files you prepared above: the submission template (authors.sbt) and the biosample attributes table (attributes.tsv).
   #. All of the assemblies you want to submit. These should be in fasta files, one per genome. Multi-segment/multi-chromosome genomes (such as Lassa virus, Influenza A, etc) should contain all segments within one fasta file.
   #. Your reference genome, as a set of fasta files, one per segment/chromosome. The fasta sequence headers should be Genbank accession numbers. This can come directly from Genbank.
   #. Your reference gene annotations, as a set of TBL files, one per segment/chromosome. These must correspond to the accessions in your reference genome. These must be presented in the same order as the reference genome fasta files, which must also be in the same order as all the sequences in all of your assembled fasta files.
   #. A genome coverage table as a two-column tabular text file (optional, but helpful).
   #. The organism name (which should match what NCBI taxonomy calls the species you are submitting for). This is a string input to the workflow, not a file.
   #. The sequencing technology used. This is a string input, not a file.
   #. The NCBI Taxonomy taxid. This is a numeric input, not a file.

#. The reference genome you provide should be annotated in the way you want your genomes annotated on NCBI. If one doesn't exist, see the addendum below about creating your own feature list.
#. Note that you will have to run the pipeline separately for each virus you are submitting AND separately for each author list.


Run the genbank submission pipeline
-----------------------------------

1. Run the :doc:`genbank` WDL workflow. Most of the metadata files described above (BioSample map, source modifier table, genome coverage table) are allowed to be a super-set of the samples you are submitting--the extra metadata will be ignored by the workflow. The samples that are included in this batch are the ones you provide to the ``assemblies_fasta`` input field. Any missing samples in the metadata inputs should not cause failures, but will produce less descriptive submission files. Viral genomes should set molType to ``cRNA``.
#. The :doc:`genbank` workflow performs the following steps: it aligns your assemblies against a Genbank reference sequence, transfers gene annotation from that Genbank reference into your assemblies' coordinate spaces, and then takes your genomes, the transferred annotations, and all of the sample metadata prepared above, and produces a zipped bundle that you send to NCBI. There are two zip bundles: ``sequins_only.zip`` is the file to email to NCBI. ``all_files.zip`` contains a full set of files for your inspection prior to submission.
#. In the ``all_files.zip`` output, for each sample, you will see a ``.sqn``, ``.gbf``, ``.val``, and ``.tbl`` file. You should also see an ``errorsummary.val`` file that you can use to check for annotation errors (or you can check the ``.val`` file for each sample individually). Ideally, your samples should be error-free before you submit them to NCBI unless you're confident enough in the genomic evidence for unusual coding effects and frameshifts. For an explanation of the cryptic error messages, see: https://www.ncbi.nlm.nih.gov/genbank/genome_validation/.
#. We currently use a bioconda wrapper of NCBI's `tbl2asn` tool called `tbl2asn-forever`. This works around some deficiencies in NCBI's tool but has the side effect of setting the submission date to Jan 1, 2019 for all submission, regardless of today's date. Unfortunately, until NCBI releases a fixed tool, you will need to search-replace the date in the SQN files in a text editor prior to submission.
#. Check your ``.gbf`` files for a preview of what your genbank entries will look like. Once you are happy with your files email the ``sequins_only.zip`` file to gb-sub@ncbi.nlm.nih.gov.
#. It often takes 2-8 weeks to receive a response and accession numbers for your samples. Do follow up if you haven’t heard anything for a few weeks!
