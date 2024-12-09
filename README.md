# CORT_Temp_learning

#### Research paper

This repository contains the final code, data, and figures used in the following manuscript:

Recio et al. 2024. 

Supplementary Materials: The supplementary materials associated with this paper is integrated within the ms.docx or ms.qmd files. 

#### How to use this repository?

Users can download a zip file of the entire repository by clicking on the green code tab at the top of the page and then clicking Download ZIP. Users who already have a GitHub account can fork the repository.

The key file in this repository is the 📄 ms.qmd. This file can be rendered in R with Quarto to reproduce the entire paper. Code chunks within the file provide the code used to reproduce figures and analyses along with supporting statements within the text. Note that inline code chunks use specific objects which are then rendered.

The 📄 ms.qmd file makes use of files within a number of folders that are identified in the code chunks. There are a number of important folders in the repository.

  📂 data folder contains all the raw data used in files. Note that there are different files, but the main one employed in our analyses is Learning.csv. In 
  📂 output/figs/ Folder contains all the figures for the paper that are read and included in the paper. See more details below (3. Figures).
  📂 R The R folder contains three files that are used to clean and process data to prepare it for use in the 📄 ms.qmd file. Note that readers do not need to open and run these files, but they are simply here to document the workflow and code used to clean up data to be used. These include:
        📄 1_data_processing.R, which is used to first download the Google Sheets used during our hackathon, process and merge with meta-data from the California Digital Library (CDL);
        📄 2_data_cleaning.R File does some additional cleaning and checking of the data and merges disparate datasets together.
        📄 3_author_affli.R, which is code to grab and process author affiliations because the lead author is a little lazy when it comes to cumbersome tasks such as these.
    📂 bib The bib folder contains:
        📄 refs.bib the references;
        📄 proceedings-of-the-royal-society-b.csl the journal formatting style file;
        📄 template.docx a template docx file to format the resulting rendered files.

3. Figures

📄 ms.qmd will rely on figures generated and stored in the 📂 output/figs/ folder. The final list of figures are as follows:

    📄 Figure 1_FINAL.png
    📄 Figure 2_FINAL.png

Note that these figures are a composite of figures patched together and the final files were modified outside of R for aesthetic reasons. These files also have associated Adobe illustrator files. The code to reproduce individual figures is provided in 📄 ms.qmd. When rendered individual files will be written to the 📂 output/figs/ folder.
4. Data

Given the project has been a major group effort we initially relied on Google Sheets to provide pathways by which all authors could contribute to the data collection process. These Google Sheets were then sourced, processed, checked and then cleaned prior to analysis. The first initial file was provided to us by the California Digital Library (CDL) team, who downloaded the relevant meta-data on the articles posted to EcoEvoRxiv as of 30 September 2023.

The original meta-data files from the CDL are located in the 📂 data folder, but the main one used is:

    📄 20231003_EER_preprints_metadata.xlsx

The column names here are:

    Preprint ID: Janeway's internal identifier for the preprint
    Preprint Title: Title of the preprint
    Preprint DOI: DOI of the preprint
    Publisher DOI: DOI of the postprint/publisher's article, if any
    Reuse Licence: Creative Commons reuse licence
    Submission Date: Date preprint was submitted to EcoEvoRxiv
    Accepted Date: Date preprint was accepted to EER
    Published Date: Date preprint was published in EER (may differ from accepted date)
    Update Date: Date preprint was last updated by an EER moderator
    Current Version: Current version now
    Version creation date: Date that version was created/submitted (may differ from update date)
    Submitting Author: Name of submitting author
    Submitting Author Email: Submitting author's email address
    Authors List: List of all authors
    Total authors: Total number of authors

The second data file used in the 📄 ms.qmd file is:

    📄 final_data2.csv

This file is the processed data file which is then merged with the 📄 20231003_EER_preprints_metadata.xlsx data. The descriptions of the data columns are as follows:

    assigned_to: Name of author assigned to checking the data.
    nr: Unique article number for each preprint/postprint.
    extractors_first_name: First name of the author who extracted data.
    extractors_last_name: Last name of the author who extracted data.
    preprint_id: Unique preprint ID from EcoEvoRxiv. These are auto-assigned numbers by the Janeway system.
    preprint_doi: Digital Object Identifier (DOI) for preprint.
    preprint_title: Title of pre/postprint
    submitting_author: Full name of submitted author.
    submitting_author_country: Country of submitting author.
    submitting_author_first_publication_year: Year of first publication for submitting author.
    taxa_being_studied: Taxa being studied in the preprint.
    data_link_preprint: Link to associated data for the paper on EcoEvoRxiv.
    code_link_preprint: Link to associated code for the paper on EcoEvoRxiv.
    number_of_citations_preprint: Number of citations to the preprint taken from Google Scholar.
    pci_recommendation_preprint: Whether or not the preprint was posted to Peer Community In (PCI).
    publication_doi: DOI for the published article.
    publication_journal: Research journal name that the article was ultimately published in.
    publication_date: Date of publication for the published version of the article. Taken as the date the article first appeared online.
    publication_title_changed: Whether or not the title of the published articles changed from the preprint.
    number_of_citations_article: Number of citations to the published article.
    current_version: The current version of the article published on EcoEvoRxiv
    data_link_article: Link to associated data for the published version of the paper.
    code_link_article: Link to associated code for the published version of the paper.
    publisher_doi: DOI for the published paper.
    preprint_published_date: Preprint/postprint publication date. The date when published on EcoEvoRxiv.
    total_authors: Total number of authors on the preprint/postprint on EcoEvoRxiv.
    time_between_preprint_and_pub_days: Days between posting the article on EcoEvoRxiv and the time the paper (if relevant) was published in a peer-reviewer research journal.
    postprint: Identifier as to whether the article posted on EcoEvoRxiv was likely a postprint or not.
    plants: Identifier on whether the study contained plants as one of the main taxa studied or not.
    algi: Identifier on whether the study contained algi as one of the main taxa studied or not.
    fungi: Identifier on whether the study contained fungi as one of the main taxa studied or not.
    microorganisms: Identifier on whether the study contained bacteria, or other microorganisms, as one of the main taxa studied or not.
    invertebrates: Identifier on whether the study contained invertebrates as one of the main taxa studied or not.
    vertebrates: Identifier on whether the study contained vertebrates as one of the main taxa studied or not.
    type_of_preprint: The type of article on EcoEvoRxiv. These could be research article, opinion, review and meta-analysis, book or book chapters, report. Anything not in these categories were classified as 'other'.
    data_link_preprint_cleaned: A cleaned up version of the data link for the article on EcoEvoRxiv.
    code_link_preprint_cleaned A cleaned up version of the code link for the article on EcoEvoRxiv.
    data_link_article_cleaned: A cleaned up version of the data link for the published article.
    code_link_article_cleaned A cleaned up version of the code link for the published article.
    check_completed: Indicator (yes/no) as to whether the paper check was complete.
    check_comment: Comments about decisions made during checking.
    journal_name: Name of research journal where the preprint/postprint was published.
    is_oa: Logical (true/false) indicating whether the paper was published in open access (true) or not (false).
    oa_status: Open access status.

5. Reporting Issues or Asking Questions
If anything is unclear or you require further detail please don't hesitate to lodge an issue.
