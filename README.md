# MitoSonar
12S metabarcoding-based Fish Taxonomy Identifier

## How to Run

### **Step one:** Make sure the working directory is set up accondingly

    Use this tree as reference:

    work-dir/
        ├── data/ <- this folder will store the generated OTU and taxonomy tables and the .fna of filtered sequences for blasting
        └── data-raw/
            ├── MiFish_all_mitogenomes.fasta <- fasta for creating blast database
            └── fastqs/ <- sequences to analyse must be placed into this folder

NOTE: the blast output table will be created under the 'data-raw' folder

### **Step two:** Install R and RStudio into your machine
You can download R here: [The R Project for Statistical Computing](https://www.r-project.org)
and RStudio here: [RStudio Desktop](https://posit.co/download/rstudio-desktop/)

### **Step three:** Make sure to install NCBI BLAST+ to your system, since the script will ask the OS to run blast commands rather than executing them inside R.
You can learn how to do so in this official tutorial by NCBI: [BLAST® Command Line Applications User Manual](https://www.ncbi.nlm.nih.gov/books/NBK569861/)

### **Step four:** Open the script MitoSonar_script.r in RStudio
In the first portion of the script, it will guide you through some package installation steps.
Be aware that you might need to install some dependencies **manually**!
*I must say I kinda feel sorry for you, but just be patient and it should work fine. (Be warned: **dada2** can be pretty tricky to set up in the first use)*

After installing dependencies, set 'work-dir' as working directy and voila! **Time to fish**