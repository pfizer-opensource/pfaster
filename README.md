# PfaSTer
pneumococcal fasta serotyping 
- *S. pneumoniae* in-silico serotype prediction from assembled genome sequences (.fasta)
- Developed by Pfizer.inc

## Installation
Clone the repository and enter the directory:

    git clone https://github.com/PfizerRD/pfaster.git
    cd pfaster

Build and activate the environment:

    conda env create -f environment.yml
    conda activate serotyper

## Usage

To predict the serotype from a genome sequence:

    python pfaster.py -f path_to_fasta_file -o output_file_directory

The -o input is optional. Results will only print to the console if no directory is specified.


### Tests

Sample genomes are available for testing. For instance:

    python pfaster.py -f tests/Pn1_test_ERR1439829.fasta

## Additional
Test genomes are sourced from the European Nucleotide Archive (ENA, https://www.ebi.ac.uk/ena)

For more information, contact jonathan.lee@pfizer.com

