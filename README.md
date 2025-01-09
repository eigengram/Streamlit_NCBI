# NCBI Data Fetcher

## Overview

The NCBI Data Fetcher is a suite of Streamlit-based applications designed to fetch and display data from various NCBI databases, including BioProject, Gene, Protein, PMC, SRA, and Nucleotide. Users can upload text files containing accession numbers, and the applications will fetch the corresponding data and display it in a user-friendly format.

## Features

- Fetches data from the following NCBI databases:
  - **BioProject**: Retrieves BioProject details using accession numbers.
  - **Gene**: Fetches gene details such as ID, length, taxonomy, and features.
  - **Protein**: Retrieves protein-specific details like sequence, source, and organism information.
  - **PMC**: Extracts publication details from the PMC database.
  - **Nucleotide**: Fetches nucleotide sequence information.
  - **SRA**: Retrieves sequencing experiment details.
- Supports bulk data fetching using a text file with multiple accession numbers.
- Displays results in a tabular format for easy analysis.

## Prerequisites

- Python 3.8 or higher
- A valid email registered with NCBI to access the Entrez API (uncomment `Entrez.email` and set your email in the respective scripts).

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/eigengram/Streamlit_NCBI.git
   cd Streamlit_NCBI

Install the required dependencies:pip install -r requirements.txt

Run all the applications using command streamlit run {name of file}.py