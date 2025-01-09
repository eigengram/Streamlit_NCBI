import streamlit as st
import xml.etree.ElementTree as ET
from Bio import Entrez
import pandas as pd

def fetch_sra_details(accession):
    handle = Entrez.efetch(db="sra", id=accession, retmode="xml")
    xml_data = handle.read()
    handle.close()
    return xml_data

def read_accession_numbers(file):
    accession_numbers = [line.decode('utf-8').strip() for line in file.readlines()]
    return accession_numbers

def main():
    # Streamlit app
    st.title("SRA Details Fetcher")
    st.write("Upload a text file containing SRA accession numbers to fetch their details.")

    # File uploader
    uploaded_file = st.file_uploader("Choose a file...", type="txt")

    if uploaded_file is not None:
        accession_numbers = read_accession_numbers(uploaded_file)

        st.write(f"{len(accession_numbers)} accession numbers found")
        display_accession_nums = ", ".join(accession_numbers)
        st.markdown(f"<p><i>{display_accession_nums}</i></p>", unsafe_allow_html=True)

        st.divider()
        st.write("Fetching data")

        progress_bar = st.progress(0)
        count_acc_nums = len(accession_numbers)
        progress = 0

        # Initialize a list to hold all SRA info from all accessions
        all_sra_info = []

        for accession_number in accession_numbers:
            xml_data = fetch_sra_details(accession_number)
            root = ET.fromstring(xml_data)
            sra_info = []

            for experiment_package in root.findall('.//EXPERIMENT_PACKAGE'):
                experiment_id = experiment_package.find('.//EXPERIMENT').get('accession', "N/A")
                title = experiment_package.find('.//TITLE').text if experiment_package.find('.//TITLE') is not None else "N/A"
                study_title = experiment_package.find('.//STUDY_TITLE').text if experiment_package.find('.//STUDY_TITLE') is not None else "N/A"
                sample_name = experiment_package.find('.//SAMPLE').find('.//SCIENTIFIC_NAME').text if experiment_package.find('.//SAMPLE/SCIENTIFIC_NAME') is not None else "N/A"

                sra_info.append({
                    'Experiment ID': experiment_id,
                    'Title': title,
                    'Study Title': study_title,
                    'Sample Name': sample_name
                })

            all_sra_info.extend(sra_info)
            progress += 100//count_acc_nums
            progress_bar.progress(progress)
        progress_bar.progress(100)

        # Convert the list of all nucleotide info to DataFrame
        df_sra_info = pd.DataFrame(all_sra_info)

        # Display the data
        st.write("### SRA Details")
        st.dataframe(df_sra_info)

# Run the app
if __name__ == "__main__":
    main()
