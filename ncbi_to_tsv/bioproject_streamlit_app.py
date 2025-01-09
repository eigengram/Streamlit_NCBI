import streamlit as st
import xml.etree.ElementTree as ET
from Bio import Entrez
import pandas as pd

# Set your email here
# Entrez.email = "your.email@example.com"

# Function to fetch BioProject details
def fetch_bioproject_details(accession):
    handle = Entrez.efetch(db="bioproject", id=accession, retmode="xml")
    xml_data = handle.read()
    handle.close()
    return xml_data

# Function to read accession numbers from a file
def read_accession_numbers(file):
    accession_numbers = [line.decode('utf-8').strip() for line in file.readlines()]
    return accession_numbers

def main():
    # Streamlit app
    st.title("BioProject Details Fetcher")
    st.write("Upload a text file containing BioProject accession numbers to fetch their details.")

    # File uploader
    uploaded_file = st.file_uploader("Choose a file...", type="txt")

    if uploaded_file is not None:
        accession_numbers = read_accession_numbers(uploaded_file)

        st.write(f"{len(accession_numbers)} accession numbers found")
        display_accession_nums = ", ".join(accession_numbers)
        st.markdown(f"<p><i>{display_accession_nums}</i></p>", unsafe_allow_html=True)

        st.divider()
        st.write("Fetching data")

        # Initialize a list to hold all BioProject info from all accessions
        all_bioproject_info = []

        progress_bar = st.progress(0)
        count_acc_nums = len(accession_numbers)
        progress = 0

        for accession_number in accession_numbers:
            xml_data = fetch_bioproject_details(accession_number)
            root = ET.fromstring(xml_data)

            for project_summary in root.findall('.//DocumentSummary'):
                uid = project_summary.get('uid', "N/A")

                archive = project_summary.find('.//ArchiveID')
                if archive is not None:
                    accession = archive.get('accession', "N/A")
                    archive_id = archive.get('id', 'N/A')
                else:
                    accession = "N/A"
                    archive_id = "N/A"

                project_desc = project_summary.find('.//ProjectDescr')
                if project_desc is not None:
                    name = project_desc.find('.//Name').text if project_desc.find('.//Name') is not None else "N/A"
                    title = project_desc.find('.//Title').text if project_desc.find('.//Title') is not None else "N/A"
                    description = project_desc.find('.//Description').text if project_desc.find('.//Description') is not None else "N/A"
                    release_date = project_desc.find('.//ProjectReleaseDate').text if project_desc.find('.//ProjectReleaseDate') is not None else "N/A"
                else:
                    name = "N/A"
                    title = "N/A"
                    description = "N/A"
                    release_date = "N/A"

                publications = []
                if project_desc is not None:
                    for publication in project_desc.findall('.//Publication'):
                        pub_id = publication.get('id', "N/A")
                        pub_status = publication.get('status', "N/A")
                        pub_reference = publication.find('.//Reference').text if publication.find('.//Reference') is not None else "N/A"

                        pub_title = publication.find('.//StructuredCitation/Title').text if publication.find('.//StructuredCitation/Title') is not None else "N/A"
                        journal_title = publication.find('.//StructuredCitation/Journal/JournalTitle').text if publication.find('.//StructuredCitation/Journal/JournalTitle') is not None else "N/A"
                        journal_year = publication.find('.//StructuredCitation/Journal/Year').text if publication.find('.//StructuredCitation/Journal/Year') is not None else "N/A"
                        
                        author_list = []
                        for author in publication.findall('.//Author'):
                            author_name = {}
                            author_name['First'] = author.find('.//First').text if author.find('.//First') is not None else 'N/A'
                            author_name['Last'] = author.find('.//Last').text if author.find('.//Last') is not None else 'N/A'
                            author_list.append(author_name)
                        
                        publications.append({
                            'ID': pub_id,
                            'Status': pub_status,
                            'Reference': pub_reference,
                            'Title': pub_title,
                            'Journal_Title': journal_title,
                            'Journal_year': journal_year,
                            'Authors': author_list
                        })

                all_bioproject_info.append({
                    'UID': uid,
                    'Accession_number': accession,
                    'Archive_id': archive_id,
                    'Name': name,
                    'Title': title,
                    'Description': description,
                    'Release_date': release_date,
                    'Publications': publications
                })

            progress += 100//count_acc_nums
            progress_bar.progress(progress)
        progress_bar.progress(100)

        # Convert the list of all BioProject info to DataFrame
        df_bioproject_info = pd.DataFrame(all_bioproject_info)

        # Display the data
        st.write("### BioProject Details")
        st.dataframe(df_bioproject_info)

# Run the app
if __name__ == "__main__":
    main()