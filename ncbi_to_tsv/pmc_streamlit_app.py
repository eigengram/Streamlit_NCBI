import streamlit as st
import xml.etree.ElementTree as ET
from Bio import Entrez
import pandas as pd

def fetch_pmc_details(accession):
    handle = Entrez.efetch(db="pmc", id=accession, retmode="xml")
    xml_data = handle.read()
    handle.close()
    return xml_data

def read_accession_numbers(file):
    accession_numbers = [line.decode('utf-8').strip() for line in file.readlines()]
    return accession_numbers

def main():
    # Streamlit app
    st.title("PMC Details Fetcher")
    st.write("Upload a text file containing PMC accession numbers to fetch their details.")

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

        # Initialize a list to hold all PMC info from all accessions
        all_pmc_info = []

        for accession_number in accession_numbers:
            xml_data = fetch_pmc_details(accession_number)
            root = ET.fromstring(xml_data)
            pmc_info = []

            for article in root.findall('.//article'):
                front = article.find('.//front')
                article_title = front.find('.//article-title').text if article.find('.//article-title') is not None else "N/A"
                journal_title = front.find('.//journal-title').text if article.find('.//journal-title') is not None else "N/A"

                journal_id = []
                for journalID in front.findall(".//journal-id"):
                    if journalID.get('journal-id-type') is not None:
                        journal_id.append({journalID.get('journal-id-type'): journalID.text if journalID.text is not None else "N/A"})
                
                article_id = []
                for articleID in front.findall(".//article-id"):
                    if articleID.get('pub-id-type') is not None:
                        article_id.append({articleID.get('pub-id-type'): articleID.text if articleID.text is not None else "N/A"})

                pub_dates = []
                for pub_date in front.findall('.//pub-date'):
                    if pub_date.get('pub-type') is not None:
                        pub_dates.append({pub_date.get('pub-type'): pub_date.find('.//year').text if pub_date.find('.//year') is not None else "N/A"})
                
                contributors = []
                for contributor in front.findall('.//contrib'):
                    if contributor.get('contrib-type') is not None:
                        name = {}
                        if contributor.find('./name/surname') is not None:
                            name['surname'] = contributor.find('./name/surname').text if contributor.find('./name/surname').text is not None else "N/A"
                        if contributor.find('./name/given-names') is not None:
                            name['given_names'] = contributor.find('./name/given-names').text if contributor.find('./name/given-names').text is not None else "N/A"
                        contributors.append({contributor.get('contrib-type'): name})

                pmc_info.append({
                    'PMC ID': accession_number,
                    'Article Title': article_title,
                    'Journal Title': journal_title,
                    'Publication Year': pub_dates,
                    'Journal_id': journal_id,
                    'Aritcle_id': article_id,
                    'Contributors': contributors
                })

            all_pmc_info.extend(pmc_info)
            progress += 100//count_acc_nums
            progress_bar.progress(progress)
        progress_bar.progress(100)

        # Convert the list of all nucleotide info to DataFrame
        df_pmc_info = pd.DataFrame(all_pmc_info)

        # Display the data
        st.write("### PMC Details")
        st.dataframe(df_pmc_info)

# Run the app
if __name__ == "__main__":
    main()