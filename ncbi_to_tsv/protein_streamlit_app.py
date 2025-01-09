import streamlit as st
from Bio import Entrez
import pandas as pd

def fetch_protein_details(accession):
    handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="xml")
    records = Entrez.read(handle)
    return records

def read_accession_numbers(file):
    accession_numbers = [line.decode('utf-8').strip() for line in file.readlines()]
    return accession_numbers

def main():
    # Streamlit app
    st.title("Protein Details Fetcher")
    st.write("Upload a text file containing Protein accession numbers to fetch their details.")

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

        # Initialize a list to hold all protein info from all accessions
        all_protein_info = []

        for accession_number in accession_numbers:
            protein_details = fetch_protein_details(accession_number)
            protein_info = []
            
            for record in protein_details:

                protein_id = record.get('GBSeq_primary-accession', "N/A")
                protein_length = record.get('GBSeq_length', "N/A")
                keywords = record.get('GBSeq_keywords')
                create_date = record.get('GBSeq_create-date')
                update_date = record.get('GBSeq_update-date', "N/A")
                definition = record.get('GBSeq_definition', "N/A")    
                locus = record.get('GBSeq_locus')
                moltype = record.get('GBSeq_moltype')
                topology = record.get('GBSeq_topology')
                division = record.get('GBSeq_division')
                accession_version = record.get('GBSeq_accession-version')
                other_seqids = record.get('GBSeq_other-seqids')
                source = record.get('GBSeq_source')
                organism = record.get('GBSeq_organism', "N/A")
                taxonomy = record.get('GBSeq_taxonomy', "N/A")
                references = record.get('GBSeq_references')
                comments = record.get('GBSeq_comment', "N/A")
                source_db = record.get('GBSeq_source_db', "N/A")
                sequence = record.get('GBSeq_sequence')

                # More detailed features can be found in the features section
                features = record.get('GBSeq_feature-table', [])

                protein_details_list = []
                for feature in features:
                    feature_key = feature['GBFeature_key']
                    feature_location = feature.get('GBFeature_location', "N/A")
                    qualifiers = feature.get('GBFeature_quals', [])
                    feature_info = {'Key': feature_key, 'Location': feature_location}

                    for qual in qualifiers:
                        feature_info[qual['GBQualifier_name'].capitalize()] = qual.get('GBQualifier_value', 'N/A')

                    protein_details_list.append(feature_info)

                protein_info.append({
                    'Protein_ID': protein_id,
                    'Protein_length': protein_length,
                    'Keywords': keywords,
                    'Definition': definition,
                    'Create_date': create_date,
                    'Update_date': update_date,
                    'Locus': locus,
                    'Moltype': moltype,
                    'Topology': topology,
                    'Division': division,
                    'Accession_version': accession_version,
                    'Other_sequids': other_seqids,
                    'Source': source,
                    'Organism': organism,
                    'Taxonomy': taxonomy,
                    'References': references,
                    'Comments': comments,
                    'Features': protein_details_list,
                    'Source_db': source_db,
                    'Sequence': sequence[:50] + '...',
                })

            all_protein_info.extend(protein_info)
            progress += 100//count_acc_nums
            progress_bar.progress(progress)
        progress_bar.progress(100)

        # Convert the list of all nucleotide info to DataFrame
        df_protein_info = pd.DataFrame(all_protein_info)

        # Display the data
        st.write("### Protein Details")
        st.dataframe(df_protein_info)

# Run the app
if __name__ == "__main__":
    main()
