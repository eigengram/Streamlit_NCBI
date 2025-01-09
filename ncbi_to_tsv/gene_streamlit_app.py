import streamlit as st
from Bio import Entrez
import pandas as pd

def fetch_gene_details(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
    records = Entrez.read(handle)
    return records

def read_accession_numbers(file):
    accession_numbers = [line.decode('utf-8').strip() for line in file.readlines()]
    return accession_numbers

def main():
    # Streamlit app
    st.title("Gene Details Fetcher")
    st.write("Upload a text file containing Gene accession numbers to fetch their details.")

    # File uploader
    uploaded_file = st.file_uploader("Choose a file...", type="txt")

    if uploaded_file is not None:
        accession_numbers = read_accession_numbers(uploaded_file)

        st.write(f"{len(accession_numbers)} accession numbers found")
        display_accession_nums = ", ".join(accession_numbers)
        st.markdown(f"<p><i>{display_accession_nums}</i></p>", unsafe_allow_html=True)

        st.divider()
        st.write("Fetching data")

        # Initialize lists to hold extracted data
        all_gene_info = []

        progress_bar = st.progress(0)
        count_acc_nums = len(accession_numbers)
        progress = 0

        for accession_number in accession_numbers:
            gene_details = fetch_gene_details(accession_number)
            for record in gene_details:
                gene_id = record.get('GBSeq_primary-accession', "N/A")
                gene_length = record.get('GBSeq_length', "N/A")
                keywords = record.get('GBSeq_keywords')
                create_date = record.get('GBSeq_create-date')
                update_date = record.get('GBSeq_update-date', "N/A")
                definition = record.get('GBSeq_definition', "N/A")    
                locus = record.get('GBSeq_locus')
                strandedness = record.get('GBSeq_strandedness')
                moltype = record.get('GBSeq_moltype')
                topology = record.get('GBSeq_topology')
                division = record.get('GBSeq_division')
                accession_version = record.get('GBSeq_accession-version')
                other_seqids = record.get('GBSeq_other-seqids')
                project = record.get('GBSeq_project')
                source = record.get('GBSeq_source')
                organism = record.get('GBSeq_organism', "N/A")
                taxonomy = record.get('GBSeq_taxonomy', "N/A")
                references = record.get('GBSeq_references')
                comments = record.get('GBSeq_comment', "N/A")
                contig = record.get('GBSeq_contig')
                xrefs = record.get('GBSeq_xrefs')

                # More detailed features can be found in the features section
                features = record.get('GBSeq_feature-table', [])

                gene_details_list = []
                for feature in features:
                    feature_key = feature['GBFeature_key']
                    feature_location = feature.get('GBFeature_location', "N/A")
                    qualifiers = feature.get('GBFeature_quals', [])
                    feature_info = {'Key': feature_key, 'Location': feature_location}

                    for qual in qualifiers:
                        feature_info[qual['GBQualifier_name'].capitalize()] = qual.get('GBQualifier_value', 'N/A')

                    gene_details_list.append(feature_info)

                all_gene_info.append({
                    'Gene_ID': gene_id,
                    'Gene_length': gene_length,
                    'Keywords': keywords,
                    'Definition': definition,
                    'Create_date': create_date,
                    'Update_date': update_date,
                    'Locus': locus,
                    'Strandedness': strandedness,
                    'Moltype': moltype,
                    'Topology': topology,
                    'Division': division,
                    'Accession_version': accession_version,
                    'Other_sequids': other_seqids,
                    'Project': project,
                    'Source': source,
                    'Organism': organism,
                    'Taxonomy': taxonomy,
                    'References': references,
                    'Comments': comments,
                    'Contig': contig,
                    'Xrefs': xrefs,
                    'Features': gene_details_list
                })
            progress += 100//count_acc_nums
            progress_bar.progress(progress)
        progress_bar.progress(100)

        # Convert the list of all Gene info to DataFrame
        df_gene_info = pd.DataFrame(all_gene_info)

        # Display the data
        st.write("### Gene Details")
        st.dataframe(df_gene_info)

# Run the app
if __name__ == "__main__":
    main()