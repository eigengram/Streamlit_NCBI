import streamlit as st

def main():
    st.title("NCBI Data Fetcher")

    st.write("Select the type of data you want to fetch:")

    # Create a dictionary with page names and their corresponding module names
    pages = {
        "BioProject": "ncbi_to_tsv.bioproject_streamlit_app",
        "Proteins": "ncbi_to_tsv.protein_streamlit_app",
        "Gene": "ncbi_to_tsv.gene_streamlit_app",
        "PMC": "ncbi_to_tsv.pmc_streamlit_app",
        "Nucleotide": "ncbi_to_tsv.nucleotide_streamlit_app",
        "SRA": "ncbi_to_tsv.sra_streamlit_app"
    }

    # Selectbox for navigation
    options = ["Select"] + list(pages.keys())
    selected_page = st.selectbox("Choose type of data", options)

    # Load the selected page
    if selected_page and selected_page != "Select":
        st.write(selected_page)
        my_module = pages[selected_page]
        exec(f"import {my_module}")
        exec(f"{my_module}.main()")

# Run the app
if __name__ == "__main__":
    main()