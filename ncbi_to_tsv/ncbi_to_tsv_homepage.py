import streamlit as st
import jwt  # from PyJWT

# Constants
SECRET_KEY="myrealnameisyash"

# Token verification function
def verify_token(token: str):
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=["HS256"])
        return payload
    except jwt.ExpiredSignatureError:
        st.error("Token has expired. Please log in again.")
    except jwt.InvalidTokenError:
        st.error("Invalid token. Please log in again.")
    return None

# Main function for NCBI Data Fetcher
def ncbi_data_fetcher():
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
        st.write(f"Fetching data for: {selected_page}")
        my_module = pages[selected_page]
        exec(f"import {my_module}")
        exec(f"{my_module}.main()")

# Main function with authentication
def main():
    # Grab the token from ?token=<JWT> in the URL
    params = st.experimental_get_query_params()
    token = params.get("token", [None])[0]

    # If no token, prompt the user to come from the auth app
    if not token:
        st.warning("No token found. Please log in from your main app.")
        st.stop()

    # Verify the token
    payload = verify_token(token)
    if not payload:
        # If invalid or expired, we've shown an error. Stop execution.
        st.stop()

    # If valid, proceed with the app logic
    st.write(f"Welcome, {payload.get('sub')}!")  # Display user info from JWT token

    # Call the NCBI Data Fetcher
    ncbi_data_fetcher()

if __name__ == "__main__":
    main()
