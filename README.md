# Streamlit Application

This guide explains how to set up and run the Streamlit application, as well as how to install the required dependencies.

---

## Prerequisites

Make sure you have the following installed on your system:

- **Python 3.8 or later**
- **pip** (Python package manager)

---

## Installation

### Step 1: Clone the Repository

Clone the project repository to your local machine:

```bash
git clone https://github.com/your-repository/streamlit-app.git
cd streamlit-app
```

### Step 2: Create a Virtual Environment (Optional)

It is recommended to create a virtual environment to manage dependencies:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### Step 3: Install Dependencies

Install the required libraries using the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

---

## Running the Application

1. Navigate to the directory where the `app.py` file is located.

2. Run the Streamlit application using the following command:

   ```bash
   streamlit run app.py
   ```

3. Open your web browser and navigate to:

   ```
   http://localhost:8501
   ```

   The application should now be running!

---

## Additional Notes

- Ensure your dataset and any required files (e.g., `metadata.csv`, `document_embeddings.npy`) are placed in the correct directory as specified in your code.
- If you encounter any issues, double-check the versions of the libraries in the `requirements.txt` file.

---

## Requirements

The dependencies required for this application are listed in `requirements.txt`. Make sure to install them using the command provided in the **Installation** section.

---

## JWT Token Verification in Streamlit

### 1. Requirements for JWT Verification

- **Streamlit** (`pip install streamlit`)
- **PyJWT** (`pip install pyjwt`)

### 2. Shared Secret Key

Your authentication system (e.g., Next.js) signs tokens with a **secret key** or private key. You need the matching **secret/public key** in your Streamlit app to verify the tokens.

Example:

```python
SECRET_KEY = "MY_SHARED_SECRET"
```

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.
