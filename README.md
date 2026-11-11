# Protein-Sequence-Logo-Generator
Protein sequence analysis and visualization tool using Biopython
# 🧬 Protein Sequence Logo Generator

A web-based bioinformatics application that analyzes aligned protein sequences and generates sequence logos to visualize amino acid frequency, conservation and information content at each position.

## 📌 Project Overview

The Protein Sequence Logo Generator is designed to simplify the analysis and visualization of multiple aligned protein sequences.

It accepts protein sequences in FASTA format, processes them using Biopython, calculates position-wise amino acid frequencies and information content, and generates a graphical sequence logo.

## 🎯 Objectives

- Parse protein sequences from FASTA files.
- Validate amino acid sequences.
- Calculate amino acid frequencies at each alignment position.
- Convert frequencies into probabilities.
- Calculate Shannon entropy and information content.
- Generate protein sequence logos.
- Allow users to download generated visualizations and frequency tables.

## ⚙️ Technologies Used

| Technology | Purpose |
|---|---|
| Python | Backend programming |
| Biopython | FASTA parsing and sequence processing |
| Streamlit | Web application interface |
| Pandas | Frequency tables and data processing |
| NumPy | Numerical computation |
| Matplotlib | Graph generation |
| Logomaker | Sequence logo visualization |

## 🧪 Methodology

### 1. FASTA Input and Parsing
Users upload aligned protein sequences in FASTA format. Biopython's `SeqIO` module parses the uploaded file and extracts the protein sequences.

### 2. Sequence Validation
The application checks:
- Whether sequences are present.
- Whether amino acid characters are valid.
- Whether all sequences have equal lengths.

### 3. Amino Acid Frequency Calculation
The occurrence of each amino acid is counted at every position across all sequences.

### 4. Probability Calculation
Frequencies are converted into probabilities:

\[
P_i(a)=\frac{f_i(a)}{N}
\]

Where:
- \(f_i(a)\) = frequency of amino acid a at position i.
- \(N\) = total number of sequences.

### 5. Shannon Entropy
The variability at each position is calculated using:

\[
H_i=-\sum_a P_i(a)\log_2 P_i(a)
\]

### 6. Information Content
Information content is calculated using a uniform background model:

\[
IC_i=\log_2(20)-H_i
\]

The maximum information content is approximately 4.3219 bits.

### 7. Sequence Logo Generation
The probability matrix is multiplied by the information content at each position. Logomaker uses this matrix to generate a graphical representation of amino acid conservation.

## ✨ Features

- FASTA file upload.
- Sample FASTA file download.
- Protein sequence validation.
- Amino acid frequency tables.
- Information content calculation.
- Custom amino acid colour scheme.
- Protein sequence logo visualization.
- Downloadable PNG sequence logo.
- Downloadable CSV frequency table.

## 🚀 Installation and Execution

Clone the repository:

```bash
git clone YOUR_GITHUB_REPOSITORY_URL
cd YOUR_REPOSITORY_NAME
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the application:

```bash
streamlit run app.py
```

## 📂 Project Structure

```text
Protein-Sequence-Logo-Generator/
│
├── app.py
├── requirements.txt
└── README.md
```

## 👩‍💻 Project Team

**Group 27**

- Aniruddh Chaturvedi
- Dithi Ray

## 📊 Limitations

- Input sequences must be pre-aligned and have equal lengths.
- Only the 20 standard amino acids are currently accepted.
- Information content uses a uniform amino acid background model.
- Small-sample correction is not currently implemented.

## 🔮 Future Scope

- Support for additional background frequency models.
- Small-sample correction for information content.
- Improved sequence alignment integration.
- Additional sequence analysis tools.
- Enhanced visualization and export options.

---

**Project:** Protein Sequence Logo Generator  
**Domain:** Bioinformatics  
**Application:** Python and Streamlit
