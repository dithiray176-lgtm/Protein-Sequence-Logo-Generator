
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import math

from Bio import SeqIO
from io import StringIO, BytesIO

st.set_page_config(
    page_title="Protein Sequence Logo Generator",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Protein Sequence Logo Generator")

st.markdown(
    "Analyze aligned protein sequences and visualize "
    "amino acid conservation using sequence logos."
)

# Sample FASTA download

sample_fasta = """>Protein1
MKTLLLTLVV
>Protein2
MKTLLLALVV
>Protein3
MKTLLLTLVA
>Protein4
MKTLLLTLVV
"""

st.download_button(
    "Download Sample FASTA",
    data=sample_fasta,
    file_name="sample_proteins.fasta",
    mime="text/plain"
)

# FASTA upload

uploaded_file = st.file_uploader(
    "Upload your protein FASTA file",
    type=["fasta", "fa", "txt"]
)

if uploaded_file is not None:

    fasta_text = uploaded_file.getvalue().decode("utf-8")

    try:

        records = list(
            SeqIO.parse(StringIO(fasta_text), "fasta")
        )

        sequences = [
            str(record.seq).upper()
            for record in records
        ]

        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

        if not sequences:
            st.error("No FASTA sequences found.")

        elif any(
            not seq or set(seq) - valid_amino_acids
            for seq in sequences
        ):
            st.error("Invalid or empty protein sequence detected.")

        elif len(set(map(len, sequences))) != 1:
            st.error("All sequences must be aligned and have equal lengths.")

        else:

            num_sequences = len(sequences)
            sequence_length = len(sequences[0])

            st.success(
                f"Successfully processed {num_sequences} sequences!"
            )

            # Sequence statistics

            col1, col2, col3 = st.columns(3)

            col1.metric("Number of Sequences", num_sequences)
            col2.metric("Sequence Length", sequence_length)
            col3.metric("Amino Acids", 20)

            amino_acids = "ACDEFGHIKLMNPQRSTVWY"

            # Frequency calculation

            frequency_matrix = pd.DataFrame(
                0,
                index=range(sequence_length),
                columns=list(amino_acids)
            )

            for position in range(sequence_length):

                for sequence in sequences:

                    aa = sequence[position]
                    frequency_matrix.loc[position, aa] += 1

            st.subheader("Amino Acid Frequency")

            display_frequency = frequency_matrix.copy()
            display_frequency.index = range(1, sequence_length + 1)

            st.dataframe(display_frequency)

            # Probability calculation

            probabilities = frequency_matrix.div(
                frequency_matrix.sum(axis=1),
                axis=0
            )

            # Shannon entropy

            entropy = []

            for _, row in probabilities.iterrows():

                h = -sum(
                    p * math.log2(p)
                    for p in row if p > 0
                )

                entropy.append(h)

            # Information content

            information = [
                max(0, math.log2(20) - h)
                for h in entropy
            ]

            # Display information content

            st.subheader("Information Content")

            information_df = pd.DataFrame({
                "Position": range(1, sequence_length + 1),
                "Information Content (bits)": information
            })

            st.dataframe(information_df, hide_index=True)

            # Logo matrix

            logo_matrix = probabilities.mul(
                information,
                axis=0
            )

            # Sequence logo visualization

            st.subheader("Generated Protein Sequence Logo")

            fig, ax = plt.subplots(
                figsize=(max(12, sequence_length * 0.8), 6)
            )

            # Custom amino acid colours

            amino_acid_colors = {
                # Hydrophobic
                'A': '#808080',
                'V': '#808080',
                'L': '#808080',
                'I': '#808080',
                'M': '#808080',
                'F': '#808080',
                'W': '#808080',
                'P': '#808080',

                # Polar
                'S': '#2CA02C',
                'T': '#2CA02C',
                'N': '#2CA02C',
                'Q': '#2CA02C',
                'C': '#2CA02C',
                'Y': '#2CA02C',
                'G': '#2CA02C',

                # Positively charged
                'K': '#1F77B4',
                'R': '#1F77B4',
                'H': '#1F77B4',

                # Negatively charged
                'D': '#D62728',
                'E': '#D62728'
            }

            logo = logomaker.Logo(
                logo_matrix,
                ax=ax,
                color_scheme=amino_acid_colors
            )

            ax.set_xlabel("Amino Acid Position")
            ax.set_ylabel("Information Content (bits)")
            ax.set_title("Protein Sequence Logo")

            ax.set_xticks(range(sequence_length))
            ax.set_xticklabels(range(1, sequence_length + 1))

            ax.set_ylim(0, math.log2(20) + 0.2)

            ax.grid(axis="y", alpha=0.2)

            plt.tight_layout()

            st.pyplot(fig)

            # Download sequence logo

            image_buffer = BytesIO()

            fig.savefig(
                image_buffer,
                format="png",
                dpi=300,
                bbox_inches="tight"
            )

            st.download_button(
                "Download Sequence Logo",
                data=image_buffer.getvalue(),
                file_name="protein_sequence_logo.png",
                mime="image/png"
            )

            plt.close(fig)

            # Download frequency table

            st.download_button(
                "Download Frequency Table",
                data=frequency_matrix.to_csv().encode("utf-8"),
                file_name="amino_acid_frequency.csv",
                mime="text/csv"
            )

    except Exception as e:

        st.error(f"Error processing FASTA file: {e}")
