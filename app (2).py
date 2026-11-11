
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import math
import numpy as np

from Bio import SeqIO
from io import StringIO, BytesIO

st.set_page_config(
    page_title="Protein Sequence Logo Generator",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Protein Sequence Logo Generator")

st.write("Upload aligned protein sequences in FASTA format.")

uploaded_file = st.file_uploader(
    "Upload FASTA file",
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
            st.error("No sequences found.")

        elif any(
            not seq or set(seq) - valid_amino_acids
            for seq in sequences
        ):
            st.error("Invalid amino acid characters detected.")

        elif len(set(map(len, sequences))) != 1:
            st.error("Sequences must have equal lengths and be aligned.")

        else:

            st.success(f"{len(sequences)} sequences processed successfully!")

            amino_acids = "ACDEFGHIKLMNPQRSTVWY"

            frequency_matrix = pd.DataFrame(
                0,
                index=range(len(sequences[0])),
                columns=list(amino_acids)
            )

            for position in range(len(sequences[0])):

                for sequence in sequences:

                    aa = sequence[position]

                    frequency_matrix.loc[position, aa] += 1

            probabilities = frequency_matrix.div(
                frequency_matrix.sum(axis=1),
                axis=0
            )

            entropy = []

            for _, row in probabilities.iterrows():

                h = -sum(
                    p * math.log2(p)
                    for p in row if p > 0
                )

                entropy.append(h)

            information = [
                math.log2(20) - h
                for h in entropy
            ]

            logo_matrix = probabilities.mul(
                information,
                axis=0
            )

            st.subheader("Amino Acid Frequency")
            st.dataframe(frequency_matrix)

            st.subheader("Information Content")

            st.dataframe(pd.DataFrame({
                "Position": range(1, len(information) + 1),
                "Information Content (bits)": information
            }))

            st.subheader("Generated Protein Sequence Logo")

            fig, ax = plt.subplots(figsize=(14, 6))

            logomaker.Logo(
                logo_matrix,
                ax=ax,
                color_scheme="chemistry"
            )

            ax.set_xlabel("Amino Acid Position")
            ax.set_ylabel("Information Content (bits)")
            ax.set_title("Protein Sequence Logo")

            st.pyplot(fig)

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

    except Exception as e:

        st.error(f"Error processing FASTA file: {e}")
