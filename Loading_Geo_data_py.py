# %%
import os
import pandas as pd
import tarfile
import gzip
from pathlib import Path
import requests
from io import StringIO, BytesIO
from IPython.display import display
import matplotlib.pyplot as plt

class GEODataAccess:
    def __init__(self, base_dir=None, series_id="GSE171485"):
        """
        Initialize GEO data access object.
            base_dir (str): Local directory containing GEO files (None for remote access)
            series_id (str): GEO Series ID (e.g., "GSE171485")
        """
        self.base_dir = base_dir
        self.series_id = series_id

    def get_file_path(self, filename):
        """Get full path to a file"""
        if self.base_dir:
            return os.path.join(self.base_dir, filename)
        else:
            raise ValueError("Base directory is not set for local file access.")

    def load_file(self, filename):
        """Load a file based on its type"""
        file_path = self.get_file_path(filename)
        if filename.endswith('.soft'):
            return self._read_soft_file(file_path)
        elif filename.endswith('.tar'):
            return self._read_tar_file(file_path)
        elif filename.endswith('.tsv'):
            return self._read_tsv_file(file_path)
        elif filename.endswith('.csv'):
            return self._read_csv_file(file_path)
        elif filename.endswith('.txt'):
            return self._read_tsv_file(file_path)  # Assuming .txt files are tab-separated
        else:
            raise ValueError(f"Unsupported file type for file: {filename}")

    def _read_soft_file(self, file_path):
        """Read SOFT format file"""
        with open(file_path, 'r') as f:
            content = f.read()
        metadata = {}
        for line in content.split('\n'):
            if line.startswith('!'):
                key, value = line[1:].split('=', 1)
                metadata[key.strip()] = value.strip()
        return {'metadata': metadata}

    def _read_tar_file(self, file_path):
        """Read tar archive"""
        with tarfile.open(file_path, 'r:*') as tar:
            files = {}
            for member in tar.getmembers():
                if member.isfile():
                    files[member.name] = tar.extractfile(member).read().decode('utf-8')
            return files

    def _read_tsv_file(self, file_path):
        """Read TSV file"""
        return pd.read_csv(file_path, sep='\t', on_bad_lines='skip', engine='python')

    def _read_csv_file(self, file_path):
        """Read CSV file"""
        return pd.read_csv(file_path)

# Primer - izpis podatkov:
if __name__ == "__main__":
    # Set the base directory where the files are located
    base_dir = r"data_file"
    geo_access = GEODataAccess(base_dir=base_dir)

    # Load specific files
    family_soft = geo_access.load_file("GSE171485_family.soft")
    print("Family SOFT metadata:")
    print(family_soft['metadata'])

    family_xml = geo_access.load_file("GSE171485_family.xml.tar")
    print("\nFamily XML contents:")
    print(family_xml.keys())

    fpkm_data = geo_access.load_file("GSE171485_norm_counts_FPKM_GRCh38.p13_NCBI.tsv")
    print("\nFPKM data preview:")
    print(fpkm_data.head())

    tpm_data = geo_access.load_file("GSE171485_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
    print("\nTPM data preview:")
    print(tpm_data.head())

    pdac_tissue_data = geo_access.load_file("GSE171485_PDAC-tissue-ExpressionMatrix.csv")
    print("\nPDAC tissue data preview:")
    print(pdac_tissue_data.head())

    raw_counts_data = geo_access.load_file("GSE171485_raw_counts_GRCh38.p13_NCBI.tsv")
    print("\nRaw counts data preview:")
    print(raw_counts_data.head())

    series_matrix = geo_access.load_file("GSE171485_series_matrix.txt")
    print("\nSeries matrix preview:")
    print(series_matrix.head())


# %%
# Comparative table display function (bolj pregledno za biologe)
if __name__ == "__main__":
    # Set the base directory where the files are located
    base_dir = r"DataFile/"
    geo_access = GEODataAccess(base_dir=base_dir)

    # Display data in a clear table format
    def display_data_tables():
        print("Family SOFT Metadata:")
        metadata_df = pd.DataFrame(list(family_soft['metadata'].items()), columns=['Key', 'Value'])
        display(metadata_df)

        print("\nFPKM Data Preview:")
        display(fpkm_data.head())

        print("\nTPM Data Preview:")
        display(tpm_data.head())

        print("\nPDAC Tissue Data Preview:")
        display(pdac_tissue_data.head())

        print("\nRaw Counts Data Preview:")
        display(raw_counts_data.head())

        print("\nSeries Matrix Preview:")
        display(series_matrix.head())

    # Call the function to display tables
    display_data_tables()



