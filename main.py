# memebrained: Ilya Adamskiy and Tal Or

import requests
import pandas as pd
import os
import re


def download_transmembrane_protein_data(type):
    url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Csequence%2Cft_transmem&format=xlsx&query=%28%28ft_transmem%3A{type}%29%29+AND+%28reviewed%3Atrue%29"
    file_path = f"transmembrane_{type}.xlsx"
    if os.path.exists(file_path):
        print("File already exists, skipping download.")
        df = pd.read_excel(file_path)
    else:
        response = requests.get(url, stream=True)

        if response.status_code == 200:
            with open(file_path, "wb") as f:
                for chunk in response.iter_content(1024):
                    f.write(chunk)
            print("File downloaded successfully!")
        else:
            print("Error downloading file:", response.status_code)
            return None  # or raise an exception, depending on your use case


def read_data_into_dataframe(file_path):
    df = pd.read_excel(file_path)
    return df


def extract_transmembrane_sequences(df):

    i = 0
    for index, row in df.iterrows():
        if i > 0:
            break
        sequence = row['Sequence']
        matches = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
        transmembrane_indexes = [(int(start), int(end)) for start, end in matches]
        print(sequence[transmembrane_indexes[0][0]-1:transmembrane_indexes[0][1]])
        i+=1


type = "helical"

download_transmembrane_protein_data(type)
transmembrane_protein_df = read_data_into_dataframe(f"transmembrane_{type}.xlsx")
extract_transmembrane_sequences(transmembrane_protein_df)
