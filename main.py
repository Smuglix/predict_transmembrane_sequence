# memebrained: Ilya Adamskiy and Tal Or

import math
import random
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
    problematic_indices = []
    transmembrane_sequences = []
    transmembrane_sequences_superlist = []
    for index, row in df.iterrows():
        try:
            sequence = row['Sequence']
            matches = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
            transmembrane_indexes = [(int(start), int(end)) for start, end in matches]
            for index in transmembrane_indexes:
                transmembrane_sequence = (sequence[index[0] - 1:index[1]])
                transmembrane_sequences.append(transmembrane_sequence)
            transmembrane_sequences_superlist.append(transmembrane_sequences)
            transmembrane_sequences = []
        except Exception as e:
            problematic_indices.append(index)
            print(f"Error processing row {index}: {e}")
        # Drop problematic rows
    df = df.drop(problematic_indices)
    df = df.assign(Transmembrane_sequences=transmembrane_sequences_superlist)
    return df


def separate_testing_data(df):
    testing_dataframe = pd.DataFrame(columns=df.columns)
    testing_indexes = []
    data_length = len(df)
    testing_data_length = int(data_length * 0.1)
    for _ in range(testing_data_length):
        print(_)
        random_num = random.randrange(0, data_length)
        testing_dataframe = pd.concat([testing_dataframe, df.iloc[[random_num]]])
        df = df.drop(df.index[random_num])
        data_length -= 1
    # print(training_dataframe)
    return testing_dataframe, df


def calc_p(training_df, amino_acid):
    total_training_sequence = ''
    for index, row in training_df.iterrows():
        total_training_sequence += (row['Sequence'])

    total_occurrence_of_specified_amino_acid = re.findall(amino_acid, total_training_sequence)

    p = len(total_occurrence_of_specified_amino_acid) / len(total_training_sequence)
    return p


def calc_q(training_df, amino_acid):
    total_transmembrane_amino_acid_sequence = ''

    for index, row in training_df.iterrows():
        for transmembrane_sequence in row['Transmembrane_sequences']:
            total_transmembrane_amino_acid_sequence += transmembrane_sequence
    total_occurrence_of_specified_transmembrane_amino_acid_occurrence = re.findall(amino_acid, total_transmembrane_amino_acid_sequence)

    q = len(total_occurrence_of_specified_transmembrane_amino_acid_occurrence) / len(total_transmembrane_amino_acid_sequence)
    return q


def calc_s(q, p):
    s = math.log(q / p)
    return s


def create_dictionary_with_s_values_for_all_amino_acids(amino_acids_list, df):
    """
    Jones DT, Taylor WR, Thornton JM.
    A model recognition approach to the prediction of all-helical membrane protein structure and topology.
    Biochemistry. 1994 Mar 15;33(10):3038-49. doi: 10.1021/bi00176a037. PMID: 8130217.
    page - 2
    """
    amino_acid_s_value_dict = {}
    for amino_acid in amino_acids_list:
        pi = calc_p(df, amino_acid)
        qi = calc_q(df, amino_acid)
        si = calc_s(qi, pi)
        amino_acid_s_value_dict[amino_acid] = si

    return amino_acid_s_value_dict


type = "beta"
amino_acids_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

download_transmembrane_protein_data(type)
transmembrane_protein_df = read_data_into_dataframe(f"transmembrane_{type}.xlsx")
training_dataframe, testing_dataframe = separate_testing_data(extract_transmembrane_sequences(transmembrane_protein_df))
calc_p(training_dataframe, 'A')
calc_q(training_dataframe, 'A')
