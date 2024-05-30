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
            df = pd.read_excel(file_path)
        else:
            print("Error downloading file:", response.status_code)
            return None
    return df

def read_data_into_dataframe(file_path):
    df = pd.read_excel(file_path)
    return df

def extract_transmembrane_sequences(df):
    problematic_indices = []
    transmembrane_sequences_superlist = []
    for index, row in df.iterrows():
        try:
            sequence = row['Sequence']
            matches = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
            transmembrane_indexes = [(int(start), int(end)) for start, end in matches]
            transmembrane_sequences = [sequence[start-1:end] for start, end in transmembrane_indexes]
            transmembrane_sequences_superlist.append(transmembrane_sequences)
        except Exception as e:
            problematic_indices.append(index)
            print(f"Error processing row {index}: {e}")
    df = df.drop(problematic_indices)
    df = df.assign(Transmembrane_sequences=transmembrane_sequences_superlist)
    return df

def separate_testing_data(df):
    testing_dataframe = pd.DataFrame(columns=df.columns)
    data_length = len(df)
    testing_data_length = int(data_length * 0.1)
    for _ in range(testing_data_length):
        random_num = random.randrange(0, data_length)
        testing_dataframe = pd.concat([testing_dataframe, df.iloc[[random_num]]])
        df = df.drop(df.index[random_num])
        data_length -= 1
    training_dataframe = df
    testing_dataframe.to_csv('testing_dataframe.csv', index=False)
    training_dataframe.to_csv('training_dataframe.csv', index=False)
    return testing_dataframe, training_dataframe

def calc_p(training_df, amino_acid):
    total_training_sequence = ''.join(training_df['Sequence'])
    total_occurrence_of_specified_amino_acid = total_training_sequence.count(amino_acid)
    p = total_occurrence_of_specified_amino_acid / len(total_training_sequence)
    return p

def calc_q(training_df, amino_acid):
    total_transmembrane_amino_acid_sequence = ''.join([''.join(seqs) for seqs in training_df['Transmembrane_sequences']])
    total_occurrence_of_specified_transmembrane_amino_acid_occurrence = total_transmembrane_amino_acid_sequence.count(amino_acid)
    q = total_occurrence_of_specified_transmembrane_amino_acid_occurrence / len(total_transmembrane_amino_acid_sequence)
    return q

def calc_s(q, p):
    s = math.log(q / p)
    return s

def create_dictionary_with_s_values_for_all_amino_acids(amino_acids_list, training_df):
    amino_acid_s_value_dict = {}
    for amino_acid in amino_acids_list:
        pi = calc_p(training_df, amino_acid)
        qi = calc_q(training_df, amino_acid)
        si = calc_s(qi, pi)
        amino_acid_s_value_dict[amino_acid] = si
    return amino_acid_s_value_dict

def predict_transmembrane_range(sequence, s_values, start_window_size):
    window_size = start_window_size
    big_sums = {}
    sequence = list(sequence)
    s_val_to_seq = (pd.Series(sequence)).map(s_values)
    list_of_seq_to_s_val = list(s_val_to_seq)
    if len(sequence) <= window_size:
        print("errm, what the sigma? ಠಿ_ಠ")
        return -1
    while window_size > 2:
        for i in range(len(sequence) - window_size):
            middle_amino_acid = int(window_size / 2) + i
            start_amino_acid = middle_amino_acid - int(window_size / 2)
            end_amino_acid = middle_amino_acid + int(window_size / 2)
            sum_of_s = sum(list_of_seq_to_s_val[start_amino_acid:end_amino_acid + 1])
            amino_acid_range = (start_amino_acid, end_amino_acid)
            big_sums[amino_acid_range] = sum_of_s
        window_size -= 1
    big_sums_sorted = dict(sorted(big_sums.items(), key=lambda item: item[1], reverse=True))
    return list(big_sums_sorted.keys())[0]

def evaluate_answer(predicted_range, ground_truth_range):
    intersection_start = max(ground_truth_range[0], predicted_range[0])
    intersection_end = min(ground_truth_range[1], predicted_range[1])
    intersection_length = max(0, intersection_end - intersection_start)
    union_start = min(ground_truth_range[0], predicted_range[0])
    union_end = max(ground_truth_range[1], predicted_range[1])
    union_length = union_end - union_start
    if union_length == 0:
        return 0.0
    iou = intersection_length / union_length
    percentage_similarity = iou * 100
    return percentage_similarity

type = "helical"
window_size = 19
amino_acids_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
df = download_transmembrane_protein_data(type)
if df is None:
    raise ValueError("Failed to download or process the data.")

transmembrane_protein_df = read_data_into_dataframe(f"transmembrane_{type}.xlsx")

if not os.path.exists('training_dataframe.csv') or not os.path.exists('testing_dataframe.csv'):
    print("Didn't find saved dataframes, so I'm going to create them.")
    training_dataframe, testing_dataframe = separate_testing_data(extract_transmembrane_sequences(transmembrane_protein_df))
else:
    training_dataframe = pd.read_csv('training_dataframe.csv')
    testing_dataframe = pd.read_csv('testing_dataframe.csv')

training_dataframe['Transmembrane_sequences'] = training_dataframe['Transmembrane_sequences'].apply(eval)
testing_dataframe['Transmembrane_sequences'] = testing_dataframe['Transmembrane_sequences'].apply(eval)

s_values = create_dictionary_with_s_values_for_all_amino_acids(amino_acids_list, training_dataframe)

# Evaluation
total_similarity = 0
similarity_count = 0
zero_similarity_count = 0

for index, row in testing_dataframe.iterrows():
    sequence = row['Sequence']
    ground_truth_ranges = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
    ground_truth_ranges = [(int(start)-1, int(end)-1) for start, end in ground_truth_ranges]
    predicted_range = predict_transmembrane_range(sequence, s_values, window_size)
    if predicted_range == -1:
        continue
    for ground_truth_range in ground_truth_ranges:
        similarity = evaluate_answer(predicted_range, ground_truth_range)
        if similarity > 0:
            total_similarity += similarity
            similarity_count += 1
        else:
            zero_similarity_count += 1
        print(f"Entry: {row['Entry']}, Predicted: {predicted_range}, Ground Truth: {ground_truth_range}, Similarity: {similarity:.2f}%")

average_similarity = total_similarity / similarity_count if similarity_count > 0 else 0
print(f"\nAverage Similarity: {average_similarity:.2f}%")
print(f"Total predictions with 0% similarity: {zero_similarity_count}")
