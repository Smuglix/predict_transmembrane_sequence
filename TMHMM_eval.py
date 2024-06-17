import csv
import random
import requests
import pandas as pd
import os
import re
from bs4 import BeautifulSoup

def download_transmembrane_protein_data(type):
    url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Csequence%2Cft_transmem&format=xlsx&query=((ft_transmem%3A{type}))+AND+(reviewed%3Atrue)"
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
            return None

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
            for t_index in transmembrane_indexes:
                transmembrane_sequence = (sequence[t_index[0] - 1:t_index[1] - 1])
                print(len(transmembrane_sequence))
                transmembrane_sequences.append(transmembrane_sequence)
            transmembrane_sequences_superlist.append(transmembrane_sequences)
            transmembrane_sequences = []
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
        print(_)
        random_num = random.randrange(0, data_length)
        testing_dataframe = pd.concat([testing_dataframe, df.iloc[[random_num]]])
        df = df.drop(df.index[random_num])
        training_dataframe = df
        data_length -= 1
    testing_dataframe.to_csv('testing_dataframe.csv')
    training_dataframe.to_csv('training_dataframe.csv')
    return testing_dataframe, training_dataframe

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

def get_transmembrane_range_from_row(row):
    matches = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
    transmembrane_indexes = [(int(start), int(end)) for start, end in matches]
    return transmembrane_indexes

def extract_tmhmm_predictions(html_file_path):
    with open(html_file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()
    soup = BeautifulSoup(html_content, 'html.parser')
    pre_tag = soup.find('pre')
    protein_info = pre_tag.get_text().strip().split('\n')
    tmhmm_predictions = {}
    for line in protein_info:
        parts = line.split('\t')
        protein_id = parts[0]
        topology = parts[5].split('=')[1]
        matches = re.findall(r'(\d+)-(\d+)', topology)
        transmembrane_ranges = [(int(start), int(end)) for start, end in matches]
        tmhmm_predictions[protein_id] = transmembrane_ranges
    return tmhmm_predictions

def mario_steals_all_your_bitches_with_tmhmm(df, tmhmm_predictions):
    sum_perc_num = 0
    total_transmembrane_evaluated = 0
    missed_transmembrane_regions_count = 0
    for index, row in df.iterrows():
        percentage_similarities = []
        print(index)
        entry_id = row['Entry']
        ground_truth_ranges = get_transmembrane_range_from_row(row)
        if entry_id in tmhmm_predictions:
            tmhmm_ranges = tmhmm_predictions[entry_id]
            for tmhmm_range in tmhmm_ranges:
                largest_perc = 0
                total_transmembrane_evaluated += 1
                for transmembrane_range in ground_truth_ranges:
                    similarity_percentage = evaluate_answer(tmhmm_range, transmembrane_range)
                    percentage_similarities.append(similarity_percentage)
                largest_perc = max(percentage_similarities) if percentage_similarities else 0
                percentage_similarities = []
                print(f"TMHMM Prediction: {tmhmm_range}, Similarity: {largest_perc}")
                sum_perc_num += largest_perc
            missed_transmembrane_regions_count += len(ground_truth_ranges) - len(tmhmm_ranges)
        else:
            print(f"No TMHMM prediction for {entry_id}")
            missed_transmembrane_regions_count += len(ground_truth_ranges)
    average_with_missed_regions = sum_perc_num / (missed_transmembrane_regions_count + total_transmembrane_evaluated)
    print(f"sum_perc_num: {sum_perc_num}")
    print(f"Evaluated {total_transmembrane_evaluated} transmembrane regions, missed {missed_transmembrane_regions_count} regions")
    average_perc = sum_perc_num / (total_transmembrane_evaluated + 1)
    print(f"YOUR AVERAGE best accuracy IS: {average_perc} %")
    print(f"YOUR AVERAGE chance that predicted transmembrane region is on spot with what should be is: {average_with_missed_regions}%")

# Main script
type = "helical"
download_transmembrane_protein_data(type)
transmembrane_protein_df = read_data_into_dataframe(f"transmembrane_{type}.xlsx")

if not os.path.exists('training_dataframe.csv') and not os.path.exists('testing_dataframe.csv'):
    print("Didn't find saved dataframes, so im gonna create em ")
    training_dataframe, testing_dataframe = separate_testing_data(
        extract_transmembrane_sequences(transmembrane_protein_df))

training_dataframe = pd.read_csv('training_dataframe.csv')
testing_dataframe = pd.read_csv('testing_dataframe.csv')

# Extract TMHMM predictions from HTML file
tmhmm_predictions = extract_tmhmm_predictions('TMHMM result_one_line.html')

# Evaluate against TMHMM predictions
mario_steals_all_your_bitches_with_tmhmm(testing_dataframe, tmhmm_predictions)