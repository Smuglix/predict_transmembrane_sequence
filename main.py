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
            for t_index in transmembrane_indexes:
                transmembrane_sequence = (sequence[t_index[0] - 1:t_index[1] - 1])
                print(len(transmembrane_sequence))
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
    total_occurrence_of_specified_transmembrane_amino_acid_occurrence = re.findall(amino_acid,
                                                                                   total_transmembrane_amino_acid_sequence)

    q = len(total_occurrence_of_specified_transmembrane_amino_acid_occurrence) / len(
        total_transmembrane_amino_acid_sequence)
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


def evaluate_answer(predicted_range, ground_truth_range):
    """
    Wikipedia contributors. (2024, April 26).
    Jaccard index. In Wikipedia, The Free Encyclopedia.
    Retrieved 09:49, May 29, 2024,
    from https://en.wikipedia.org/w/index.php?title=Jaccard_index&oldid=1220812875
    """

    # Calculate intersection
    intersection_start = max(ground_truth_range[0], predicted_range[0])
    intersection_end = min(ground_truth_range[1], predicted_range[1])
    intersection_length = max(0, intersection_end - intersection_start)

    # Calculate union
    union_start = min(ground_truth_range[0], predicted_range[0])
    union_end = max(ground_truth_range[1], predicted_range[1])
    union_length = union_end - union_start

    # Calculate IoU
    if union_length == 0:  # To avoid division by zero
        return 0.0
    iou = intersection_length / union_length

    # Convert IoU to percentage
    percentage_similarity = iou * 100
    return percentage_similarity


"""def predict_transmembrane_range(sequence, s_values, window_size):
    sequence = list(sequence)
    print(sequence)
    s_val_to_seq = (pd.Series(sequence)).map(s_values)
    list_of_seq_to_s_val = list(s_val_to_seq)
    print(list_of_seq_to_s_val)
    averaged_sequence = list_of_seq_to_s_val
    if len(sequence) <= window_size:
        print("errm, what the sigma? ಠಿ_ಠ")
        return -1
    for i in range(len(sequence) - window_size):
        middle_amino_acid = int(window_size / 2) + i
        start_amino_acid = middle_amino_acid - int(window_size / 2)
        end_amino_acid = middle_amino_acid + int(window_size / 2)

        sum_of_s = sum(list_of_seq_to_s_val[start_amino_acid:end_amino_acid+1])
        averaged_sequence[middle_amino_acid] = sum_of_s/window_size
    print(averaged_sequence)"""


"""def predict_transmembrane_range(sequence, s_values, start_window_size):
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
            amino_acid_range = (int(start_amino_acid), int(end_amino_acid))
            big_sums[amino_acid_range] = sum_of_s
        window_size -= 1
    big_sums_sorted = dict(sorted(big_sums.items(), key=lambda item: item[1], reverse=True))
    print(big_sums_sorted)"""

"""def predict_and_check_best_result(df, start_window_size, s_values):
    for index, row in df.iterrows():
        sequence = row['Sequence']
        matches = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
        transmembrane_indexes = [(int(start), int(end)) for start, end in matches]
        prediction = predict_transmembrane_range(sequence, s_values, start_window_size)"""


def filter_overlapping_ranges(big_sums_sorted):
    selected_ranges = []
    for current_range, current_value in big_sums_sorted.items():
        overlap = False
        for selected_range in selected_ranges:
            if ranges_overlap(current_range, selected_range):
                overlap = True
                break
        if not overlap:
            selected_ranges.append(current_range)
    return selected_ranges


def ranges_overlap(range1, range2):
    return not (range1[1] < range2[0] or range1[0] > range2[1])


def predict_transmembrane_range(sequence, s_values, start_window_size):
    window_size = start_window_size
    big_sums = {}
    sequence = list(sequence)
    s_val_to_seq = (pd.Series(sequence)).map(s_values)
    list_of_seq_to_s_val = list(s_val_to_seq)
    if len(sequence) <= window_size:
        print("errm, what the sigma? ಠಿ_ಠ")
        return -1
    while window_size > 17:
        for i in range(len(sequence) - window_size):
            middle_amino_acid = int(window_size / 2) + i
            start_amino_acid = middle_amino_acid - int(window_size / 2)
            end_amino_acid = middle_amino_acid + int(window_size / 2)

            sum_of_s = sum(list_of_seq_to_s_val[start_amino_acid:end_amino_acid + 1])
            amino_acid_range = (int(start_amino_acid), int(end_amino_acid))
            big_sums[amino_acid_range] = sum_of_s
        window_size -= 1
    big_sums_sorted = dict(sorted(big_sums.items(), key=lambda item: item[1], reverse=True))

    filtered_ranges = filter_overlapping_ranges(big_sums_sorted)
    for item in filtered_ranges:
        print(f"{item}: {big_sums_sorted[item]}")
    return filtered_ranges


# Example usage
type = "helical"
window_size = 40
amino_acids_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
download_transmembrane_protein_data(type)
transmembrane_protein_df = read_data_into_dataframe(f"transmembrane_{type}.xlsx")

if not os.path.exists('training_dataframe.csv') and not os.path.exists('testing_dataframe.csv'):
    print("Didn't find saved dataframes, so im gonna create em ")
    training_dataframe, testing_dataframe = separate_testing_data(
        extract_transmembrane_sequences(transmembrane_protein_df))

training_dataframe = pd.read_csv('training_dataframe.csv')
testing_dataframe = pd.read_csv('testing_dataframe.csv')

s_values = create_dictionary_with_s_values_for_all_amino_acids(amino_acids_list, training_dataframe)
print(s_values)
print(training_dataframe['Entry'].iloc[30560], training_dataframe['Sequence'].iloc[30560])
predict_transmembrane_range(training_dataframe['Sequence'].iloc[30560], s_values, window_size)
