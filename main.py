# memebrained: Ilya Adamskiy and Tal Or

import csv
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
        training_dataframe = df
        data_length -= 1
    # print(training_dataframe)
    testing_dataframe.to_csv('testing_dataframe.csv')
    training_dataframe.to_csv('training_dataframe.csv')
    return testing_dataframe, training_dataframe


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

    # Check if 's_values.csv' already exists and load it if so
    if os.path.exists('s_values.csv'):
        print("s values already calculated, loading from s_values.csv")
        with open('s_values.csv', newline='') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)  # Read the header row
            for row in reader:
                amino_acid = row[0]
                s_value = float(row[1])
                amino_acid_s_value_dict[amino_acid] = s_value
        return amino_acid_s_value_dict

    # Calculate s values for each amino acid
    for amino_acid in amino_acids_list:
        pi = calc_p(df, amino_acid)
        qi = calc_q(df, amino_acid)
        si = calc_s(qi, pi)
        amino_acid_s_value_dict[amino_acid] = si

    # Save the dictionary to 's_values.csv'
    with open('s_values.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Amino Acid', 'S Value'])  # Write the header
        for key, value in amino_acid_s_value_dict.items():
            writer.writerow([key, value])  # Write each key-value pair as a row

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
    for current_range in big_sums_sorted.keys():
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


print("Ilgyatt adamskibidi")


def average_sequence(list_of_seq_to_s_val, averaging_window_size):
    averaged_sequence = list_of_seq_to_s_val
    if len(list_of_seq_to_s_val) <= averaging_window_size:
        print("errm, what the sigma? ಠಿ_ಠ")
        return -1
    for i in range(len(list_of_seq_to_s_val) - averaging_window_size):
        middle_amino_acid = int(averaging_window_size / 2) + i
        start_amino_acid = middle_amino_acid - int(averaging_window_size / 2)
        end_amino_acid = middle_amino_acid + int(averaging_window_size / 2)

        sum_of_s = sum(list_of_seq_to_s_val[start_amino_acid:end_amino_acid + 1])
        averaged_sequence[middle_amino_acid] = sum_of_s / averaging_window_size
    return averaged_sequence


def predict_transmembrane_range(sequence, s_values, start_window_size, end_window_size, averaging_window_size, row):
    window_size = start_window_size
    big_sums = {}
    filtered_ranges = {}
    sequence = list(sequence)
    s_val_to_seq = (pd.Series(sequence)).map(s_values)
    list_of_seq_to_s_val = list(s_val_to_seq)

    averaged_sequence = average_sequence(list_of_seq_to_s_val, averaging_window_size)

    if len(sequence) <= window_size:
        print("errm, what the sigma? ಠಿ_ಠ")
        window_size = len(sequence)

    while window_size > end_window_size:
        for i in range(len(sequence) - window_size):
            middle_amino_acid = int(window_size / 2) + i
            start_amino_acid = middle_amino_acid - int(window_size / 2)
            end_amino_acid = middle_amino_acid + int(window_size / 2)

            sum_of_s = sum(averaged_sequence[start_amino_acid:end_amino_acid + 1])
            amino_acid_range = (int(start_amino_acid), int(end_amino_acid))
            big_sums[amino_acid_range] = sum_of_s

        window_size -= 1

    big_sums_sorted = dict(sorted(big_sums.items(), key=lambda item: item[1], reverse=True))

    selected_ranges = filter_overlapping_ranges(big_sums_sorted)
    for selected_range in selected_ranges:
        filtered_ranges[selected_range] = big_sums_sorted[selected_range]
    filtered_ranges_dataframe = pd.DataFrame.from_dict(filtered_ranges, orient='index')
    os.makedirs('filtered_ranges', exist_ok=True)
    filtered_ranges_dataframe.to_csv(f'filtered_ranges/{row["Entry"]}.csv')
    return filtered_ranges


def get_transmembrane_range_from_row(row):
    matches = re.findall(r"TRANSMEM\s(\d+)..(\d+)", row['Transmembrane'])
    transmembrane_indexes = [(int(start), int(end)) for start, end in matches]
    return transmembrane_indexes


def mario_steals_your_best_liver(df, s_values, start_window_size, end_window_size, averaging_window_size):
    sum_perc = 0
    for index, row in df.iterrows():
        percentage_similarities = []
        print(index)
        sequence = row['Sequence']
        prediction = predict_transmembrane_range(sequence, s_values, start_window_size, end_window_size,
                                                 averaging_window_size, row)
        ground_truth_ranges = get_transmembrane_range_from_row(row)

        try:
            max_prediction = max(prediction.items(), key=lambda x: x[1])
            for transmembrane_range in ground_truth_ranges:
                similarity_percentage = evaluate_answer(max_prediction[0], transmembrane_range)
                percentage_similarities.append(similarity_percentage)
            largest_perc = max(percentage_similarities)
            print(f"{max_prediction[1]}: {largest_perc}")
            sum_perc += largest_perc
        except ValueError as e:
            print(f"Skipping index {index} due to error: {e}")

    average_perc = sum_perc / (index + 1)
    print(f"YOUR AVERAGE cock size IS: {average_perc} milimeters")


def mario_steals_all_your_bitches(df, s_values, start_window_size, end_window_size, averaging_window_size, prediction_score_minimum):
    sum_perc = 0
    transmembrane_evaluated = 0
    for index, row in df.iterrows():
        percentage_similarities = []
        print(index)
        sequence = row['Sequence']
        predictions = predict_transmembrane_range(sequence, s_values, start_window_size, end_window_size,
                                                 averaging_window_size, row)
        ground_truth_ranges = get_transmembrane_range_from_row(row)

        try:
            for transmembrane_prediction_range, transmembrane_prediction_score in predictions.items():
                largest_perc = 0
                if transmembrane_prediction_score < prediction_score_minimum:
                    continue
                transmembrane_evaluated += 1
                for transmembrane_range in ground_truth_ranges:
                    similarity_percentage = evaluate_answer(transmembrane_prediction_range, transmembrane_range)
                    percentage_similarities.append(similarity_percentage)
                largest_perc = max(percentage_similarities)
                percentage_similarities = []
                print(f"{transmembrane_prediction_score}: {largest_perc}")
                sum_perc += largest_perc
        except ValueError as e:
            print(f"Skipping index {index} due to error: {e}")

    average_perc = sum_perc / (transmembrane_evaluated + 1)
    print(f"YOUR AVERAGE cock size IS: {average_perc} milimeters")


type = "helical"
window_size = 40
end_window_size = 19
averaging_window_size = 7
prediction_score_minimum = 0

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
mario_steals_all_your_bitches(testing_dataframe, s_values, window_size, end_window_size, averaging_window_size, prediction_score_minimum)
