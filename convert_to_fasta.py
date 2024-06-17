import pandas as pd

def csv_to_fasta(input_csv, output_fasta):
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Open the output FASTA file
    with open(output_fasta, 'w') as fasta_file:
        # Iterate over each row in the DataFrame
        for index, row in df.iterrows():
            # Write the entry and sequence in FASTA format
            fasta_file.write(f">{row['Entry']}\n")
            fasta_file.write(f"{row['Sequence']}\n")


input_csv = 'testing_dataframe.csv'
output_fasta = 'output.fasta'
csv_to_fasta(input_csv, output_fasta)