import os
import pandas as pd

# Set the directory containing the subdirectories
parent_dir = '/home/ocanal/Escritorio/Trombosi'

# Initialize an empty list to store the data
all_data = []

# Iterate over the subdirectories
for subdir in os.listdir(parent_dir):
    subdir_path = os.path.join(parent_dir, subdir)

    # Check if the current item is a directory
    if os.path.isdir(subdir_path):
        # Check if the subdirectory contains an "Excel" directory
        excel_dir_path = os.path.join(subdir_path, 'excels')
        if os.path.isdir(excel_dir_path):
            # Iterate over the Excel files in the excels directory
            for file in os.listdir(excel_dir_path):
                # Check if the file is an Excel file
                if file.endswith('.xlsx'):
                    file_path = os.path.join(excel_dir_path, file)
                    # Read the Excel file into a pandas DataFrame
                    df = pd.read_excel(file_path)
                    # Append the data to the list
                    all_data.append(df)

#create a joined_excels folder, if it haven't already exists

output_directory = f"{parent_dir}/joined_excels"
dir_exists = os.path.exists(output_directory)
if dir_exists:
	pass 

else:
	os.mkdir(output_directory)

# Concatenate all the data into a single DataFrame
all_data_df = pd.concat(all_data, ignore_index=True)
all_data_sorted = all_data_df.sort_values('Uploaded_variation')
# Save the DataFrame to an Excel file
all_data_df.to_excel(f'{output_directory}/variants_joined.xlsx', index=False)
all_data_sorted.to_excel(f'{output_directory}/variants_joined_sorted.xlsx', index = False)
