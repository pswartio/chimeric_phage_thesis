import re
import os
import pandas as pd 
import csv
import numpy as np
import math

input_dir_path = "/Users/pimswart/alignment_hittables/"
hit_genome_dir = '/Users/pimswart/alignmentsfrom_hittables_collected/'
output_dir_path = '/Users/pimswart/hittables_prophageregion_collected_new/'

for file_name in os.listdir(input_dir_path):
    if file_name.endswith(".csv"):
        path = os.path.join(input_dir_path, file_name)
        filee = os.path.basename(path)
        with open(path, 'r') as f:

            # Create a CSV reader object
            reader = csv.reader(f)

            # Initialize a dictionary to store the collections
            collections = {}

            # Loop over each line in the file
            for row in reader:
                # Skip empty rows
                if not row:
                    continue

                # Extract the second column value
                key = row[1] if len(row) > 1 else None

                # Extract the last column value
                value = float(row[-1]) if row[-1] else None

                # Check if the collection exists in the dictionary
                if key in collections:
                    # Add the row to the collection
                    collections[key].append(row)
                else:
                    # Create a new collection with the first row
                    collections[key] = [row]

            # Initialize a list to store the filtered rows
            filtered_rows = []

            # Loop over each collection
            for key in collections:
                # Sort the collection based on the last column value
                sorted_collection = sorted(collections[key], key=lambda x: float(x[-1]) if x[-1] else None, reverse=True)

                # Add the two rows with the highest values to the filtered rows list
                filtered_rows.extend(sorted_collection[:2])

            # Now you have a list of the filtered rows
            # You can then extract the data you need from these rows
            for row in filtered_rows:
                column_2 = row[1] if len(row) > 1 else None # second column
                column_6 = row[5] if len(row) > 5 else None # sixth column
                column_7 = row[6] if len(row) > 6 else None # seventh column

        # print(filtered_rows)

        # loop over each sublist
        for sublist in filtered_rows:
            # remove first entry of sublist
            del sublist[0]

        # print(filtered_rows)
        first_entries = []

        # loop over each sublist
        for sublist in filtered_rows:
            # check if first entry is already in list
            if sublist[0] not in first_entries:
                # add first entry to list
                first_entries.append(sublist[0])

        # print(first_entries)

        new_dict = {}

        # loop over each first entry to create dictionary keys
        for entry in first_entries:
            # create empty list for values
            values_list = []
            # loop over each sublist to find matching strings
            for sublist in filtered_rows:
                if sublist[0] == entry:
                    values_list.append(sublist)
            # add key-value pair to dictionary
            new_dict[entry] = values_list

        # print(new_dict)

        df1 = pd.DataFrame(columns=['hit_genome', 'coord1', 'coord2'])

        for key, value in new_dict.items():
            
            # Loop through the two lists in the value
            for sublist in value:
                
                # Add a new row to the dataframe
                new_row = {'hit_genome': sublist[0], 'coord1': sublist[7], 'coord2': sublist[8], 'overlap?': np.nan}
                df1 = df1.append(new_row, ignore_index=True)

        # add 10 new columns to the dataframe
        for i in range(1):
          	 col_name = f'prophage_regions_{i}'
          	 df1[col_name] = np.nan


        # Show the resulting dataframe
        # print('testdf1',df1)
        df1['prophage_regions_0'] = df1['prophage_regions_0'].astype('object')

        appended_data = []
        for filename in os.listdir(hit_genome_dir+filee[:-4]):
        	if ".txt" in filename:
        		# print(filename)
        		# Open the file for reading
        		f = os.path.join(hit_genome_dir+filee[:-4], filename)
        		filename = filename.removesuffix('_summary.txt')
        		# print(filename)

        		# Open the file for reading
        		with open(f, 'r') as filea:

        		    # Read the contents of the file into a list
        		    lines = filea.readlines()

        		    # Initialize empty lists for the extracted data
        		    filenames = []
        		    strings = []
        		    numbers = []

        		    # Loop through each line in the list
        		    for line in lines:

        		        # Use regular expressions to search for the two desired string patterns
        		        pattern1 = re.search(r'\d+-\d+', line)
        		        pattern2 = re.search(r'(intact|questionable|incomplete)\(\d+\)', line)

        		        # If both patterns are found, extract the matched strings and append them to the lists
        		        if pattern1 and pattern2:
        		            filenames.append(filename)
        		            strings.append(pattern2.group(1))
        		            numbers.append(pattern1.group())

        		    # Create a dictionary with the extracted data
        		    data = {'ref_hit_genome': filenames, 'region_status': strings, 'Numbers': numbers}

        		    # Create a pandas DataFrame with the data
        		    df = pd.DataFrame(data)

        			# split the 'Numbers' column into two new columns 'Number_1' and 'Number_2'
        		    df[['coord3', 'coord4']] = df['Numbers'].str.split('-', expand=True)

        			# drop the original 'Numbers' column
        		    df = df.drop('Numbers', axis=1)
        		    appended_data.append(df)
        			# print the resulting dataframe
        		    # print(df.head())
        		    # print('dfffff',df)
        df3 = pd.concat(appended_data)
        print('final',df3.to_string())


        for i, row1 in df1.iterrows():
        	for j, row2 in df3.iterrows():
        		roww = row2.tolist()
        		# print('roww',roww)
        		# print('tesss',row1['hit_genome'],row2['ref_hit_genome'])
        		if row1['hit_genome'] == row2['ref_hit_genome']:
        			if isinstance(df1.at[i,'prophage_regions_0'],list):
        				# print(True)
        				# print('df1hitgenom',df1['hit_genome'])
        				df1.at[i,'prophage_regions_0'].extend(roww)
        			else:
        				# print(False)
        				df1.at[i,'prophage_regions_0'] = roww






        			# print('lala',row1['hit_genome'],row2['ref_hit_genome'])
        			# df1.at[i,'prophage_regions_0'] = roww
        		else:
        			continue

        # Change the NaN's in column to empty lists
        isna = df1['prophage_regions_0'].isna()
        df1.loc[isna, 'prophage_regions_0'] = pd.Series([[]] * isna.sum()).values


        for index, row in df1.iterrows():
            # check if the list in column `prophage_regions_0` is empty
            if not row['prophage_regions_0']:
                continue
            # create a new list of lists containing every four items in the original list
            new_list = [row['prophage_regions_0'][i:i+4] for i in range(0, len(row['prophage_regions_0']), 4)]
            for sublist in new_list:
            	sublist[2]=int(sublist[2])
            	sublist[3]=int(sublist[3])
            # update the value in the dataframe
            df1.at[index, 'prophage_regions_0'] = new_list

        print('lastdf1')
        print(df1.to_string())

        # create an empty list to hold the row dictionaries
        new_rows = []

        # loop through the rows of the dataframe
        for index, row in df1.iterrows():
            
            # check if the value in the 'prophage_regions_0' column is a nested list
            if isinstance(row['prophage_regions_0'], list) and any(isinstance(sublist, list) for sublist in row['prophage_regions_0']):
                
                # if so, create a dictionary with the values for each new column
                new_row_dict = {}
                for i, sublist in enumerate(row['prophage_regions_0']):
                    new_col_name = f"prophage_regions_0_sublist_{i}"
                    new_row_dict[new_col_name] = sublist
                
                # append the new row dictionary to the list
                new_rows.append(new_row_dict)
            
            else:
                # if the value is not a nested list, append an empty dictionary to the list
                new_rows.append({})
                
        # create a new dataframe with the updated columns
        split_df = pd.concat([df1, pd.DataFrame(new_rows)], axis=1)
        print(split_df.to_string())

        split_df.drop('prophage_regions_0', axis=1,inplace=True)
        print('split_df')
        print(df1.to_string())
        for col in split_df.columns:
            # check if the column name contains the string "sublist"
            if "sublist" in col:
                # get the index of the column
                col_index = split_df.columns.get_loc(col)
                # insert a new empty column to the right of the current column
                split_df.insert(col_index+1, f"{col}_match?", pd.Series(dtype='object'))

        split_df['coord1'] = split_df['coord1'].astype(int)
        split_df['coord2'] = split_df['coord2'].astype(int)

        # Define a function to calculate the overlap between two ranges
        def overlap(a, b):
            return max(0, min(a[1], b[1]) - max(a[0], b[0]))

        # Define a function to check if there is more than 10% overlap between two ranges
        def overlap_gt_10_percent(a, b):
            return overlap(a, b) / (a[1] - a[0]) >= 0.1 and overlap(a, b) / (b[1] - b[0]) >= 0.1

        # Loop through each row in the dataframe
        for i, row in split_df.iterrows():
            # Loop through each column in the row
            for col in row.index:
                # Check if the column name contains the word 'sublist' and the cell value is a list
                if 'sublist' in col and isinstance(row[col], list):
                    sublist = row[col]
                    # Check if the sublist has at least 4 entries
                    if len(sublist) >= 4:
                        # Get the coordinates from the sublist
                        coord_sublist = sublist[2:4]
                        # Get the corresponding coordinates from the dataframe
                        coord_df = [row['coord1'], row['coord2']]
                        # Check if there is more than 10% overlap between the coordinates
                        if overlap_gt_10_percent(coord_sublist, coord_df):
                            # Write 'yes' into the column to the right of the current column
                            col_index = row.index.get_loc(col)
                            split_df.iloc[i, col_index+1] = 'yes'



        # print(split_df.to_string())
        split_df.to_csv(output_dir_path+f'{filee}_hittables_prophageregion_2.tsv', sep='\t', index=False)

        for col in split_df.columns:
        	print(col)
