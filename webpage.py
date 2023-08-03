import csv
import webbrowser
import requests
import os


path = '/Users/pimswart/alignment_hittables/P10.bd320ESW563013-Alignment-HitTable.csv'
filename = os.path.basename(path)

try:
    with open(path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader) # skip header row
        unique_strings = set()
        filename_withoutcsv = filename[:-4]
        os.mkdir(filename_withoutcsv)

        for row in csv_reader:
            if len(row) >= 2: # check if row has at least 2 columns
                unique_strings.add(row[1])
            else:
                print(f"Invalid row format: {row}") # handle invalid rows

        for unique_string in unique_strings:
            url = f"https://phaster.ca/jobs/{unique_string}/summary.txt"
            r = requests.get(url, allow_redirects=True)
            open(f'/Users/pimswart/{filename_withoutcsv}/{unique_string}_summary.txt', 'wb').write(r.content)
except Exception as e:
    print(f"An error occurred: {e}")
