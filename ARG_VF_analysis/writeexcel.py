# compile the .tab files from each sample into a singular excel sheet with tabs for each sample, which can be consolidated into a single sheet using the sortexcel.py script

import os
import csv
import xlsxwriter

# Set the directory path to where your .tab files are located for a singular database
tab_file_dir = '/path/to/results'

# Create an excel file, and set the output Excel file path and name
excel_file_path = '/path/to/excel/compiled_data.xlsx'

# Create an Excel workbook object
workbook = xlsxwriter.Workbook(excel_file_path)

# Loop through all .tab files in the specified directory
for filename in os.listdir(tab_file_dir):
    if filename.endswith('.tab'):
        # Get the file name without the extension
        sheet_name = os.path.splitext(filename)[0]
        
        # Create a new worksheet in the Excel workbook
        worksheet = workbook.add_worksheet(sheet_name)
        
        # Open the .tab file and read its contents
        with open(os.path.join(tab_file_dir, filename), 'r') as tab_file:
            reader = csv.reader(tab_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                for col_num, value in enumerate(row):
                    worksheet.write(row_num, col_num, value)

# Close the Excel workbook
workbook.close()
