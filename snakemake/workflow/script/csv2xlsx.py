import sys
import csv
from openpyxl import Workbook

# for frontend only
csvfilename = sys.argv[1]
excelfilename = sys.argv[2]

wb = Workbook()
ws = wb.active

with open(csvfilename) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        ws.append(row)
wb.save(excelfilename)




#Reference: https://www.ashwoodcomputer.com/uncategorized/short-python-code-to-create-an-excel-file-from-a-csv-file/
