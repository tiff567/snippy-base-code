import pandas as pd
import pyfaidx
import math
import csv
import numpy as np
import sys

from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

from UtilityFunctions import sent_email_notification

# procedure step 1: obtain csv
fields = ['Chr', 'PamNum', 'OptimalNum', 'gRNA']
df = pd.read_csv('Output.csv', skipinitialspace=True, usecols=fields)

rowMaxNum = df.shape[0]     # get the row num of data
gRNA_Str = ""
PamNum = 0
optimalNum = 0
rownum = 0

gRNA_List = []

while (rownum <= rowMaxNum - 1):
    Chr = df.loc[rownum, 'Chr']

    if not math.isnan(Chr):  # check not nan string
        PamNum = int(df.loc[rownum, 'PamNum'])
        optimalNum = int(df.loc[rownum, 'OptimalNum'])
        gRNA_Str = str(df.loc[rownum, 'gRNA'])

        if (gRNA_Str != ""):   # not an empty string
            gRNA_StrList = gRNA_Str.split()

            for gRNA in gRNA_StrList:
                gRNA_Seq = Seq(gRNA)
                if(str(gRNA_Seq) != "nan"):
                    print ("Sequence to be searched: " + str(rownum) + "\t" + str(gRNA_Seq))
                    gRNA_List.append(str(gRNA_Seq))
    rownum += 1

print("Searching NCBI BLAST sequence database. Please wait ....")



result_handle_list = []

length = len(gRNA_List)
for i in range(length):
    result_handle = NCBIWWW.qblast("blastn", "nr", gRNA_List[i])
    step = int(((i+1)/length) * 100)
    print(f'sequence to be searched: {gRNA_List[i]:s} {step:3d} % Done!...')
    result_handle_list.append(result_handle)

print("Searching was Finished!")
# send email notification to the users
#sent_email_notification("tiffashheart@gmail.com", "tfang_99@yahoo.com", "BLAST search done", "BlastDoneNote.txt")

# save the searching result to local XML file
for i in range(length):
    gRNA_Str = gRNA_List[i]
    fileName = gRNA_Str + "_blast.xml"
    with open(fileName, "w") as out_handle:
        result_handle = result_handle_list[i]
        out_handle.write(result_handle.read())

print("Results has been written to XML file")

# retrieve the searching results from saved file
for i in range(length):
    gRNA_Str = gRNA_List[i]
    fileName = gRNA_Str + "_blast.xml"
    result_handle = open(fileName)
    blast_record = NCBIXML.read(result_handle)

    E_VALUE_THRESH = 4.0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("e value:", hsp.expect)
            if hsp.expect < E_VALUE_THRESH:
                print("\nAlignment for " + gRNA_Str)
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("score:", hsp.score)
                print("gaps:", hsp.gaps)
                print("e value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")




