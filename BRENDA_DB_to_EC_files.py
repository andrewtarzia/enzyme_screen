# script to convert full BRENDA_database.txt file into one file per EC.
import os

file_loc = '/home/atarzia/psp/brenda_details/'
os.chdir(file_loc)

file = 'brenda_download.txt'


switch = 0 
curr_list = []
with open(file, 'r') as f:
    for line in f:
        if 'ID\t' in line and switch == 1:
            # output
            with open(filename, 'w') as a:
                a.write('\n'.join(curr_list))
            switch = 0
            curr_list = []
            #import sys
            #sys.exit()
        if 'ID\t' in line and switch == 0:
            l = line.rstrip().split('\t')
            pt2 = l[1]
            if '(' in pt2:
                EC = pt2.split(" ")[0]
            else:
                EC = pt2
            filename = file.replace(".txt", "_"+EC.replace(".", "_")+".txt")
            switch = 1
            print(line, l, filename)
        if switch == 1:
            # implies writing to list
            curr_list.append(line)