# python
import os
import sys
import subprocess
import pipes
from Bio import SeqIO

# transform the rar into fasta
def decomplessing_file(file_name_row):
    folder_name = file_name_row.split('.')[0]
    cp1 = subprocess.Popen('rm -rf \'{}\''.format(folder_name), shell=True)
    cp1.wait()
    print 'rm -rf \'{}\''.format(folder_name)
    cp2 = subprocess.Popen('mkdir \'{}\''.format(folder_name), shell=True)
    cp2.wait()
    if file_name_row[-3:] == "rar":
        cp3 = subprocess.Popen('rar x \'{}\' \'{}\''.format(file_name_row, folder_name), shell=True)
    if file_name_row[-3:] == "zip":
        cp3 = subprocess.Popen('unzip -q -n \'{}\' -d \'{}\''.format(file_name_row, folder_name), shell=True)
    if file_name_row[-6:] == "tar.gz":
        cp3 = subprocess.Popen('tar -zcxf \'{}\' -C \'{}\''.format(file_name_row, folder_name), shell=True)
    cp3.wait()
    return folder_name

def parse_ab1(ab1_file):
    handle = open(ab1_file, "rb")
    for record in SeqIO.parse(handle, "abi"):
        return record.id
    handle.close()

def parse_seq(seq_file):
    handle = open(seq_file, "r")
    return handle.read().replace('^M$', '')
    handle.close()

file_name_row = sys.argv[1]
folder_name = decomplessing_file(file_name_row)

fas_file = open(folder_name + ".fas", "w")

for (path, subdirs, files) in os.walk(folder_name):
    for name in files:
        file_path = os.path.join(path, name)
        if name[-4:] == ".abl" or name[-4:] == ".ab1" or name[-4:] == ".abi":
            parse_ab1(file_path)
        else:
            fas_file.write('>' + name[0:-4] + '\n')
            fas_file.write(parse_seq(file_path) + '\n')
            subprocess.Popen('rm \'{}\''.format(file_path), shell=True)

fas_file.close()

#os.system('rm -rf ./' + folder_name)
