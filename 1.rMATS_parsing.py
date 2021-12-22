#! /data/home/Minkyu/anaconda3/bin/python

##### ENVIRONMENTS #####
import os, sys

##### PATHS #####
input_dic_Path = sys.argv[1]
pos = sys.argv[2]
out_path = sys.argv[3]
file_lst = os.listdir(input_dic_Path)
input_lst = [i for i in file_lst if i.endswith("JC.txt")]

##### USAGE #####
def Usage():
    print(
        "python 1.rMATS_parsing.py <<input file directory Path>> <<reference position file>> <<output file directory Path>>"
    )


##### SCRIPTS #####
def parsing(file, file2):
    txt = open(file, "r")
    outp = open(out_path + "select_" + file, "w")
    # number = open(file2, "r")
    header = txt.readline()
    outp.write(header)
    for line in txt:
        ln = line.strip().split("\t")
        if not ln[0].startswith("ID"):
            number = open(file2, "r")
            if len(header.strip().split("\t")) == 23:
                for ln2 in number:
                    if (
                        int(min(ln[5:11])) <= int(ln2) <= int(max(ln[5:11]))
                        and float(ln[19]) < 0.05
                    ):
                        outp.write(line)
                        break
            else:
                for ln2 in number:
                    if (
                        int(min(ln[5:13])) <= int(ln2) <= int(max(ln[5:13]))
                        and float(ln[21]) < 0.05
                    ):
                        outp.write(line)
                        break
    number.close()
    txt.close()
    outp.close()


def main():
    for ip_fi in input_lst:
        parsing(ip_fi, pos)


if __name__ == "__main__":
    main()
