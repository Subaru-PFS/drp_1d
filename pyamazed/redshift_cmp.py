"""
File: redshift_cmp.py

Created on: 22/01/19
Author: CeSAM
"""

import argparse
import csv

REDSHIFT_PRECISION = 1e-12
MERIT_PRECISION = 1e-12

cmpDict = {
    "#Spectrum": lambda x,y : x == y,
    "ProcessingID": lambda x,y : True,
    "Redshift": lambda x,y : abs(float(x)-float(y)) < REDSHIFT_PRECISION,
    "Merit": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
    "Template": lambda x,y : x == y,
    "Method": lambda x,y : x == y,
    "Deltaz": lambda x,y : int(x) == int(y),
    "Reliability": lambda x,y : x == y,
    "snrHa": lambda x,y : float(x) == float(y),
    "lfHa": lambda x,y : float(x) == float(y),
    "snrOII": lambda x,y : float(x) == float(y),
    "lfOII": lambda x,y : float(x) == float(y),
    "Type": lambda x,y : x == y
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1',
                        type=str,
                        help='Enter the path to the first file to compare')
    parser.add_argument('file2',
                        type=str,
                        help='Enter the path to the second file to compare')
    args = parser.parse_args()
    return run(args)


def read_csv(args):
    with open(args) as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        dict_list = {}
        header = next(reader)
        for line in reader:
            res = dict(zip(header, line))
            dict_list.update({res['#Spectrum']: res})
    return dict_list


def run(args):
    reader1 = read_csv(args.file1)
    reader2 = read_csv(args.file2)

    for spectrum in reader1.keys():
        for k in reader1[spectrum].keys():
            assert cmpDict[k](reader1[spectrum][k],reader2[spectrum][k])


if __name__ == '__main__':
    main()
