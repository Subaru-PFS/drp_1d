"""
File: redshiftresult_cmp.py

Created on: 28/01/19
Author: CeSAM
"""

import argparse
import csv

# GLOBAL VARIABLES FOR PRECISION COMPARISION
FLOAT_PRECISION = 1e-12
REDSHIFT_PRECISION = 1e-12
MERIT_PRECISION = 1e-12
SIGME_PRECISION = 1e-12
INTGPROBA_PRECISION = 1e-12
GAUSSAMP_PRECISION = 1e-12
GAUSSAMPERR_PRECISION = 1e-12
GAUSSSIGMA_PRECISION = 1e-12
GAUSSSIGMAER_PRECISION = 1e-12

def switch(arg):
    switcher = {
        1 = redshift_cmp,
        2 = redshiftresult_cmp,
        3 = candidatesresult_cmp,
        4 = classificationresult_cmp
        }
    print switcher.get(arg)

dl = [redshift_cmp,
      redshiftresult_cmp,
      candidatesresult_cmp,
      classificationresult_cmp]

redshift_cmp = {
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
    # "Type": lambda x,y : x == y
    "Type": True
}

redshiftresult_cmp = {
    "#Redshift": lambda x,y : abs(float(x)-float(y)) < REDSHIFT_PRECISION,
    "Merit": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
    "TemplateRatio": lambda x,y : x == y,
    "TemplateContinuum": lambda x,y : x == y,
    "method": lambda x,y : x == y,
    "sigma": lambda x,y : abs(float(x)-float(y)) < SIGME_PRECISION
}

candidatesresult_cmp = {
    "#rank": lambda x,y : int(x)-int(y),
    "redshift": lambda x,y : abs(float(x)-float(y)) < REDSHIFT_PRECISION,
    "intgProba": lambda x,y : abs(float(x)-float(y)) < INTGPROBA_PRECISION,
    "gaussAmp": lambda x,y : abs(float(x)-float(y)) < GAUSSAMP_PRECISION,
    "gaussAmpErr": lambda x,y : abs(float(x)-float(y)) < GAUSSAMPERR_PRECISION,
    "gaussSigma": lambda x,y : abs(float(x)-float(y)) < GAUSSSIGMA_PRECISION,
    "gaussSigmaErr": lambda x,y : abs(float(x)-float(y)) < GAUSSSIGMAER_PRECISION
}

classificationresult_cmp = {
    "#Type": lambda x,y : x == y,
    "Merit": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
    "EvidenceG": lambda x,y : abs(float(x)-float(y)) < FLOAT_PRECISION,
    "EvidenceS": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
    "EvidenceQ": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir',
                        type=str,
                        help='Path to the "output" directory.')
    parser.add_argument('spectrum_list',
                        type=str,
                        help='Path to spectrum.list file.')
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


    for spectrum in reader1.keys():
        for k in reader1[spectrum].keys():
            assert cmpDict[k](reader1[spectrum][k],reader2[spectrum][k])


if __name__ == '__main__':
    main()
