"""
File: redshiftresult_cmp.py

Created on: 28/01/19
Author: CeSAM
"""

import argparse
import csv

# GLOBAL VARIABLES FOR COMPARISION PRECISION
FLOAT_PRECISION = 1e-12
REDSHIFT_PRECISION = 1e-12
MERIT_PRECISION = 1e-12
SIGME_PRECISION = 1e-12
INTGPROBA_PRECISION = 1e-12
GAUSSAMP_PRECISION = 1e-12
GAUSSAMPERR_PRECISION = 1e-12
GAUSSSIGMA_PRECISION = 1e-12
GAUSSSIGMAER_PRECISION = 1e-12

# ------------------------------------------------------------------------------


class ResultsComparator():
    def factory(type):
        if type == "redshift.csv":
            return Redshift():
        elif type == "candidatesresult.csv":
            return Candidate():
        elif type == "redshiftresult.csv":
            return RedshiftResult():
        elif type == "classificationresult.csv":
            return ClassificationResult():

# ------------------------------------------------------------------------------


class Redshift(ResultsComparator):
    def __init__():
        self.referencePath = "/home/ffauchier/amazed/"

        self.dict = {
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

# ------------------------------------------------------------------------------


class RedshiftResult(ResultsComparator):
    def __init__():
        self.referencePath = "/home/ffauchier/amazed/"

        self.dict = {
        "#Redshift": lambda x,y : abs(float(x)-float(y)) < REDSHIFT_PRECISION,
        "Merit": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
        "TemplateRatio": lambda x,y : x == y,
        "TemplateContinuum": lambda x,y : x == y,
        "method": lambda x,y : x == y,
        "sigma": lambda x,y : abs(float(x)-float(y)) < SIGME_PRECISION
        }

# ------------------------------------------------------------------------------


class Candidate(ResultsComparator):
    def __init__():
        self.referencePath = "/home/ffauchier/amazed/"

        self.dict = {
        "#rank": lambda x,y : int(x)-int(y),
        "redshift": lambda x,y : abs(float(x)-float(y)) < REDSHIFT_PRECISION,
        "intgProba": lambda x,y : abs(float(x)-float(y)) < INTGPROBA_PRECISION,
        "gaussAmp": lambda x,y : abs(float(x)-float(y)) < GAUSSAMP_PRECISION,
        "gaussAmpErr": lambda x,y : abs(float(x)-float(y)) < GAUSSAMPERR_PRECISION,
        "gaussSigma": lambda x,y : abs(float(x)-float(y)) < GAUSSSIGMA_PRECISION,
        "gaussSigmaErr": lambda x,y : abs(float(x)-float(y)) < GAUSSSIGMAER_PRECISION
        }

# ------------------------------------------------------------------------------


class ClassificationResult(ResultsComparator):
    def __init__():
        self.referencePath = "/home/ffauchier/amazed/"

        self.dict = {
        "#Type": lambda x,y : x == y,
        "Merit": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
        "EvidenceG": lambda x,y : abs(float(x)-float(y)) < FLOAT_PRECISION,
        "EvidenceS": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION,
        "EvidenceQ": lambda x,y : abs(float(x)-float(y)) < MERIT_PRECISION
        }

# ------------------------------------------------------------------------------


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

# ------------------------------------------------------------------------------


def read_csv(args):
    with open(args) as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        dict_list = {}
        header = next(reader)
        for line in reader:
            res = dict(zip(header, line))
            dict_list.update({res['#Spectrum']: res})
    return dict_list

# ------------------------------------------------------------------------------


def run(args):
    # TODO: definir le type en fonction du fichier donné
    # TODO: construire les path avec les args OUTPUT

    readerRef = read_csv(args.file1) # ref
    reader2 = read_csv(args.file2) # res

    obj = Redshift.factory(type)
    dictRes = obj.dict
    dict_cmp = obj.referencePath

    for spectrum in readerRef.keys():
        for k in readerRef[spectrum].keys():
            assert cmpDict[k](readerRef[spectrum][k],dictRes[spectrum][k])

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
