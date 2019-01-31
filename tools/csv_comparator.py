"""
File: redshiftresult_cmp.py

Created on: 28/01/19
Author: CeSAM
"""

import os
import csv
import argparse

# GLOBAL VARIABLES FOR COMPARISION PRECISION
FLOAT_PRECISION = 1e-12
# MERIT_PRECISION = 1e-12
# SIGMA_PRECISION = 1e-12
# REDSHIFT_PRECISION = 1e-12
# GAUSSAMP_PRECISION = 1e-12
# INTGPROBA_PRECISION = 1e-12
# GAUSSSIGMA_PRECISION = 1e-12
# GAUSSAMPERR_PRECISION = 1e-12
# GAUSSSIGMAER_PRECISION = 1e-12


strCmp = str.__eq__
trueCmp = lambda x, y: True

def float_cmp(precision):
    return lambda x,y: abs(float(x)-float(y)) < precision


def int_cmp():
    return lambda x,y: int(x) == int(y)


def read_spectrumlist(spectrumListFile):
    spectrumList = []
    for e in open(spectrumListFile, 'r'):
        spectrumList.append(e.split(' ')[2].strip())

    return spectrumList


class ResultsComparator():


    def compare(self, path1, path2):
        pass


class OutputDirComparator(ResultsComparator):


    def compare(self, path1, path2, spectrumListFile):

        spectrumList = read_spectrumlist(spectrumListFile)
        result = Redshift().compare(path1, path2)
        for spectrum in spectrumList:
            for method in [CandidateResult,
                           RedshiftResult,
                           ClassificationResult]:
                r = method().compare(os.path.join(path1, spectrum),
                                     os.path.join(path2, spectrum))
                if r:
                    result.append(r)
        return result


class CSVComparator():
    tests = {}
    filename = ""
    key = None
    headerMark = None

    def _read_csv(self, args):
        csvFile = open(args)
        for l in csvFile:
            if l.startswith(self.headerMark):
                header = l.split()
                break
        reader = csv.reader(csvFile, delimiter='\t')
        dict_list = {}

        for line in reader:
            if not line:
                #strip empty lines
                continue
            res = dict(zip(header, line))
            dict_list.update({res[self.key]: res})
        return dict_list


    def compare(self, path1, path2):

        csv1 = self._read_csv(os.path.join(path1, self.filename))
        csv2 = self._read_csv(os.path.join(path2, self.filename))
        result = []


        for spectrum in csv1.keys():
            for k in csv1[spectrum].keys():
                try:
                    if not self.tests[k](csv1[spectrum][k], csv2[spectrum][k]):
                        result.append("{} : [{}] != [{}]".format(k, csv1[spectrum][k], csv2[spectrum][k]))
                except Exception as e :
                    print("error on key {}".format(k))
                    raise
        return result


class Redshift(CSVComparator):
    filename = "redshift.csv"
    key = '#Spectrum'
    headerMark = '#Spectrum'
    tests = {
    "#Spectrum": strCmp,
    "ProcessingID": strCmp,
    "Redshift": float_cmp(FLOAT_PRECISION),
    "Merit": float_cmp(FLOAT_PRECISION),
    "Template": strCmp,
    "Method": strCmp,
    "Deltaz": float_cmp(FLOAT_PRECISION),
    "Reliability": strCmp,
    "snrHa": float_cmp(FLOAT_PRECISION),
    "lfHa": float_cmp(FLOAT_PRECISION),
    "snrOII": float_cmp(FLOAT_PRECISION),
    "lfOII": float_cmp(FLOAT_PRECISION),
    "Type": strCmp
    }


class RedshiftResult(CSVComparator):
    filename = "redshiftresult.csv"
    key = '#Redshifts'
    headerMark = '#Redshifts'
    tests = {
    "#Redshifts": float_cmp(FLOAT_PRECISION),
    "Merit": float_cmp(FLOAT_PRECISION),
    "TemplateRatio": strCmp,
    "TemplateContinuum": strCmp,
    "method": strCmp,
    "sigma": float_cmp(FLOAT_PRECISION)
    }


class CandidateResult(CSVComparator):
    filename = "candidatesresult.csv"
    key = '#rank'
    headerMark = '#rank'
    tests = {
    "#rank": int_cmp(),
    "redshift": float_cmp(FLOAT_PRECISION),
    "intgProba": float_cmp(FLOAT_PRECISION),
    "gaussAmp": float_cmp(FLOAT_PRECISION),
    "gaussAmpErr": float_cmp(FLOAT_PRECISION),
    "gaussSigma": float_cmp(FLOAT_PRECISION),
    "gaussSigmaErr": float_cmp(FLOAT_PRECISION)
    }


class ClassificationResult(CSVComparator):
    filename = "classificationresult.csv"
    key = "#Type"
    headerMark = "#Type"
    tests = {
    "#Type": strCmp,
    "Merit": float_cmp(FLOAT_PRECISION),
    "EvidenceG": float_cmp(FLOAT_PRECISION),
    "EvidenceS": float_cmp(FLOAT_PRECISION),
    "EvidenceQ": trueCmp
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('referenceDir',
                        type=str,
                        help='Path to the output reference directory.')
    parser.add_argument('outputDir',
                        type=str,
                        help='Path to the output tests directory.')
    parser.add_argument('spectrumListFile',
                        type=str,
                        help='Path to spectrum.list file.')
    args = parser.parse_args()

    r = OutputDirComparator().compare(args.referenceDir,
                                      args.outputDir,
                                      args.spectrumListFile)
    if r: print(r)


if __name__ == '__main__':
    main()
