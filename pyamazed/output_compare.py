"""
File: output_compare.py

Created on: 28/01/19
Author: CeSAM
"""

import os
import csv
import argparse
import numpy as np
from pprint import pprint

# GLOBAL VARIABLES FOR COMPARISION PRECISION
FLOAT_PRECISION = 1e-12
PDF_PRECISION = 1e-7
# MERIT_PRECISION = 1e-12
# SIGMA_PRECISION = 1e-12
# REDSHIFT_PRECISION = 1e-12
# GAUSSAMP_PRECISION = 1e-12
# INTGPROBA_PRECISION = 1e-12
# GAUSSSIGMA_PRECISION = 1e-12
# GAUSSAMPERR_PRECISION = 1e-12
# GAUSSSIGMAER_PRECISION = 1e-12


str_cmp = str.__eq__


def true_cmp(x, y):
    return True


def float_cmp(precision):
    return lambda x, y: abs(float(x)-float(y)) < precision


def int_cmp():
    return lambda x, y: int(x) == int(y)


def read_spectrumlist(spectrumListFile):
    spectrumList = []
    for e in open(spectrumListFile, 'r'):
        if e and not e.startswith('#'):
            spectrumList.append(e.split()[2].strip())

    return spectrumList


class ResultsComparator():

    def compare(self, path1, path2):
        pass


class OutputDirComparator(ResultsComparator):

    def compare(self, path1, path2):

        result = Redshift().compare(path1, path2)
        spectrumlist = read_spectrumlist(os.path.join(path1,
                                                      'input.spectrumlist'))
        for spectrum in spectrumlist:
            for method in [CandidateResult,
                           RedshiftResult,
                           ClassificationResult,
                           zPDFComparator]:
                try:
                    r = method().compare(os.path.join(path1, spectrum),
                                         os.path.join(path2, spectrum))
                    if r:
                        result.append(r)
                except Exception as e:
                    result.append('{}/{} : {}'.format(spectrum,
                                                      method.__name__, e))
        return result


class CSVComparator(ResultsComparator):
    tests = {}
    filename = ""
    key = None
    header_mark = None

    def _read_csv(self, args):
        csvfile = open(args)
        for l in csvfile:
            if l.startswith(self.header_mark):
                header = l.split()
                break
        reader = csv.reader(csvfile, delimiter='\t')
        dict_list = {}

        for line in reader:
            if not line:
                # strip empty lines
                continue
            res = dict(zip(header, line))
            dict_list.update({res[self.key]: res})
        return dict_list

    def compare(self, path1, path2):

        try:
            csv1 = self._read_csv(os.path.join(path1, self.filename))
            csv2 = self._read_csv(os.path.join(path2, self.filename))
        except Exception as e:
            return ["can't read file : {}".format(e)]

        result = []

        for spectrum in csv1.keys():
            for k in csv1[spectrum].keys():
                try:
                    if not self.tests[k](csv1[spectrum][k], csv2[spectrum][k]):
                        result.append('{}/{} : '
                                      '[{}] != [{}]'.format(spectrum, k,
                                                            csv1[spectrum][k],
                                                            csv2[spectrum][k]))
                except Exception:
                    print("error on key {}".format(k))
                    raise
        return result


class Redshift(CSVComparator):
    filename = "redshift.csv"
    key = '#Spectrum'
    header_mark = '#Spectrum'
    tests = {
        "#Spectrum": str_cmp,
        "ProcessingID": str_cmp,
        "Redshift": float_cmp(FLOAT_PRECISION),
        "Merit": float_cmp(FLOAT_PRECISION),
        "Template": str_cmp,
        "Method": str_cmp,
        "Deltaz": float_cmp(FLOAT_PRECISION),
        "Reliability": str_cmp,
        "snrHa": float_cmp(FLOAT_PRECISION),
        "lfHa": float_cmp(FLOAT_PRECISION),
        "snrOII": float_cmp(FLOAT_PRECISION),
        "lfOII": float_cmp(FLOAT_PRECISION),
        "Type": str_cmp
    }


class RedshiftResult(CSVComparator):
    filename = "redshiftresult.csv"
    key = '#Redshifts'
    header_mark = '#Redshifts'
    tests = {
        "#Redshifts": float_cmp(FLOAT_PRECISION),
        "Merit": float_cmp(FLOAT_PRECISION),
        "TemplateRatio": str_cmp,
        "TemplateContinuum": str_cmp,
        "method": str_cmp,
        "sigma": float_cmp(FLOAT_PRECISION)
    }


class CandidateResult(CSVComparator):
    filename = "candidatesresult.csv"
    key = '#rank'
    header_mark = '#rank'
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
    header_mark = "#Type"
    tests = {
        "#Type": str_cmp,
        "Merit": float_cmp(FLOAT_PRECISION),
        "EvidenceG": float_cmp(FLOAT_PRECISION),
        "EvidenceS": float_cmp(FLOAT_PRECISION),
        "EvidenceQ": true_cmp
    }


class zPDFComparator(ResultsComparator):
    filename = 'logposterior.logMargP_Z_data.csv'

    @staticmethod
    def _read_pdf(path):
        with open(os.path.join(path, 'zPDF',
                               zPDFComparator.filename), "r") as f:
            header = f.readline()
            assert header.startswith('#z_tested')
            evidence = float(f.readline().split('=')[1])
            pdf = []
            for l in f:
                pdf.append(l.split())
        return evidence, np.array(pdf, dtype=float)

    def compare(self, path1, path2):
        result = []
        ev1, pdf1 = self._read_pdf(path1)
        ev2, pdf2 = self._read_pdf(path2)
        if not float_cmp(FLOAT_PRECISION)(ev1, ev2):
            result.append('zPDF EvidenceLog: [{}] != [{}]'.format(ev1, ev2))
        i1 = np.trapz(np.exp(pdf1[:, 1]), x=pdf1[:, 0])
        i2 = np.trapz(np.exp(pdf2[:, 1]), x=pdf2[:, 0])
        if not float_cmp(FLOAT_PRECISION)(i1, i2):
            result.append('zPDF : [{}] != [{}]'.format(i1, i2))
        if not float_cmp(PDF_PRECISION)(i1, 1.0):
            result.append('PDF integrate should be 1.0+/{}'
                          '[{}]'.format(PDF_PRECISION, i1))
        return result


def main():
    parser = argparse.ArgumentParser(description='AMAZED results comparison tool')
    parser.add_argument('referencedir',
                        type=str,
                        help='Path to the output reference directory.')
    parser.add_argument('outputdir',
                        type=str,
                        help='Path to the output tests directory.')
    args = parser.parse_args()

    r = OutputDirComparator().compare(args.referencedir,
                                      args.outputdir)
    if r:
        pprint(r)


if __name__ == '__main__':
    main()
