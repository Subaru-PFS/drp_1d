"""
File: output_compare.py

Created on: 28/01/19
Author: CeSAM
"""

import os
import csv
import json
import argparse
import numpy as np
from pprint import pprint
from astropy.io import fits
from xml.etree import ElementTree

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


def str_cmp(x, y):
    return str(x) == str(y)


def true_cmp(x, y):
    """
    Temporary function used in ClassificationResult() to return True, due to no
    value of EvidenceQ in the csv.

    :return: True
    """
    return True


def float_cmp(precision):
    """
    Check the precision of the difference between two float.

    :param precision: The precision of the difference.
    :return: Return True if the precision is lower than the precision given in
    argument.
    """
    return lambda x, y: abs(float(x)-float(y)) < precision


def int_cmp(x, y):
    """
    Check that the two integers are equal.

    :return: Return True if the two integers are equal.
    """
    return int(x) == int(y)


def array_cmp(element_cmp):
    """
    Check that the two arrays are equal.

    Arrays are compared element-wise with element_cmp

    :return: Return True if the two arrays are equal.
    """
    return lambda X, Y: np.all([element_cmp(x, y) for (x, y) in zip(X, Y)])


def read_spectrumlist(spectrum_list_file):
    """
    Read the spectrum list file and get all the process id for each spectrum.

    :param spectrum_list_file: The "spectrum.list" file.
    :return: Return the list of all process_id.
    """
    spectrumList = []
    for e in open(spectrum_list_file, 'r'):
        if e and not e.startswith('#'):
            spectrumList.append(e.split()[2].strip())

    return spectrumList


class ResultsComparator:
    def compare(self, path1, path2):
        pass


class FitsComparator(ResultsComparator):
    """Compare two FITS files"""

    @staticmethod
    def compare_columns(hdu, name, comparator, col1, col2):
        result = []
        for i, (v1, v2) in enumerate(zip(col1, col2)):
            if not comparator(v1, v2):
                if len(result) < 3:
                    result.append('HDU {}, column {}, line {} : '
                                  '{} != {}'.format(hdu, name, i, v1, v2))
                else:
                    result.append('...')
                    break
        return result

    def compare_hdu(self, name, hdul1, hdul2):
        result = []
        hdu1 = hdul1[name].data
        hdu2 = hdul2[name].data
        for col, comparator in self.tests[name].items():
            if col == '__key':
                continue
            r = self.compare_columns(name, col, comparator,
                                     hdu1[col], hdu2[col])
            if r:
                result.extend(r)
        return result

    def compare(self, path1, path2):
        result = []
        with fits.open(path1) as hdul1, fits.open(path2) as hdul2:
            for name in self.tests.keys():
                r = self.compare_hdu(name, hdul1, hdul2)
                if r:
                    result.extend(r)
        return result


class SpeClassificationCatalogComparator(FitsComparator):
    """Compare two Euclid SpeClassificationCatalog"""
    tests = {
        'SPE_CLASSIFICATION_CAT': {'__key': 'OBJECT_ID',
                                   'OBJECT_ID': int_cmp,
                                   'SPE_STAR_PROB': float_cmp(FLOAT_PRECISION),
                                   'SPE_GAL_PROB': float_cmp(FLOAT_PRECISION),
                                   'SPE_QSO_PROB': float_cmp(FLOAT_PRECISION)}}


class SpeZCatalogComparator(FitsComparator):
    """Compare two Euclid SpeZCatalog"""
    tests = {
        'SPE_REDSHIFT_CAT': {'__key': 'SPE_RANK',
                             'OBJECT_ID': int_cmp,
                             'SPE_RANK': int_cmp,
                             'SPE_Z':  float_cmp(FLOAT_PRECISION),
                             'SPE_Z_REL': float_cmp(FLOAT_PRECISION)}}


class SpePdfCatalogComparator(FitsComparator):
    """Compare two Euclid SpePdfCatalog"""
    tests = {
        'SPE_PDF_CAT': {'__key': 'OBJECT_ID',
                        'OBJECT_ID': int_cmp,
                        'SPE_PDF_CLASS': str_cmp,
                        'SPE_PDF': array_cmp(float_cmp(FLOAT_PRECISION))}}


class SpeLineCatalogComparator(FitsComparator):
    """Compare two Euclid SpeLineCatalog"""
    tests = {
        'SPE_LINE_FEATURES_CAT': {
            '__key': 'OBJECT_ID',
            'OBJECT_ID': int_cmp,
            'SPE_RANK': int_cmp,
            'SPE_LINE_ID': int_cmp,
            'SPE_LINE_NAME': str_cmp,
            'SPE_LINE_CENTRAL_WL': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_CENTRAL_WL_ERR': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_FLUX': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_FLUX_ERR': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_EW': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_EW_ERR': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_FWHM': float_cmp(FLOAT_PRECISION),
            'SPE_LINE_FWHM_ERR': float_cmp(FLOAT_PRECISION)}}


class SpeAbsorptionCatalogComparator(FitsComparator):
    """Compare two Euclid SpeAbsorptionCatalog"""
    tests = {
        'SPE_SPECTRAL_FEATURES_CAT': {
            '__key': 'OBJECT_ID',
            'OBJECT_ID': int_cmp,
            'SPE_RANK': int_cmp,
            'SPE_ABSORPTION_NAME': str_cmp,
            'SPE_ABSORPTION': float_cmp(FLOAT_PRECISION),
            'SPE_ABSORPTION_ERR': float_cmp(FLOAT_PRECISION)}}


class SpeParametersCatalogComparator(FitsComparator):
    """Compare two Euclid SpeParametersCatalog"""
    tests = {
        'SPE_REST_FRAME_PARAMETERS_CAT': {
            '__key': 'OBJECT_ID',
            'OBJECT_ID': int_cmp,
            'SPE_RANK': int_cmp,
            'SPE_SFR': float_cmp(FLOAT_PRECISION)}}


def read_xml(xml):

    xml_tree = ElementTree.parse(xml)
    root = xml_tree.getroot()

    type = ['SpeParametersCatalog', 'SpeAbsorptionCatalog', 'SpeLineCatalog',
            'SpePdfCatalog', 'SpeZCatalog', 'SpeClassificationCatalog']

    product_dict = {}
    for element in type:
        basename = root.findall('./Data/{}/DataContainer/'
                                'FileName'.format(element))
        filename = basename[0].text
        product_dict[element] = filename

    return product_dict


class EuclidXMLComparator(ResultsComparator):

    def compare(self, xml1, xml2):
        """
        Compares XML data product files.

        :param xml1: Path to the reference XML.
        :param xml2: Path to the XML being compared.
        :return: Return the result of the comparison.
        """
        map_dict = {
            'SpeParametersCatalog': SpeParametersCatalogComparator,
            'SpeAbsorptionCatalog': SpeAbsorptionCatalogComparator,
            'SpeLineCatalog': SpeLineCatalogComparator,
            'SpePdfCatalog': SpePdfCatalogComparator,
            'SpeZCatalog': SpeZCatalogComparator,
            'SpeClassificationCatalog': SpeClassificationCatalogComparator
        }
        result = []
        dict1 = read_xml(xml1)
        dict2 = read_xml(xml2)

        for key, value in dict1.items():
            try:
                r = map_dict[key]().compare(os.path.join(os.path.dirname(xml1),
                                                         value),
                                            os.path.join(os.path.dirname(xml2),
                                                         dict2.get(key)))
                if r:
                    result.append(r)
            except Exception as e:
                result.append('{}/{} : {}'.format(value, key, e))

        return result


class OutputDirComparator(ResultsComparator):

    def compare(self, path1, path2):
        """
        Calls a comparator for each type of file by iterating if them.

        :param path1: Path to the reference file.
        :param path2: Path to the file being compared.
        :return: Return the result of the comparison.
        """
        global_tests = [(Redshift, []),
                        (JSONComparator, ['parameters.json']),
                        (JSONComparator, ['config.json', ['output_folder']]),
                        (LinecatalogComparator, [])]

        result = []
        for method, args in global_tests:
            r = method(*args).compare(path1, path2)
            if r:
                result.append(r)

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


class JSONComparator(ResultsComparator):

    def __init__(self, filename, ignore_keys=None):
        self.filename = filename
        self.ignore_keys = ignore_keys if ignore_keys is not None else []

    def compare(self, path1, path2):
        result = []

        with open(os.path.join(path1, self.filename), 'r') as f:
            j1 = json.load(f)
        with open(os.path.join(path2, self.filename), 'r') as f:
            j2 = json.load(f)

        s1 = set(j1.keys())
        s2 = set(j2.keys())
        if s1 != s2:
            result.append('{} : keys mismatch [{}]'.format(
                self.filename, s1 ^ s2))
            self.ignore_keys.extend(s1 ^ s2)

        for k in j1.keys():
            if k in self.ignore_keys:
                continue
            if j1[k] != j2[k]:
                result.append('{}/{} : value mismatch [{}] / [{}]'.format(
                    self.filename, k, j1[k], j2[k]))
        return result


class CSVComparator(ResultsComparator):
    tests = {}
    filename = ""
    key = None
    header_mark = None

    def _read_csv(self, args):
        """
        Create dictionnaries with the content of the CSV files given in
        parameters.

        :param args: Path to the csv to read.
        :return: Return the dictionnaries lists
        """
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
        """
        :param path1: Path to the first file to compare.
        :param path2: Path to the second file to compare.
        :return: Return list of all differences
        """
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
        "#rank": int_cmp,
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


class LinecatalogComparator(CSVComparator):
    filename = "linecatalog.txt"
    key = "name"
    header_mark = "#lambda"
    tests = {
        "#lambda": float_cmp(FLOAT_PRECISION),
        "name": str_cmp,
        "type": str_cmp,
        "force": str_cmp,
        "profile": str_cmp,
        "amp_group": str_cmp,
        "nominal_ampl": str_cmp,
        "vel_group": str_cmp,
    }


class zPDFComparator(ResultsComparator):
    filename = 'logposterior.logMargP_Z_data.csv'

    @staticmethod
    def _read_pdf(path):
        """
        :param path: Path to the csv to read.
        :return: Return the list of all differences.
        """
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
        """
        :param path1: Path to the first file to compare.
        :param path2: Path to the second file to compare.
        :return: Return the list of all differences.
        """
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


def compare_amazed(args):
    r = OutputDirComparator().compare(args.firstdir,
                                      args.seconddir)
    return r


def compare_euclid(args):
    r = EuclidXMLComparator().compare(args.firstxml,
                                      args.secondxml)
    return r


def main():
    parser = argparse.ArgumentParser(
        description='AMAZED results comparison tool'
    )
    subparsers = parser.add_subparsers(dest='command')

    amazed_parser = subparsers.add_parser('amazed',
                                          help='Compare AMAZED results')
    amazed_parser.add_argument('referencedir',
                               type=str,
                               help='Path to the first output directory.')
    amazed_parser.add_argument('outputdir',
                               type=str,
                               help='Path to the second output directory.')
    amazed_parser.set_defaults(handler=compare_amazed)

    euclid_parser = subparsers.add_parser('euclid',
                                          help='Compare EUCLID results')
    euclid_parser.add_argument('firstxml',
                               type=str,
                               help='Path to the first XML file.')
    euclid_parser.add_argument('secondxml',
                               type=str,
                               help='Path to the second XML file.')
    euclid_parser.set_defaults(handler=compare_euclid)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return

    r = args.handler(args)

    if r:
        pprint(r)


if __name__ == '__main__':
    main()
