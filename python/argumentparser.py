import argparse

parser = argparse.ArgumentParser(description='CPF-redshift client tool.')

parser.add_argument('--parameters', '-p', dest='parameters_file', metavar='FILE', type=str,
                    default='parameters.json', help='Parameters file')
parser.add_argument('--datapath', '-d', dest='data_path', metavar='DIR', type=str,
                    help='Data root-path')
parser.add_argument('--input', '-i', dest='input_file', metavar='FILE', type=str,
                    default='input.spectrumlist',
                    help='Input file containing the spectrums list file to process')
parser.add_argument('--error', '-e', dest='error_file', metavar='FILE', type=str,
                    help='Path to a spectrum file containing error data to be used if a specific fits file was used as input')
parser.add_argument('--spectrum_dir', '-f', dest='spectrum_dir', metavar='DIR', type=str,
                    default='spectrum',
                    help='Specify directory in which input spectrum files are stored')
parser.add_argument('--template_dir', '-t', dest='template_dir', metavar='DIR', type=str,
                    default='./calibration/templates/ExtendedTemplatesJan2017_v3',
                    help='Specify directory in which input templates files are stored')
parser.add_argument('--linecatalog', '-y', dest='linecatalog', metavar='FILE', type=str,
                    default='linecatalog.txt',
                    help='Path to the rest lines catalog file')
parser.add_argument('--calibration_dir', '-k', dest='calibration_dir', metavar='DIR', type=str,
                    default='./calibration',
                    help='Specify directory in which calibration files are stored')
parser.add_argument('--zclassifier_dir', '-z', dest='zclassifier_dir', metavar='DIR', type=str,
                    default='',
                    help='Specify directory in which zClassifier files are stored')

