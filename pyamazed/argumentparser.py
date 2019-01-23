import argparse
from .redshift import CLog, get_version


def log_level(lvl):
    levels = {
        'critical': CLog.nLevel_Critical,
        'error': CLog.nLevel_Error,
        'warning': CLog.nLevel_Warning,
        'info': CLog.nLevel_Info,
        'detail': CLog.nLevel_Detail,
        'debug': CLog.nLevel_Debug,
        'none': CLog.nLevel_None
    }
    if lvl.lower() not in levels:
        print('Warning: unknown log lever "%s"' % lvl)
        return CLog.nLevel_Info
    else:
        return levels[lvl.lower()]


parser = argparse.ArgumentParser(description='CPF-redshift client tool.')

parser.add_argument('--parameters', '-p', dest='parameters_file',
                    metavar='FILE', type=str, help='Parameters file')
parser.add_argument('--datapath', '-d', dest='data_path', metavar='DIR',
                    type=str, help='Data root-path')
parser.add_argument('--output', '-o', dest='output_folder', metavar='DIR',
                    type=str,
                    help='Directory where all generated files are '
                    'going to be stored')
parser.add_argument('--input', '-i', dest='input_file', metavar='FILE',
                    type=str,
                    help='Input file containing the spectrums list '
                    'file to process')
parser.add_argument('--error', '-e', dest='error_file', metavar='FILE',
                    type=str,
                    help='Path to a spectrum file containing error data to be '
                    'used if a specific fits file was used as input')
parser.add_argument('--spectrum_dir', '-f', dest='spectrum_dir', metavar='DIR',
                    type=str,
                    help='Specify directory in which input spectrum files are '
                    'stored')
parser.add_argument('--template_dir', '-t', dest='template_dir', metavar='DIR',
                    type=str,
                    help='Specify directory in which input templates files '
                    'are stored')
parser.add_argument('--linecatalog', '-y', dest='linecatalog', metavar='FILE',
                    type=str,
                    help='Path to the rest lines catalog file')
parser.add_argument('--calibration_dir', '-k', dest='calibration_dir',
                    metavar='DIR', type=str,
                    help='Specify directory in which calibration files are '
                    'stored')
parser.add_argument('--zclassifier_dir', '-z', dest='zclassifier_dir',
                    metavar='DIR', type=str,
                    help='Specify directory in which zClassifier files are '
                    'stored')
parser.add_argument('--config', '-c', dest='config', metavar='FILE', type=str,
                    help='Json configuration file giving all these command '
                    'line parameters.')
parser.add_argument('--log_level', '-l', dest='log_level', metavar='LEVEL',
                    type=log_level,
                    help='Verbosity level. Either "none", "debug", "detail",'
                    '"info", "warning", "error" or "critical".')
parser.add_argument('--linecatalog_convert', '-a', action='store_true',
                    help='Convert the line catalog from Vacuum to Air')
parser.add_argument('--version', '-v', action='version', version=get_version(),
                    help='Print version and exit.')
