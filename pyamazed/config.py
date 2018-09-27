import json
import os.path
from .redshift import CLog

defaults = {
    'parameters_file': 'parameters.json',
    'config': None,
    'data_path': '',
    'output_folder': './output',
    'input_file': 'input.spectrumlist',
    'error_file': None,
    'spectrum_dir': 'spectrum',
    'template_dir': 'calibration/templates/ExtendedTemplatesJan2017_v3',
    'linecatalog':'linecatalog.txt',
    'calibration_dir': './calibration',
    'zclassifier_dir': '',
    'log_level': CLog.nLevel_Warning
    }

class Config(object):

    def __init__(self, args):

        # set up defaults
        for k, v in defaults.items():
            setattr(self, k, v)

        # load config file
        if args.config is not None:
            with open(os.path.expanduser(args.config), 'r') as f:
                cfg = json.load(f)
            for k, v in cfg.items():
                if k not in defaults:
                    raise AttributeError('Invalid config file parameter {}'.format(k))
                setattr(self, k, v)


        # override with command line parameters
        for arg in vars(args):
            if arg not in defaults:
                raise AttributeError('Invalid command line parameter {}'.format(arg))
            if getattr(args, arg) is not None:
                setattr(self, arg, getattr(args, arg))
