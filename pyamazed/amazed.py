#!/usr/bin/env python
import os.path
from argumentparser import parser
from redshift import *
from config import Config


def datapath(config, *path):
    return os.path.expanduser(os.path.join(config.data_path, *path))

def calibrationpath(config, *path):
    return os.path.expanduser(os.path.join(config.calibration_dir, *path))

def spectrumpath(config, *path):
    return os.path.expanduser(os.path.join(config.spectrum_dir, *path))

def update_paramstore(param, config):
    param.Set('calibrationDir', config.calibration_dir)

def amazed():

    args = parser.parse_args()
    config = Config(args)

    zlog = CLog()
    logConsoleHandler = CLogConsoleHandler( zlog )
    logConsoleHandler.SetLevelMask ( zlog.nLevel_Info )

    param = CParameterStore()
    param.Load(os.path.expanduser(config.parameters_file))
    update_paramstore(param, config)

    classif = CClassifierStore()

    if config.zclassifier_dir:
        classif.Load(calibrationpath(config, config.zclassifier_dir))

    spectrumList = []
    with open(os.path.expanduser(config.input_file), 'r') as f:
        for spectrum in f:
            if not spectrum.startswith('#'):
                spectrumList.append(spectrum.strip())

    retcode, medianRemovalMethod = param.Get_String("continuumRemoval.method",
                                                    "IrregularSamplingMedian")
    assert retcode

    retcode, opt_medianKernelWidth = param.Get_Float64("continuumRemoval.medianKernelWidth", 75)
    assert retcode

    retcode, opt_nscales = param.Get_Float64("continuumRemoval.decompScales", 8)
    assert retcode

    retcode, dfBinPath = param.Get_String("continuumRemoval.binPath",
                                          "absolute_path_to_df_binaries_here")
    assert retcode

    template_catalog = CTemplateCatalog(medianRemovalMethod, opt_medianKernelWidth,
                                        opt_nscales, dfBinPath)
    print("Loading %s" % config.template_dir)

    try:
        template_catalog.Load(calibrationpath(config, config.template_dir))
    except Exception as e:
        print(f"Can't load template : {e}")
        raise

    line_catalog = CRayCatalog()
    print("Loading %s" % config.linecatalog)
    line_catalog.Load(calibrationpath(config, config.linecatalog))

    print(spectrumList)

    for line in spectrumList:
        spectrum_path, noise_path, proc_id = line.split()
        spectrum = CSpectrum()
        try:
            spectrum.LoadSpectrum(spectrumpath(config, spectrum_path),
                                  spectrumpath(config, noise_path))
        except Exception as e:
            print("Can't load spectrum : {}".format(e))
            continue

        try:
            ctx = CProcessFlowContext()
            ctx.Init(spectrum,
                     proc_id,
                     template_catalog,
                     line_catalog,
                     param,
                     classif)
        except Exception as e:
            print("Can't init process flow : {}".format(e))
            continue

        pflow=CProcessFlow()
        try:
            pflow.Process(ctx)
        except Exception as e:
            print("Can't process : {}".format(e))
            continue

        ctx.GetDataStore().SaveRedshiftResult(config.output_folder)
        #ctx.GetDataStore().SaveReliabilityResult('/tmp/bar')
        ctx.GetDataStore().SaveAllResults(config.output_folder, 'all')

if __name__ == '__main__':
    amazed()
