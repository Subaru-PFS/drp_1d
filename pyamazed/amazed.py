#!/usr/bin/env python
import os.path
from .argumentparser import parser
from .redshift import CSpectrum_default, CProcessFlowContext, \
    CProcessFlow, CLog, CParameterStore, CClassifierStore, \
    CLogFileHandler, CRayCatalog, CTemplateCatalog, get_version
from .catalog import FitsTemplateCatalog
from .config import Config
import json


zlog = None

def datapath(config, *path):
    return os.path.expanduser(os.path.join(config.data_path, *path))


def calibrationpath(config, *path):
    return os.path.expanduser(os.path.join(config.calibration_dir, *path))


def spectrumpath(config, *path):
    return os.path.expanduser(os.path.join(config.spectrum_dir, *path))


def update_paramstore(param, config):
    param.Set_String('calibrationDir',
                     os.path.expanduser(config.calibration_dir))


def process_spectrum(spectrum_path, noise_path, proc_id,
                     template_catalog, line_catalog,
                     param, classif):
    """Process a spectrum.

    Return a CProcessFlowContext object"""

    spectrum = CSpectrum_default()
    try:
        spectrum.LoadSpectrum(spectrum_path, noise_path)
    except Exception as e:
        zlog.LogWarning("Can't load spectrum : {}".format(e))

    ctx = CProcessFlowContext()
    ctx.Init(spectrum, proc_id, template_catalog, line_catalog, param, classif)

    pflow = CProcessFlow()
    pflow.Process(ctx)

    return ctx


class Logger(CLog):
    pass


def amazed():
    args = parser.parse_args()
    config = Config(args)

    zlog = CLog()
    logConsoleHandler = CLogFileHandler(zlog,
                                        os.path.join(config.output_folder,
                                                     'amazed.log'))
    logConsoleHandler.SetLevelMask(config.log_level)

    param = CParameterStore()
    param.Load(os.path.expanduser(config.parameters_file))
    update_paramstore(param, config)
    opt_saveIntermediateResults = param.Get_String('SaveIntermediateResults',
                                                   'all')

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

    retcode, opt_medianKernelWidth = param.Get_Float64("continuumRemoval.medianKernelWidth")
    assert retcode

    retcode, opt_nscales = param.Get_Float64("continuumRemoval.decompScales",
                                             8.0)
    assert retcode

    retcode, dfBinPath = param.Get_String("continuumRemoval.binPath",
                                          "absolute_path_to_df_binaries_here")
    assert retcode

    template_catalog = CTemplateCatalog(medianRemovalMethod,
                                        opt_medianKernelWidth,
                                        opt_nscales, dfBinPath)
    # template_catalog = FitsTemplateCatalog(medianRemovalMethod,
    #                                       opt_medianKernelWidth,
    #                                       opt_nscales, dfBinPath)
    zlog.LogInfo("Loading %s" % config.template_dir)

    try:
        template_catalog.Load(calibrationpath(config, config.template_dir))
    except Exception as e:
        zlog.LogCritical("Can't load template : {}".format(e))
        raise

    line_catalog = CRayCatalog()
    zlog.LogInfo("Loading %s" % config.linecatalog)
    line_catalog.Load(calibrationpath(config, config.linecatalog))
    if config.linecatalog_convert:
        line_catalog.ConvertVacuumToAir()

    for line in spectrumList:
        spectrum_path, noise_path, proc_id = line.split()
        try:
            ctx = process_spectrum(spectrumpath(config, spectrum_path),
                                   spectrumpath(config, noise_path),
                                   proc_id, template_catalog,
                                   line_catalog, param, classif)
        except Exception as e:
            zlog.LogError("Can't process spectrum "
                          "{} / {} : {}".format(spectrum_path,
                                                noise_path, e))
            continue

        try:
            ctx.GetDataStore().SaveRedshiftResult(config.output_folder)
        except Exception as e:
            zlog.LogWarning("Can't save redshift for "
                            "{} / {} : {}".format(spectrum_path,
                                                  noise_path, e))

        else:
            ctx.GetDataStore().SaveAllResults(os.path.join(config.output_folder,
                                                           proc_id),
                                              'all')
        # ctx.GetDataStore().SaveReliabilityResult('/tmp/bar')

    # save cpf-redshift version in output dir
    with open(os.path.join(config.output_folder, 'version.json'), 'w') as f:
        json.dump({'cpf-redshift-version': get_version()}, f)


if __name__ == '__main__':
    amazed()
