#!/usr/bin/env python
import os.path
from .argumentparser import parser
from .redshift import CSpectrum_default, CProcessFlowContext, \
    CProcessFlow, CLog, CParameterStore, CClassifierStore, \
    CLogFileHandler, CRayCatalog, CTemplateCatalog, get_version
from .catalog import FitsTemplateCatalog
from .config import Config
import json


def absolutepath(config, *path):
    return os.path.expanduser(os.path.join(*path))


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
    spectrum.LoadSpectrum(spectrum_path, noise_path)

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

    if os.path.exists(config.output_folder):
        message = "Output directory {} already exists.".format(config.output_folder)
        zlog.LogError(message)
        raise Exception(message)

    param = CParameterStore()
    param.Load(os.path.expanduser(config.parameters_file))
    update_paramstore(param, config)
    opt_saveIntermediateResults = param.Get_String('SaveIntermediateResults',
                                                   'all')

    classif = CClassifierStore()

    if config.zclassifier_dir:
        classif.Load(absolutepath(config, config.zclassifier_dir))

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
        template_catalog.Load(absolutepath(config, config.template_dir))
    except Exception as e:
        zlog.LogError("Can't load template : {}".format(e))
        raise

    line_catalog = CRayCatalog()
    zlog.LogInfo("Loading %s" % config.linecatalog)
    line_catalog.Load(absolutepath(config, config.linecatalog))
    if config.linecatalog_convert:
        line_catalog.ConvertVacuumToAir()

    for line in spectrumList:
        spectrum_path, noise_path, proc_id = line.split()
        target_path = os.path.join(config.output_folder, proc_id)
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

            ctx.GetDataStore().SaveAllResults()

        try:
            ctx.GetDataStore().SaveRedshiftResult(target_path)
            ctx.GetDataStore().SaveCandidatesResult(target_path);
            ctx.GetDataStore().SaveReliabilityResult(target_path);
            ctx.GetDataStore().SaveStellarResult(target_path);
            ctx.GetDataStore().SaveClassificationResult(target_path);
        except Exception as e:
            zlog.LogWarning("Can't save redshift for "
                            "{} / {} : {}".format(spectrum_path,
                                                  noise_path, e))

        if config.save_intermediate_results != 'no':
            ctx.GetDataStore().SaveAllResults(target_path,
                                              config.save_intermediate_results)

    # save cpf-redshift version in output dir
    with open(os.path.join(config.output_folder, 'version.json'), 'w') as f:
        json.dump({'cpf-redshift-version': get_version()}, f)


def main():
    try:
        amazed()
    except Exception as e:
        print("Critical error: {}".format(e))

if __name__ == '__main__':
    main()
