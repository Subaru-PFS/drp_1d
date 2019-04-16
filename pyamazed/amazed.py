#!/usr/bin/env python
import os.path
from collections import namedtuple
from .argumentparser import parser
from .redshift import CSpectrum_default, CProcessFlowContext, \
    CProcessFlow, CLog, CParameterStore, CClassifierStore, \
    CLogFileHandler, CRayCatalog, CTemplateCatalog, get_version
from .catalog import FitsTemplateCatalog
from .config import Config
from .workers import Local, PBS
import math
import json
import traceback
from tempfile import NamedTemporaryFile

Observation = namedtuple('Observation', ['spectrum', 'noise', 'id_'])
RunSetup = namedtuple('RunSetup', ['logger', 'template_catalog',
                                   'line_catalog', 'param', 'classif'])


def absolutepath(config, *path):
    return os.path.expanduser(os.path.join(*path))


def spectrumpath(config, *path):
    return os.path.expanduser(os.path.join(config.spectrum_dir, *path))


def update_paramstore(param, config):
    param.Set_String('calibrationDir',
                     os.path.expanduser(config.calibration_dir))


def _output_observations(observations):
    return ''.join(['{} {} {}\n'.format(obs.spectrum, obs.noise, obs.id_)
                    for obs in observations])


def _load_observations(path):
    observations = []
    with open(path, 'r') as f:
        for spectrum in f:
            if not spectrum.startswith('#'):
                observations.append(Observation(*(spectrum.strip().split())))
    return observations


def process_spectrum(config, observation,
                     template_catalog, line_catalog,
                     param, classif):
    """Process a spectrum.

    Return a CProcessFlowContext object"""

    target_path = os.path.join(config.output_dir, observation.id_)

    spectrum = CSpectrum_default()
    spectrum.LoadSpectrum(spectrumpath(config, observation.spectrum),
                          spectrumpath(config, observation.noise))

    ctx = CProcessFlowContext()
    ctx.Init(spectrum, observation.id_,
             template_catalog, line_catalog, param, classif)

    pflow = CProcessFlow()
    pflow.Process(ctx)

    try:
        ctx.GetDataStore().SaveRedshiftResult(config.output_dir)
        ctx.GetDataStore().SaveCandidatesResult(config.output_dir)
        ctx.GetDataStore().SaveReliabilityResult(config.output_dir)
        ctx.GetDataStore().SaveStellarResult(config.output_dir)
        ctx.GetDataStore().SaveQsoResult(config.output_dir)
        ctx.GetDataStore().SaveClassificationResult(config.output_dir)
    except Exception as e:
        raise Exception("Can't save redshift for "
                        "{} / {} : {}".format(observation.spectrum,
                                              observation.noise, e))

    if config.save_intermediate_results != 'no':
        ctx.GetDataStore().SaveAllResults(target_path,
                                          config.save_intermediate_results)


def _setup_pass(config):
    """Setup template catalog, line catalog, parameterstore and zclassifier"""

    zlog = CLog()
    logConsoleHandler = CLogFileHandler(zlog,
                                        os.path.join(config.output_dir,
                                                     'log.txt'))
    logConsoleHandler.SetLevelMask(config.log_level)

    param = CParameterStore()
    param.Load(os.path.expanduser(config.parameters_file))
    update_paramstore(param, config)
    saveIntermediateResults = param.Get_String('SaveIntermediateResults',
                                               'all')

    classif = CClassifierStore()

    if config.zclassifier_dir:
        classif.Load(absolutepath(config, config.zclassifier_dir))

    param.Set_String('linemeascatalog', config.linemeascatalog)

    medianRemovalMethod = param.Get_String("continuumRemoval.method",
                                           "IrregularSamplingMedian")
    medianKernelWidth = param.Get_Float64("continuumRemoval.medianKernelWidth")
    nscales = param.Get_Float64("continuumRemoval.decompScales", 8.0)
    dfBinPath = param.Get_String("continuumRemoval.binPath",
                                 "absolute_path_to_df_binaries_here")

    _template_dir = absolutepath(config, config.template_dir)
    if os.path.isfile(config.template_dir):
        # template_dir is actually a file: load templates from FITS
        template_catalog = FitsTemplateCatalog(medianRemovalMethod,
                                               medianKernelWidth,
                                               nscales, dfBinPath)
    elif os.path.isdir(config.template_dir):
        # template_dir a directory: load templates from calibration dirs
        template_catalog = CTemplateCatalog(medianRemovalMethod,
                                            medianKernelWidth,
                                            nscales, dfBinPath)
    zlog.LogInfo("Loading %s" % _template_dir)

    try:
        template_catalog.Load(_template_dir)
    except Exception as e:
        zlog.LogError("Can't load template : {}".format(e))
        raise

    line_catalog = CRayCatalog()
    zlog.LogInfo("Loading %s" % config.linecatalog)
    line_catalog.Load(absolutepath(config, config.linecatalog))
    if config.linecatalog_convert:
        line_catalog.ConvertVacuumToAir()

    return RunSetup(zlog, template_catalog, line_catalog, param, classif)


def process_spectra(config, observations):
    """Process a bunch of spectra using worker"""

    os.makedirs(config.output_dir)

    setup = _setup_pass(config)

    setup.logger.LogInfo('processing {}'.format(observations))
    for observation in observations:
        try:
            process_spectrum(config, observation, setup.template_catalog,
                             setup.line_catalog, setup.param, setup.classif)
        except Exception as e:
            print("Can't process spectrum "
                  "{} / {} : {}".format(observation.spectrum,
                                        observation.noise, e))

    # save cpf-redshift version in output dir
    with open(os.path.join(config.output_dir, 'version.json'), 'w') as f:
        json.dump({'cpf-redshift-version': get_version()}, f)

    # copy input files to output dir
    with open(os.path.join(config.output_dir,
                           'input.spectrumlist'), 'w') as f:
        f.write(_output_observations(observations))

    setup.param.Save(os.path.join(config.output_dir, 'parameters.json'))

    config.save(os.path.join(config.output_dir, 'config.json'))
    setup.line_catalog.Save(os.path.join(config.output_dir,
                                         'linecatalog.txt'))


def bunch(bunch_size, iterable):
    """Split iterable in bunches of "bunch_size" lists."""
    if bunch_size <= 0:
        return iterable

    start = 0
    for _ in range(int(math.ceil(len(iterable) / bunch_size))):
        yield iterable[start:start+bunch_size]
        start += bunch_size


def create_worker(worker, config):
    if worker.lower() == 'pbs':
        return PBS(config)
    else:
        return Local(config)


def amazed(config):
    """Split spectra into bunch and process them"""

    observations = _load_observations(os.path.expanduser(config.input_file))

    worker = create_worker(config.worker, config)

    if os.path.exists(config.output_dir):
        raise Exception("Output directory {} already exists.".format(
            config.output_dir))
    os.makedirs(config.output_dir)

    config_json = NamedTemporaryFile(dir=config.output_dir,
                                     prefix='config_', delete=False)
    config.save(config_json.name)

    observation_files = []

    for i, spectra in enumerate(bunch(int(config.bunch_size), observations)):
        spectralist = NamedTemporaryFile(mode='w', dir=config.output_dir,
                                         prefix='spectralist_', delete=False)
        spectralist.write(_output_observations(spectra))
        spectralist.flush()

        # keep a reference on temp file to avoid premature deletion
        observation_files.append(spectralist)

        args = ['process_spectra_worker',
                '--config={}'.format(config_json.name),
                '--output={}'.format(os.path.join(config.output_dir,
                                                  '{:04}'.format(i))),
                '--input={}'.format(spectralist.name)]
        worker.run(args)

    worker.wait_all()


def process_spectra_worker():
    args = parser.parse_args()
    config = Config(args)

    print('loading {}'.format(os.path.expanduser(config.input_file)))
    observations = _load_observations(os.path.expanduser(config.input_file))
    try:
        process_spectra(config, observations)
    except Exception as e:
        traceback.print_exc()
        print("Critical error: {}".format(e))


def main():
    try:
        args = parser.parse_args()
        config = Config(args)
        amazed(config)
    except Exception as e:
        traceback.print_exc()
        print("Critical error: {}".format(e))


if __name__ == '__main__':
    main()
