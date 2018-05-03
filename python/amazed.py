import os.path

from redshift import *
from argumentparser import parser
from config import Config

def fullpath(config, *path):
    return os.path.expanduser(os.path.join(config.data_path, *path))

class TestReader(CSpectrumIOGenericReader):

    def Read(self, path, spectrum):
        print('reading {} into {}'.format(path,spectrum))
        res = super().Read(path, spectrum)
        print('here')
        return res

def amazed():

    args = parser.parse_args()
    config = Config(args)

    zlog = CLog()
    logConsoleHandler = CLogConsoleHandler( zlog )
    logConsoleHandler.SetLevelMask ( zlog.nLevel_Info )

    param = CParameterStore()
    param.Load(fullpath(config, config.parameters_file))

    classif = CClassifierStore()

    if config.zclassifier_dir:
        classif.Load(fullpath(config, config.zclassifier_dir))

    spectrumList = []
    with open(fullpath(config, config.input_file), 'r') as f:
        for spectrum in f:
            if not spectrum.startswith('#'):
                spectrumList.append(fullpath(config, config.spectrum_dir, spectrum.strip()))

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
    template_catalog.Load(fullpath(config, config.template_dir))

    line_catalog = CRayCatalog()
    print("Loading %s" % config.linecatalog)
    line_catalog.Load(fullpath(config, config.linecatalog))

    reader = TestReader()

    print(spectrumList)

    for line in spectrumList:
        spectrum, noise, proc_id = line.split()
        ctx=CProcessFlowContext()
        ctx.Init(fullpath(config, config.spectrum_dir, spectrum),
                 fullpath(config, config.spectrum_dir, noise),
                 reader,
                 reader,
                 proc_id,
                 template_catalog,
                 line_catalog,
                 param,
                 classif)

        pflow=CProcessFlow()
        pflow.Process(ctx)
        ctx.GetDataStore().SaveRedshiftResult(config.output_folder)
        #ctx.GetDataStore().SaveReliabilityResult('/tmp/bar')
        ctx.GetDataStore().SaveAllResults(config.output_folder, 'all')

if __name__ == '__main__':
    amazed()
