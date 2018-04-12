import os.path

from redshift import *
from argumentparser import parser

class TestReader(CSpectrumIOGenericReader):

    def Read(self, path, spectrum):
        print('there')
        print('reading {} into {}'.format(path,spectrum))
        res = super().Read(path, spectrum)
        print('here')
        return res


def amazed():

    args = parser.parse_args()
    print(args)
    zlog = CLog()
    logConsoleHandler = CLogConsoleHandler( zlog )
    logConsoleHandler.SetLevelMask ( zlog.nLevel_Info )

    param = CParameterStore()
    param.Load(os.path.join(args.data_path, args.parameters_file))

    classif = CClassifierStore()

    if args.zclassifier_dir:
        classif.Load(args.zclassifier_dir)

    spectrumList = []
    with open(os.path.join(args.data_path, args.input_file), 'r') as f:
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
    print("Loading %s" % args.template_dir)
    template_catalog.Load(args.template_dir)

    line_catalog = CRayCatalog()
    print("Loading %s" % args.linecatalog)
    line_catalog.Load(args.linecatalog)

    reader = TestReader()

    print(spectrumList)

    for line in spectrumList:
        spectrum, noise, proc_id = line.split()
        ctx=CProcessFlowContext()
        ctx.Init(os.path.join(args.spectrum_dir, spectrum),
                 os.path.join(args.spectrum_dir, noise),
                 reader,
                 reader,
                 proc_id,
                 template_catalog,
                 line_catalog,
                 param,
                 classif)

        pflow=CProcessFlow()
        pflow.Process(ctx)
        ctx.GetDataStore().SaveRedshiftResult(args.output_folder)
        #ctx.GetDataStore().SaveReliabilityResult('/tmp/bar')
        ctx.GetDataStore().SaveAllResults(args.output_folder, 'all')

if __name__ == '__main__':
    amazed()
