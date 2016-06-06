# -*- coding: utf-8 -*-
"""
"""
import time
import sys
import os
import shutil
import math
import random
import inspect

subfolder = "../aview"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],subfolder)))
if cmd_subfolder not in sys.path:
    #print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
    sys.path.insert(0, cmd_subfolder)
import spectrum
import modelresult

import utilbins

#parameters
enablePlot = False
optMission = "euclid"
#optMission = "pfs"

if optMission == "euclid":
    n_count_per_bin = 100
    lambda_obs_min = 12500.0
    lambda_obs_max = 18500.0
    dlambda = 13.5
    R=250.0;
    optISM = "fixed0.0"
    optLineAmplitudes = "file"
    #optLineAmplitudes = "pseudorandom"
    optVelocity = "fixed"
    optLineWidth = "nispsim2016"
    optRules = "no"
elif optMission == "pfs":
    n_count_per_bin = 10
    lambda_obs_min = 3800.0
    lambda_obs_max = 12600.0
    dlambda = 0.6
    R=2350.0;
    optISM = "fixed0.25"
    optLineAmplitudes = "pseudorandom"
    optVelocity = "random"
    optLineWidth = "combined"
    optRules = "balmer"


temp_directory = "./tmp/"
amazed_directory = "./amazed/"
output_directory = "./simulation_{}_output/".format(optMission)
refz_fname = os.path.join(output_directory, "refz.txt")
if os.path.exists(refz_fname):
    print("Deleting existing refz file: {}".format(refz_fname))
    os.remove(refz_fname)
if os.path.exists(output_directory):
    shutil.rmtree(output_directory) 
os.mkdir(output_directory) 

text_file = open(refz_fname, "w")
data = "#id   Z   MAG   TplID   E(B-V)   SFR_COEFF   ISM_COEFF  SigmaE   SigmaA\n"
text_file.write("{}".format(data))
text_file.close() 

print("\n_\nCreate simulation set by z-mag-sfr bins")
ubins = utilbins.UtilBins(n_count_per_bin, optMission)
WarningCheckBins = raw_input("\n\nWARNING: Bins have been created: Please CHECK!\nPress any key to continue...".format())
        

allfull = False
while not allfull:
    
    print("\n\nNew Realization:")
    
    #prepare z, mag and template random variables
    z, mag, sfr = ubins.getRandZMagSfr_fromOpenBin()
    if z==-1:
        allfull = True
        print("z-mag-sfr bins are now full ! leaving processing loop...")
        break
    print("RANDOM: z = {}, mag = {}, sfr = {}".format(z, mag, sfr))
    
    tpl_name = ubins.getRandTemplate()
    print("RANDOM: tpl name = {}".format(tpl_name))
  
    #prepare the template spc with the given mag
    tplPath = os.path.join(ubins.templates_path, tpl_name)
    tpl = spectrum.Spectrum(tplPath ,'template', snorm=False)
    tpl.extendWavelengthRangeRed(lambda_obs_max)
    tpl.extendWavelengthRangeBlue(lambda_obs_min)
    
    tpl.applyLyaExtinction(redshift=z)
    tpl.applyRedshift(z)    
    tpl.setMagIAB(mag)
    tpl.applyLambdaCrop(lambda_obs_min, lambda_obs_max)
    tpl.interpolate(dx=0.1) #high sampling for the synthesis process

    
    print(tpl)
    exported_fits_path = tpl.exportFits(temp_directory, "spectrum_simu_tmp", addNoise=False, exportNoiseSpectrum=False)
    exported_fits_name = os.path.split(exported_fits_path)[1]
    if enablePlot:
        tpl.plot()

    #prepare the spectrumlist file for the amazed pipeline            
    dest_spclist_fname = os.path.join(temp_directory, "tmp.spectrumlist".format())
    text_file = open(dest_spclist_fname, "w")
    text_file.write("{}".format(exported_fits_name))
    text_file.close()    

    #prepare the config file for the amazed pipeline          
    # (reshift range and wavelength range)           
    config_tpl_fname = os.path.join(amazed_directory, "config_template.txt".format())
    with open(config_tpl_fname, 'r') as myfile:
        data=myfile.read()
    data = data.replace("lambdarange=", "lambdarange={} {}".format(lambda_obs_min, lambda_obs_max))
    data = data.replace("redshiftrange=", "redshiftrange={:.4f} {:.4f} {}".format(z, z, 0.0001))
    data = data.replace("__tmp_path__", "{}".format(os.path.abspath(temp_directory)))
    data = data.replace("__amazed_path__", "{}".format(os.path.abspath(amazed_directory)))
    #print("config file generated content is :\n{}".format(data))
    config_fname = os.path.join(amazed_directory, "config.txt".format())
    text_file = open(config_fname, "w")
    text_file.write("{}".format(data))
    text_file.close()
    

    #SigmaA = 311.0
    #SigmaE = 111.0 
    if optVelocity == "random":
        SigmaE = ubins.getRandVelocity()
        SigmaA = ubins.getRandVelocity()        
    elif optVelocity == "fixed":
        SigmaA = 175.0
        SigmaE = 175.0
        
    
    #prepare the parameters file for the amazed pipeline
    params_tpl_fname = os.path.join(amazed_directory, "parameters/parameters_template.json".format())
    with open(params_tpl_fname, 'r') as myfile:
        data=myfile.read()
    data = data.replace("lambdaRangeMinxxx", "{}".format(lambda_obs_min))
    data = data.replace("lambdaRangeMaxxxx", "{}".format(lambda_obs_max))
    data = data.replace("redshiftRangeMinxxx", "{}".format(z))
    data = data.replace("redshiftRangeMaxxxx", "{}".format(z))
    data = data.replace("velocityemissionxxx", "{}".format(SigmaE))
    data = data.replace("velocityabsorptionxxx", "{}".format(SigmaA))
    data = data.replace("instresolutionionxxx", "{}".format(R))
    data = data.replace("linewidthtypexxx", "{}".format(optLineWidth))
    data = data.replace("rulesxxx", "{}".format(optRules))
    
    #print("parameters file generated content is :\n{}".format(data))
    params_fname = os.path.join(amazed_directory, "parameters/parameters.json".format())
    text_file = open(params_fname, "w")
    text_file.write("{}".format(data))
    text_file.close()   

    #prepare the linemodel amplitude fit file for the linemodel
    coeffE_SFR = sfr/100.0
    if optISM == "fixed0.25":
        coeffA_ISM = 0.25
    elif optISM == "fixed0.0":
        coeffA_ISM = 0.0
    if optISM == "random":
        coeffA_ISM = 0.5*random.random()
    ism = coeffA_ISM
    print("USING: ism = {}".format(ism))
        
    mAmpsOutputpath = os.path.join(amazed_directory, "linecatalogs/linemodelsolve.linemodel_fit_extrema_0.csv".format())
    if optLineAmplitudes == "pseudorandom":        
        #open linemodel template for catalog linecatalogamazedvacuum_B9D_LyaAbsSYMprofile
        mAmpsTplpath = os.path.join(amazed_directory, "linecatalogs/linemodelsolve.linemodel_fit_extrema_0_template.csv".format())
        model_amps = modelresult.ModelResult(mAmpsTplpath)
        model_amps.randomAmplitudes(coeffE=coeffE_SFR, coeffA=coeffA_ISM)
    elif optLineAmplitudes == "file":
        tplmodelresultfilepath = "/home/aschmitt/gitlab/cpf-redshift/tools/simulation/templates_lines/201605022133_id0_z1.12688_mag20.7_sfr10.0_linemodel_amps_2thirdAmps.csv"   
        model_amps = modelresult.ModelResult(tplmodelresultfilepath)     
    model_amps.save(mAmpsOutputpath)
    
    #run the linemodel with the random line amplitude fitting option, with rules
    #current_path = os.getcwd()
    #os.chdir('./amazed')
    bashCommand = "~/gitlab/amazed/bin/amazed-0.0.0 -c {}".format(config_fname)
    os.system(bashCommand) 
    #os.chdir(current_path)

    print("Amazed pipeline finished...")            
     
    #open linemodel generated
    mresultpath = os.path.join(amazed_directory, "output/{}/linemodelsolve.linemodel_fit_extrema_0.csv".format(os.path.splitext(exported_fits_name)[0]))
    model_result = modelresult.ModelResult(mresultpath)
    mean_amp_emission_lines = model_result.getAmplitudeMeanNhighest()

    #sfr = mean_amp_emission_lines/coeff_mag_from_pfs  
    print("Amazed generated sfr = {}".format(sfr))
      
    
    ret = ubins.add(z, mag, sfr)
    if not ret ==1:
        print("unable to add the current spectrum... skip")
    else:
                
        #save simulation spectrum and noise
        modelPath = os.path.join(amazed_directory, "output/{}/linemodelsolve.linemodel_spc_extrema_0.csv".format(os.path.splitext(exported_fits_name)[0]))
        model = spectrum.Spectrum(modelPath ,'template', snorm=False)
        model.correctZeros()
        if enablePlot:
            model.plot()
        spc_simu_name = ubins.getSpcName(z, mag, sfr)
        
        if optMission=='pfs':
            model.interpolate(dx=dlambda)
            model.exportFits(output_directory, "{}".format(spc_simu_name), addNoise=False, exportNoiseSpectrum=False)
            model.exportFits(output_directory, "{}".format(spc_simu_name), addNoise=True, exportNoiseSpectrum=True)
        elif optMission=='euclid':
            _spc_name = spc_simu_name + "_raw"
            model.exportFits(output_directory, "{}".format(_spc_name), addNoise=False, exportNoiseSpectrum=False)
            model.exportFits(output_directory, "{}".format(_spc_name), addNoise=True, exportNoiseSpectrum=True)
            
            delta_rolls = [0.0, 0.25, 0.7, 0.5]
            for i,d in enumerate(delta_rolls):
                _spc_name = spc_simu_name + "_roll{}".format(i)
                model_roll = model.copy()
                model_roll.interpolate(dx=dlambda, offsetx=dlambda*d)
                model_roll.exportFits(output_directory, "{}".format(_spc_name), addNoise=False, exportNoiseSpectrum=False)
                model_roll.exportFits(output_directory, "{}".format(_spc_name), addNoise=True, exportNoiseSpectrum=True)
            
        #save the linemodel input amplitude fits csv file  
        model_amps_savePath = os.path.join(output_directory, "{}_linemodel_amps.csv".format(spc_simu_name))        
        model_amps.save(model_amps_savePath)
          
        #save simulation refz line 
        data = "{}\t{}\t{:.1f}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.1f}\t{:.1f}\n".format(spc_simu_name, z, mag, tpl_name, -1, sfr, ism, SigmaE, SigmaA)
        text_file = open(refz_fname, "a")
        text_file.write("{}".format(data))
        text_file.close()
        
        p = ubins.getPercentageFull()
        print("Added SPECTRUM: matrix is now : {:.3f} % full".format(p))
        if p==100.0:
            allfull = True
            print("z-mag-sfr bins are 100.0 full ! leaving processing loop...")
            
        #export bin count reports
        print("\nINFO: export bin count reports...")
        bin_count_path = os.path.join(output_directory, "report_bin_count.txt".format())    
        ubins.exportBinCount(bin_count_path)
        
        #export bin count plot : warning, this takes a few seconds to process... deactivate if time is an issue.
        if 0: 
            bin_count_path = os.path.join(output_directory, "plot_bin_count".format()) 
            if not os.path.exists(bin_count_path):
                os.mkdir(bin_count_path) 
            ubins.exportBinCountPlot(bin_count_path)


    #break
    #time.sleep(2)
    
    
print("END of createSimu.py script.")
