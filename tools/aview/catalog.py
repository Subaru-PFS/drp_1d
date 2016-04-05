# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 14:39:53 2015

@author: aschmitt
"""
import os
import re
from astropy.io import fits

import matplotlib.pyplot as pp
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from scipy.optimize import fsolve

import modelresult

class Catalog(object):
    def __init__(self, spath, ctype='undef'):
        self.logTagStr = "Catalog"
        self.spath = spath
        self.name = os.path.basename(spath)
        self.ctype = ctype
        self.n = -1
        self.linelambda = -1
        self.linename = -1
        self.lineforce = -1
        self.linetype = -1
        self.lineprofile = -1
        self.linegroup = -1
        self.linenominalamp = -1
        
        #self.conversion = "toAir"
        #self.conversion = "toVacuum"#"toAir"
        self.conversion = ""
        
        self.load()
        
    def load(self):
        filename = self.spath
        name = []
        linelambda = []
        force = []
        linetype = []
        
        lineprofile = []
        linegroup = []
        linenominalamp = []
        
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                #data = lineStr.split("\t| ")
                data = re.split("\t| ", lineStr)
                data = [r for r in data if r != '']
                #print len(data)
                if(len(data) >=4):
                    # fill the list                
                    #print data[0]
                    #print data[1]
                    #_name  = data[1].replace(',', '_')
                    _name  = data[1]
                    name.append(_name)
                    linelambda.append(float(data[0]))
                    force.append(data[3])
                    linetype.append(data[2])
                    if(len(data) >=5):
                        lineprofile.append(data[4])
                    else:
                        lineprofile.append("-1")
                    if(len(data) >=7):
                        linegroup.append(data[5])
                        linenominalamp.append(float(data[6]))
                    else:
                        linegroup.append("-1")
                        linenominalamp.append(-1) 
                    
        f.close()
        self.n = len(linelambda)
        #print('len wave = {0}'.format(self.n))
        #---- default xaxis index array
        self.linelambda = range(0,self.n)
        self.linename = range(0,self.n)
        self.lineforce = range(0,self.n)
        self.linetype = range(0,self.n)  
        self.lineprofile = range(0,self.n)  
        self.linegroup = range(0,self.n)  
        self.linenominalamp = range(0,self.n)  
        for x in range(0,self.n):
            if not self.conversion=="":
                if self.conversion=="toAir":
                    lval = self.vacuumToAir(linelambda[x])
                    self.ctype = "{}_{}".format(self.ctype, "convertedToAir")
                    print("WARNING: converting to AIR !!")
                elif self.conversion=="toVacuum":                
                    lval = self.airToVacuum(linelambda[x])
                    self.ctype = "{}_{}".format(self.ctype, "convertedToVacuum")
                    print("WARNING: converting to VACUUM !!")
            else:
                lval = linelambda[x]
            self.linelambda[x] = lval
            self.linename[x] = name[x]
            self.lineforce[x] = force[x]
            self.linetype[x] = linetype[x]
            self.lineprofile[x] = lineprofile[x]
            self.linegroup[x] = linegroup[x]
            self.linenominalamp[x] = linenominalamp[x]
     
    def save(self, outputPath):
        text_file = open(outputPath, "w")
        data = "#version:0.3.0\n"
        text_file.write("{}".format(data))
        data = "#lambda	Name	type	force	profile	group	nominal_ampl"
        text_file.write("{}\n".format(data))
                
        for k in range(self.n):
            data = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.linelambda[k], self.linename[k], self.linetype[k], self.lineforce[k], self.lineprofile[k], self.linegroup[k], self.linenominalamp[k])
            text_file.write("{}".format(data))
        
        text_file.close()            
        
    def vacuumToAir(self, vacuumVal):
        s = (1e-4)/vacuumVal;
        coeff = 1 + 8.34254*1e-5 + (2.406147*1e-2)/(130-s*s) + (1.5998*1e-4)/(38.9-s*s)

        airVal = vacuumVal/coeff
        return airVal

    def shiftedVacuumToAir(self, vacuumVal, targetAirVal):
        return self.vacuumToAir(vacuumVal)-targetAirVal
        
    def airToVacuum(self, airVal):
        vacuumVal = fsolve(self.shiftedVacuumToAir, airVal, airVal)
        return vacuumVal
                
        
    def __str__(self):
        a = "\nCatalog: {0}\n".format(self.name)
        a = a + ("    type = {0}\n".format(self.ctype))
        a = a + ("    n = {0}\n".format(self.n))
        a = a + ("\n")

        
        for x in range(0,self.n):
            a = a + "{:<20}\t{}\t{}\n".format(self.linename[x], self.linetype[x], self.linelambda[x])
        
        return a
        
    def plot(self): 
        shiftedctlg = c.getShiftedCatalog(0.0, "A", -1)
        #print(shiftedctlg)
        self.linesx = shiftedctlg['lambda']
        self.linesxrest = shiftedctlg['lambdarest']
        self.linesname = shiftedctlg['name']
        self.linesforce = shiftedctlg['force']    
        
        
        #do the plotting
        pp.figure(figsize=(10,8))
        for k in range(len(self.linesx)):
            #x = self.linesx[k]*(1+self.z)
            x = self.linesx[k]
            pp.plot((x, x), (-1000,1000) , 'r-', label=self.linesname[k] )
            #pp.text(x, self.ymax*0.75, '{0}'.format(self.linesname[k]))
            pp.text(x, 900, '{0}'.format(self.linesname[k]))
        
        pp.grid(True) # Affiche la grille
        #pp.legend(('spectrum','shifted template'), 'lower left', shadow = True)
        pp.xlabel('Angstrom')
        pp.ylabel('y')
        pp.title(self.name) # Titre
        #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
        pp.show()
    
    def plotInZplane(self):
        if 0:
            #EUCLID
            obs_lambda_min = 12500.0
            obs_lambda_max = 17500.0
        if 1:        
            #PFS
            obs_lambda_min = 3800.0
            obs_lambda_max = 12600.0
        if 0:        
            #PFS
            obs_lambda_min = 3600.0
            obs_lambda_max = 9400.0

        filter_type = "E"
        filter_force = -1
        ctlg_rest = self.getShiftedCatalog(0.0, filter_type, filter_force)       
        nlines = len(ctlg_rest['lambda']) 
        print("nlines={}".format(nlines))
        print("catalog rest = {}".format(ctlg_rest))
        linesaxisticks = ["{} - {}A".format(ctlg_rest['name'][a], ctlg_rest['lambdarest'][a]) for a in range(len(ctlg_rest['name']))]
        
        zmin = 0.0
        zmax = 10.0
        zstep = 0.01
        nz = int((zmax-zmin)/zstep)
        zaxis = np.linspace(zmin, zmax, nz)
        decimateZAxisN = int(nz/25.0)
        zaxisticks = ["{:.1f}".format(a) for a in zaxis[::decimateZAxisN]]
        zaxisticksinds = [int(a) for a in range(nz)[::decimateZAxisN]]
        print zaxisticksinds
                
        print("zaxis : n={}".format(len(zaxis)))
        matrix = np.zeros((nz, nlines))

        for iz,z in enumerate(zaxis):
            ctlg = self.getShiftedCatalog(z, filter_type, filter_force)
            print("Processing for z={}".format(z))
            #print("catalog = {}".format(ctlg))
            for ic in range(nlines):
                a = ctlg['lambda'][ic]
                #print("lambda for line={} and z={} is: lambda={}".format(ctlg['name'][ic],z,ctlg['lambda'][ic])) 
                if a >= obs_lambda_min and a<=obs_lambda_max:
                    matrix[iz,ic] = 1.0
                else:
                    matrix[iz,ic] = -0.50
        
        pp.clf()
        pp.close()
        #fig = pp.figure('catalog z map', figsize=(9, 8))
        #ax = fig.add_subplot(111)
        #i = ax.matshow(np.transpose(matrix), interpolation='nearest', aspect='equal', cmap=cmap)
        
        cmap = pp.get_cmap('Greys')
        pp.matshow(np.transpose(matrix), interpolation='nearest',aspect='auto', cmap=cmap)
        pp.xticks(zaxisticksinds, zaxisticks, rotation='horizontal',verticalalignment='bottom')
        pp.yticks(range(nlines), linesaxisticks, rotation='horizontal',verticalalignment='bottom')
        
        #pp.subplots_adjust(bottom=0.15)
        
        pp.xlabel('z')
        pp.ylabel('LINE')
        name1 = "Lines presence = f(z) for observed spectrum in [{:.1f}A - {:.1f}A]\nblack=present, white=absent".format(obs_lambda_min, obs_lambda_max)
        pp.title(name1)
        
        
        cmin = -0.5
        cmax = 1.75
        pp.clim(cmin,cmax)

        #ml = MultipleLocator(5)
        #pp.axes().yaxis.set_minor_locator(ml)

        pp.grid(True,'major',linewidth=1)
        #pp.grid(True,'minor', linewidth=1)
        pp.show()

#    def getFilteredCatalog(self, ltype='E', lforce='S'):
#        inds = []
#        for x in range(0,self.n):
#            if (self.linetype[x]==ltype or ltype==-1) and (self.lineforce[x]==lforce or lforce==-1):
#                inds.append(x)
#            
#        _linelambdarest = [x for i,x in enumerate(self.linelambda) if i in inds]
#        _linename = [x for i,x in enumerate(self.linename) if i in inds]
#        _lineforce = [x for i,x in enumerate(self.lineforce) if i in inds]
#        _linetype = [x for i,x in enumerate(self.linetype) if i in inds]
#        
#        return {"lambda":_linelambda, "lambdarest":_linelambdarest, "name":_linename,
#                "force":_lineforce, "type":_linetype} 
    
    def getShiftedCatalog(self, z, lineTypeFilter=-1, lineForceFilter=-1):
        coeff= 1.0+z
        linelambda = range(0,self.n)
        for x in range(0,self.n):
            linelambda[x] = self.linelambda[x]*coeff
            
        if lineTypeFilter == -1 and lineForceFilter == -1:
            inds = [i for i,x in enumerate(self.linetype)]
            #print("getShiftedCatalog, case 1")
        elif lineForceFilter == -1 and not lineTypeFilter == -1:
            inds = [i for i,x in enumerate(self.linetype) if x == lineTypeFilter]
            #print("getShiftedCatalog, case 2")
        elif lineTypeFilter == -1 and not lineForceFilter == -1:
            inds = [i for i,x in enumerate(self.lineforce) if x == lineForceFilter]
            #print("getShiftedCatalog, case 3")
        else:
            inds = [i for i,x in enumerate(self.lineforce) if x == lineForceFilter and self.linetype[i]==lineTypeFilter]
            #print("getShiftedCatalog, case 4")
            
        _linelambda = [x for i,x in enumerate(linelambda) if i in inds]
        _linelambdarest = [x for i,x in enumerate(self.linelambda) if i in inds]
        _linename = [x for i,x in enumerate(self.linename) if i in inds]
        _lineforce = [x for i,x in enumerate(self.lineforce) if i in inds]
        _linetype = [x for i,x in enumerate(self.linetype) if i in inds]
        
        return {"lambda":_linelambda, "lambdarest":_linelambdarest, "name":_linename,
                "force":_lineforce, "type":_linetype} 

    def getAmbiguityZ(self, zreference, zmin, zmax, lineTypeFilter=-1, lineForceFilterReference=-1, lineForceFilterFail=-1):
        print("getting ambiguity for type = {}, refForce = {}, failforce = {}".format(lineTypeFilter, lineForceFilterReference, lineForceFilterFail))
        ambiguityZs = []
        shiftedCatalogReference = self.getShiftedCatalog(0.0, lineTypeFilter, lineForceFilterReference)
        shiftedCatalogFail = self.getShiftedCatalog(0.0, lineTypeFilter, lineForceFilterFail)
        
        for x in range(0, len(shiftedCatalogReference["lambda"])):
            for y in range(0, len(shiftedCatalogFail["lambda"])):
                if shiftedCatalogReference["type"][x] == shiftedCatalogFail["type"][y]:
                    zfail = shiftedCatalogReference["lambda"][x]*(1+zreference)/shiftedCatalogFail["lambda"][y] -1
                    if not abs(zfail - zreference) < 0.0005 and zfail>zmin and zfail<zmax:
                        ambiguityZs.append({"z":zfail, "refname":shiftedCatalogReference["name"][x]\
                        , "reftype":shiftedCatalogReference["type"][x], "failname":shiftedCatalogFail["name"][y]\
                        , "failtype":shiftedCatalogFail["type"][y]})
                else:
                    #print("ambiguity type ref = {}, type fail ={}".format(shiftedCatalogReference["type"][x], shiftedCatalogFail["type"][y]))
                    pass
        return ambiguityZs
        
    def getRedmineTableString(self, ):
        outStr = ""      
        
        outStr = outStr + "|_.lambda|_.Name|_.type|_.force|_.profile|_.group|_.nominal_ampl|" + "\n"
        for x in range(0,self.n):
            nominalamp = ""
            if not self.linenominalamp[x] == -1:
               nominalamp = str(self.linenominalamp[x]) 
            group = ""
            if not self.linegroup[x] == "-1":
               group = self.linegroup[x]
            outStr = outStr + "|{}|{}|{}|{}|{}|{}|{}|".format(self.linelambda[x], self.linename[x], self.linetype[x], self.lineforce[x], self.lineprofile[x], group, nominalamp) + "\n"
      
        return outStr
    
            
    def getLatexTableString(self, ):
        outStr = ""      
        
        outStr = outStr + "\\begin{tabular}{ | c | c | c | c | c | c | c | } \n"
        outStr = outStr + "\\hline \n"
        outStr = outStr + "lambda & Name & type & force & profile & group & nominal_ampl \\\\ \n"
        outStr = outStr + "\\hline \n"
        
        for x in range(0,self.n):
            nominalamp = ""
            if not self.linenominalamp[x] == -1:
               nominalamp = str(self.linenominalamp[x]) 
            group = ""
            if not self.linegroup[x] == "-1":
               group = self.linegroup[x]
            outStr = outStr + "{}&{}&{}&{}&{}&{}&{} \\\\ \n".format(self.linelambda[x], self.linename[x], self.linetype[x], self.lineforce[x], self.lineprofile[x], group, nominalamp) + "\n"
        
        outStr = outStr + "\\hline \n"
        outStr = outStr + "\end{tabular}\n" 

        #correcting underscores for lateX
        outStr = outStr.replace("_", "\\_")      
        
        return outStr

    def applyTemplateShapeFromLinemodelFitResult(self, lmFitResFilePath):
        mres = modelresult.ModelResult(lmFitResFilePath)
        
        #
        for x in range(0,self.n):        
            a = mres.getAmplitude( self.linename[x], self.linetype[x])
            self.linegroup[x] = "tplshape"
            self.linenominalamp[x] = a
        #path = "/home/aschmitt/code/python/linemodel_tplshape/amazed/linecatalogs"
        #name = "linecatalogamazedvacuum_B9D_LyaAbsSYMprofile_tplShape.txt"
        #outpath = os.path.join(path,name)
        #self.save(outpath)
        
     
            
if __name__ == '__main__':
    path = "/home/aschmitt/gitlab/cpf-redshift/tools/simulation/amazed/linecatalogs"
    name = "linecatalogamazedvacuum_B10C.txt"
    cpath = os.path.join(path,name)
    print('using full path: {0}'.format(cpath))
    c = Catalog(cpath, ctype="vacuum")
    print(c) 
    #print(c.getShiftedCatalog(1.0, "E"))
    
    #c.plot()
    #c.plotInZplane()  
    
    #print("the REDMINE (copy/paste) generated table is:\n{}".format(c.getRedmineTableString()))
    print("the LATEX (copy/paste) generated table is:\n{}".format(c.getLatexTableString()))
    
    #lmResPath = "/home/aschmitt/code/python/linemodel_tplshape/amazed/output/spectrum_tpl_NEW_Im_extended.dat_TF/linemodelsolve.linemodel_fit_extrema_0.csv"
    #c.applyTemplateShapeFromLinemodelFitResult(lmResPath)