# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 14:39:53 2015

@author: aschmitt
"""
import os
import re
from astropy.io import fits

import matplotlib.pyplot as pp
import numpy as np
from scipy.optimize import fsolve

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
        
        self.conversion = "" #"toAir" #"toVacuum"#"toAir"
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
                    lineprofile.append(data[4])
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
      
            
if __name__ == '__main__':
    path = "/home/aschmitt/data/muse/muse1_20160126/amazed/linecatalogs"
    name = "linecatalogamazedvacuum_B8B.txt"
    cpath = os.path.join(path,name)
    print('using full path: {0}'.format(cpath))
    c = Catalog(cpath, ctype="vacuum")
    print(c) 
    #print(c.getShiftedCatalog(1.0, "E"))
    c.plot()
    
    print("the REDMINE (copy/paste) generated table is:\n{}".format(c.getRedmineTableString()))