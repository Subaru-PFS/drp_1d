# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 22:17:42 2015

@author: aschmitt
"""


import os
import re


class Catalog(object):
    def __init__(self, spath, ctype='detected'):
        self.logTagStr = "Catalog"
        self.spath = spath
        self.name = os.path.basename(spath)
        self.ctype = ctype
        self.n = -1
        
        self.linelambda = -1
        self.linename = -1
        self.lineforce = -1
        self.linecut = -1
        self.linetype = -1
        self.linewidth = -1
        
        self.load()
        
    def load(self):
        filename = self.spath
        name = []
        linelambda = []
        force = []
        cut = []
        width = []
        linetype = []
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                #data = lineStr.split("\t| ")
                data = re.split("\t| ", lineStr)
                data = [r for r in data if r != '']
                #print len(data)
                if(len(data) >=5):
                    # fill the list                
                    #print data[0]
                    #print data[1]
                    #_name  = data[1].replace(',', '_')
                    _name  = data[0]
                    name.append(_name)
                    linelambda.append(float(data[1]))
                    force.append(data[2])
                    cut.append(float(data[3]))
                    width.append(float(data[4]))
                    linetype.append('detected')
                    
        f.close()
        self.n = len(linelambda)
        #print('len wave = {0}'.format(self.n))
        #---- default xaxis index array
        self.linelambda = range(0,self.n)
        self.linename = range(0,self.n)
        self.linecut = range(0,self.n)
        self.lineforce = range(0,self.n)
        self.linetype = range(0,self.n) 
        self.linewidth =  range(0,self.n)
        for x in range(0,self.n):
            self.linelambda[x] =linelambda[x]
            self.linename[x] = name[x]
            self.lineforce[x] = force[x]
            self.linecut[x] = cut[x]
            self.linewidth[x] = width[x]
            self.linetype[x] = linetype[x]
                
        
    def __str__(self):
        a = "\nDetectedCatalog: {0}\n".format(self.name)
        a = a + ("    type = {0}\n".format(self.ctype))
        a = a + ("    n = {0}\n".format(self.n))
        a = a + ("\n")
        
        return a
        
      
      
            
if __name__ == '__main__':
    path = "/home/aschmitt/data/pfs/pfs_reallyjustline/amazed/res_20150826_linematching2_cut1/EZ_fits-W-F_0"
    name = "linematching2solve.raycatalog.csv"
    cpath = os.path.join(path,name)
    print('using full path: {0}'.format(cpath))
    c = Catalog(cpath)
    print(c) 
    #print(c.getShiftedCatalog(1.0, "E"))
    #c.plot()