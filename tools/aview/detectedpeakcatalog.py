# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:31:43 2015

@author: aschmitt
"""

import os
import re

class Peak(object):
    def __init__(self, name, inf, sup):
        self.inf = inf
        self.sup = sup
        self.name = name

class DetectedPeakCatalog(object):
    def __init__(self, spath, ctype='peak'):
        self.logTagStr = "PeakCatalog"
        self.spath = spath
        self.name = os.path.basename(spath)
        self.ctype = ctype
        
        self.n = -1
        self.list = []
        
        self.load()
        
    def load(self):
        filename = self.spath
       
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                #data = lineStr.split("\t| ")
                data = re.split("\t| ", lineStr)
                data = [r for r in data if r != '']
                #print len(data)
                if(len(data) >=3):
                    # fill the list                
                    #print data[0]
                    #print data[1]
                    #_name  = data[1].replace(',', '_')
                    _name  = data[0]
                    _inf = float(data[1])
                    _sup = float(data[2])
                    self.list.append(Peak(_name, _inf, _sup))
                    self.n += 1
        f.close()
              
        
    def __str__(self):
        a = "\nDetectedCatalog: {0}\n".format(self.name)
        a = a + ("    type = {0}\n".format(self.ctype))
        a = a + ("    n = {0}\n".format(self.n))
        a = a + ("\n")
        
        return a
        
      
      
            
if __name__ == '__main__':
    path = "/home/aschmitt/data/pfs/pfs_reallyjustline/amazed/res_20150831_tf_linematching2_tol_analysis/res_20150831_tf_linematching2_tol_0.0001_withRules/EZ_fits-W-TF_13"
    name = "linematching2solve.peakdetection.csv"
    cpath = os.path.join(path,name)
    print('using full path: {0}'.format(cpath))
    c = DetectedPeakCatalog(cpath)
    print(c) 
    #print(c.getShiftedCatalog(1.0, "E"))
    #c.plot()