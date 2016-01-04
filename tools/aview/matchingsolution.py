# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 14:42:06 2015

@author: aschmitt
"""


import os
import re

class SolutionSet(object):
    def __init__(self, matchN, avgZ):
        self.matchN = matchN
        self.avgZ = avgZ
    

class MatchingSolution(object):
    def __init__(self, spath, ctype='solution'):
        self.logTagStr = "MatchingSolution"
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
                if(len(data) >=2):
                    matchnum = int(data[0])
                    ndata = int(len(data))
                    avgZ = float(data[ndata-1])
                    solset = SolutionSet(matchnum, avgZ)
                    self.list.append(solset)
                    self.n += 1
                    #print("matchingsolution imported : n={}, z={}".format(solset.matchN, solset.avgZ))
        f.close()
        print("loaded {} match sol. results".format(self.n))

    def getMatchNumMax(self):
        maxi = 0
        
        for a in range(self.n):
            if maxi<self.list[a].matchN:
                maxi = self.list[a].matchN
        return maxi
                   
    def getNCorrectSolution(self, zref, thres = 0.0025):
        matchnummax = self.getMatchNumMax()
        nOK = 0
        nFAIL = 0
        for a in range(self.n):
            if True:#self.list[a].matchN > 1 :
                if abs(self.list[a].avgZ - zref) < thres:
                    nOK +=1 
                else:
                    nFAIL +=1
        print [nOK, nFAIL]
        return [nOK, nFAIL]
        
      
    def __str__(self):
        a = "\nMatchingSolution: {0}\n".format(self.name)
        a = a + ("    type = {0}\n".format(self.ctype))
        a = a + ("    n = {0}\n".format(self.n))
        a = a + ("\n")
        
        return a
        
      
      
            
if __name__ == '__main__':
    path = "/home/aschmitt/data/pfs/pfs_reallyjustline/amazed/res_20150826_linematching2_cut_analaysis/res_20150826_linematching2_cut1.5/EZ_fits-W-F_0"
    name = "linematching2solve.raymatching.csv"
    mspath = os.path.join(path,name)
    print('using full path: {0}'.format(mspath))
    c = MatchingSolution(mspath)
    print(c) 
    #print(c.getShiftedCatalog(1.0, "E"))
    #c.plot()