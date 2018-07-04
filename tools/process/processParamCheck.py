#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:24:40 2018

@author: aschmitt
"""

import os
import sys
import time
import argparse
import json


def findkeys(node, kv):
    if isinstance(node, list):
        for i in node:
            for x in findkeys(i, kv):
               yield x
    elif isinstance(node, dict):
        if kv in node:
            yield node[kv]
        for j in node.values():
            for x in findkeys(j, kv):
                yield x
                
class processParamCheck(object):
    def __init__(self, parampath):
        self.logTagStr = "processParamCheck"
        self.parampath = parampath
                    
    def checkAll(self):
        """
        return 0 if no errors detected
        """
        ret_errors = 0
        
        ret_errors += self.checkParamJsonForDeprecatedKeywords()
        ret_errors += self.checkParamJsonForBestUseInPractice()
        ret_errors += self.checkParamJsonForDeprecatedParamTypes()
        
        
        print("\n****************************************")
        if ret_errors==0:
            print(" JSON check RESULT: no errors found")
        else:
            print("ERROR: some errors were found")
        print("****************************************")
            
            
        return ret_errors
        
    def checkParamJsonForDeprecatedKeywords(self):
        """
        check the presence of deprecated keyword in the json
        return 1 if a deprecated is found
        return 0 if no deprecated is found
        """
        deprecated_keywords = ["lambdaRange", "redshiftRange", "redshiftStep"]
        print("Info: checking for the following deprecated keywords: \n{}\n".format(deprecated_keywords))
        
        f = open(self.parampath, 'r')
        for line in f:
            for w in deprecated_keywords:
                if w in line:
                    print("WARNING: found deprecated keyword in json param file ! ({} found in {})".format(w, line))
                    return 1
        
        return 0 
       

    
    def checkParamJsonForDeprecatedParamTypes(self, verbose=False):
        """
        check some parameters types. 
        In that test, it was a string initially, then has been modified to be a float value
        return 1 (error) if a string is found
        return 0 (success) if a float is found
        """
        target_keywords = ["stronglinesprior"]
        print("Info: checking for the following deprecated types: \n{}\n".format(target_keywords))
        
        f = open(self.parampath, 'r')
        data = json.load(f)
        for w in target_keywords:
            occurences = list(findkeys(data, w))
            if verbose:
                print("key {}: found n={} times".format(w, len(occurences)))
            
            for k, o in enumerate(occurences):
                if verbose:
                    print("INFO: occurence = {}".format(o))
                    print("INFO: type occurence = {}".format(type(o)))
                try:
                    float(o)
                except:
                    print("WARNING: found deprecated value type (STRING that should be FLOAT now) in json param file ! ({} found with string value = {})".format(w, o))
                    return 1
        return 0
        
        
    def checkParamJsonForBestUseInPractice(self):
        """
        This will inform if some keywords were found in the json. The keywords searched for here should probably better be defined in the config file or as cmd line args in the best classical use of amazed
        """
        somewhereelse_keywords = ["SaveIntermediateResults", "calibrationDir"]
        print("Info: checking for the following 'better-defined-somewhere-else' keywords: \n{}\n".format(somewhereelse_keywords))
        
        f = open(self.parampath, 'r')
        for line in f:
            for w in somewhereelse_keywords:
                if w in line:
                    print("WARNING: found best-defined-somewhere-else (in the config file ?) keyword in json param file ! ({} found in {})".format(w, line))
                    return 1
        
        return 0

        
        
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-p", "--paramPath", dest="paramPath", default="",
                    help="path to the param (JSON) file to be checked") 
                    
   
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.paramPath) :
        rpath = os.path.abspath(options.paramPath)
        print('INFO: using param full path: {}'.format(rpath))

                
        pcheck = processParamCheck(parampath=rpath)
        pcheck.checkAll()
            
    else :
        print("Error: invalid arguments")
        exit()
    
    
    
def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print(".")
    Main( sys.argv )
    
