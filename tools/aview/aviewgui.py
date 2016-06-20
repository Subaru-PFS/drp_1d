# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 11:21:31 2015

@author: aschmitt
"""
 
from PyQt4 import QtGui, QtCore
import sys
import os
import inspect
import glob

import aview
import resultstat
import resparser

subfolder = "../stats"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],subfolder)))
if cmd_subfolder not in sys.path:
    #print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
    sys.path.insert(0, cmd_subfolder)
     
import processAmazedOutputStats
 
class ObjectAViewGui(object):
    def __init__(self):
        self.value = []
 
 
class AViewGui(QtGui.QMainWindow):
    def __init__(self, obj):
        super(AViewGui, self).__init__()
        self.obj = obj
        wdg = QtGui.QWidget()
        layout = QtGui.QGridLayout(wdg)     
        
        #Add the result dir. setup ctrls
        self.lblResdir = QtGui.QLabel(' AMAZED Output Result Dir. ', wdg)
        layout.addWidget(self.lblResdir, 0, 0, 1, 1)
        
        self.leResDir = QtGui.QLineEdit(wdg)
        self.leResDir.setFixedWidth(500) 
        layout.addWidget(self.leResDir, 0, 1, 1, 10)
        
        self.btnBrowseResdir = QtGui.QPushButton(' Browse ', wdg)
        self.btnBrowseResdir.setToolTip('Browse to select the AMAZED output directory...')
        self.btnBrowseResdir.clicked.connect(self.bt_setResultDir)
        layout.addWidget(self.btnBrowseResdir, 0, 2, 1, 1)
        
        self.btnBrowseResdirPrevious = QtGui.QPushButton(' < ', wdg)
        self.btnBrowseResdirPrevious.setToolTip('Back to the previous ResDir...')
        self.btnBrowseResdirPrevious.clicked.connect(self.bt_setResultDirPrevious)
        layout.addWidget(self.btnBrowseResdirPrevious, 0, 3, 1, 1)
        
        #Add the failure threshold filter ctrls
        self.lblDiffThres = QtGui.QLabel(' Filter by z diff. threshold: ', wdg)
        layout.addWidget(self.lblDiffThres, 1, 0, 1, 1)
        
        self.leDiffThres = QtGui.QLineEdit(wdg)
        self.leDiffThres.setToolTip('Enter threshold value to filter the dataset by the z error value...')
        self.leDiffThres.setFixedWidth(200) 
        layout.addWidget(self.leDiffThres, 1, 1, 1, 10)
        
        #Add the spectrum filter ctrls
        self.lblSpcFilter = QtGui.QLabel(' Filter by spectrum name: ', wdg)
        layout.addWidget(self.lblSpcFilter, 2, 0, 1, 1)
        
        self.leSpcFilter = QtGui.QLineEdit(wdg)
        self.leSpcFilter.setToolTip('Enter a tag to filter the dataset by the spectrum name...')
        self.leSpcFilter.setFixedWidth(300) 
        layout.addWidget(self.leSpcFilter, 2, 1, 1, 10)
        
        #Add the show button
        self.btnLoad = QtGui.QPushButton('Load Result List', wdg)
        self.btnLoad.setFixedWidth(500)
        self.btnLoad.setFixedHeight(50)
        self.btnLoad.clicked.connect(self.bt_loadresults)
        layout.addWidget(self.btnLoad, 4, 1, 1, 1)

        #Add the result list section separator
        self.lblResultsSection = QtGui.QLabel('----------', wdg)
        layout.addWidget(self.lblResultsSection, 5, 0, 1, 1)   
        
        #Add the result list N
        self.lblResultsN = QtGui.QLabel('N Results loaded:', wdg)
        layout.addWidget(self.lblResultsN, 6, 0, 1, 1) 
        self.leResultsN = QtGui.QLineEdit(wdg)
        self.leResultsN.setFixedWidth(100) 
        layout.addWidget(self.leResultsN, 6, 1, 1, 10) 
        self.leResultsN.setEnabled(False)
        
        #Add the result index
        self.lblResultIndex = QtGui.QLabel('Result #', wdg)
        layout.addWidget(self.lblResultIndex, 7, 0, 1, 1) 
        self.leResultIndex = QtGui.QLineEdit(wdg)
        self.leResultIndex.setFixedWidth(100) 
        layout.addWidget(self.leResultIndex, 7, 1, 1, 2) 

        self.btnSetResultPrevious = QtGui.QPushButton(' < ', wdg)
        self.btnSetResultPrevious.setToolTip('Got to the previous Result in the list...')
        self.btnSetResultPrevious.clicked.connect(self.bt_previousResultIndex)
        layout.addWidget(self.btnSetResultPrevious, 7, 2, 1, 1)
        self.btnSetResultNext = QtGui.QPushButton(' > ', wdg)
        self.btnSetResultNext.setToolTip('Go to the next Result in the list...')
        self.btnSetResultNext.clicked.connect(self.bt_nextResultIndex)
        layout.addWidget(self.btnSetResultNext, 7, 3, 1, 1)
        #Add the result name 
        self.lblResultName = QtGui.QLabel('        - name', wdg)
        layout.addWidget(self.lblResultName, 8, 0, 1, 1) 
        self.leResultSpcName = QtGui.QLineEdit(wdg)
        self.leResultSpcName.setFixedWidth(500) 
        layout.addWidget(self.leResultSpcName, 8, 1, 1, 2) 
        self.leResultSpcName.setEnabled(False)
        #Add the result zdiff
        self.lblResultZdiff = QtGui.QLabel('        - zdiff', wdg)
        layout.addWidget(self.lblResultZdiff, 9, 0, 1, 1) 
        self.leResultSpcZdiff = QtGui.QLineEdit(wdg)
        self.leResultSpcZdiff.setFixedWidth(100) 
        layout.addWidget(self.leResultSpcZdiff, 9, 1, 1, 2) 
        self.leResultSpcZdiff.setEnabled(False)
        #Add the result zcalc 
        self.lblResultZcalc = QtGui.QLabel('        - zcalc', wdg)
        layout.addWidget(self.lblResultZcalc, 10, 0, 1, 1) 
        self.leResultSpcZcalc = QtGui.QLineEdit(wdg)
        self.leResultSpcZcalc.setFixedWidth(100) 
        layout.addWidget(self.leResultSpcZcalc, 10, 1, 1, 2) 
        self.leResultSpcZcalc.setEnabled(False)
        #Add the result zref
        self.lblResultZref = QtGui.QLabel('        - zref', wdg)
        layout.addWidget(self.lblResultZref, 11, 0, 1, 1) 
        self.leResultSpcZref = QtGui.QLineEdit(wdg)
        self.leResultSpcZref.setFixedWidth(100) 
        layout.addWidget(self.leResultSpcZref, 11, 1, 1, 2) 
        self.leResultSpcZref.setEnabled(False)  
      
        #Add the show button
        self.btn = QtGui.QPushButton('Show', wdg)
        self.btn.setFixedWidth(500)
        self.btn.setFixedHeight(50)
        self.btn.clicked.connect(self.bt_showResult)
        self.btn.setToolTip('Display the Chi2/spectrum/fitted template/linemodel results successively...')
        layout.addWidget(self.btn, 15, 1, 1, 1)
        
        self.setCentralWidget(wdg)
        self.setWindowTitle('Aview')
        

        #Set path from settings
        self.settings = QtCore.QSettings("amazed", "aviewgui")
        _resDir = self.settings.value("resDir").toString()
        self.leResDir.setText(_resDir)      
        if _resDir=="":
            self.enableCtrls(False)
            
        #Set spc filter
        self.leSpcFilter.setText("Not implement yet...")
        self.leSpcFilter.setEnabled(False)
            
        #Set zdiff Threshold
        self.leDiffThres.setText("")
            
            
        self.show()
 
    def bt_loadresults(self):
        self.setCurrentDir()
    
        _resDir = str(self.leResDir.text())
        print("using _resDir: {}".format(_resDir))
        
        statsPath = os.path.join(_resDir, "stats")
        diffPath = os.path.join(statsPath, "diff.txt")
        if not os.path.exists(diffPath):
            print("diff file not found... computing the diff file is necessary: please select the reference redshift file list in order continue !")
            calcFile = os.path.join(_resDir, "redshift.csv") 
            
            _refzfilepathDefault = self.settings.value("refzfilepath").toString()
            _refzfilepath = str(QtGui.QFileDialog.getOpenFileName(self, "Select Reference Redshift list file", _refzfilepathDefault))
            if os.path.exists(_refzfilepath) :
                if os.path.isdir(statsPath)==False:
                    os.mkdir( statsPath, 0755 );
                self.settings.setValue("refzfilepath", _refzfilepath)
                processAmazedOutputStats.setPFSRefFileType()
                processAmazedOutputStats.ProcessDiff( _refzfilepath, calcFile, diffPath , reftype='pfs')
            else:
                print("Selected zref file does not exist...")
        
        _spcName = str(self.leSpcFilter.text()) #todo: filter the results by name tag
        _diffthres = str(self.leDiffThres.text())
        if _diffthres=="":
            _diffthres = -1
        
        self.resList = resultstat.ResultList(_resDir, diffthreshold=float(_diffthres), opt='brief') 
        self.leResultsN.setText(str(self.resList.n))
        #init result index to 1 (first index)
        self.leResultIndex.setText(str(1))
        self.refreshResultDetails()

    def refreshResultDetails(self):
        current_index = int(self.leResultIndex.text() )-1
        self.leResultSpcName.setText(self.resList.list[current_index].name)
        self.leResultSpcZdiff.setText(str(self.resList.list[current_index].zdiff))
        self.leResultSpcZcalc.setText(str(self.resList.list[current_index].zcalc))
        self.leResultSpcZref.setText(str(self.resList.list[current_index].zref))
        
    
    def bt_nextResultIndex(self):
        current_index = int(self.leResultIndex.text() )
        if current_index>self.resList.n-1:
            return
        self.leResultIndex.setText(str(current_index+1))  
        self.refreshResultDetails()  

    def bt_previousResultIndex(self):
        current_index = int(self.leResultIndex.text() )
        if current_index<2:
            return
        self.leResultIndex.setText(str(current_index-1))
        self.refreshResultDetails()
        
        
    #deprecated...
    def bt_showDialog(self):
        #self.setCurrentDir()
        _resDir = str(self.leResDir.text())
        print("using _resDir: {}".format(_resDir))
        
        statsPath = os.path.join(_resDir, "stats")
        diffPath = os.path.join(statsPath, "diff.txt")
        if not os.path.exists(diffPath):
            print("diff file not found... computing the diff file is necessary: please select the reference redshift file list in order continue !")
            calcFile = os.path.join(_resDir, "redshift.csv") 
            
            _refzfilepathDefault = self.settings.value("refzfilepath").toString()
            _refzfilepath = str(QtGui.QFileDialog.getOpenFileName(self, "Select Reference Redshift list file", _refzfilepathDefault))
            if os.path.exists(_refzfilepath) :
                if os.path.isdir(statsPath)==False:
                    os.mkdir( statsPath, 0755 );
                self.settings.setValue("refzfilepath", _refzfilepath)
                processAmazedOutputStats.setPFSRefFileType()
                processAmazedOutputStats.ProcessDiff( _refzfilepath, calcFile, diffPath )
            else:
                print("Selected zref file does not exist...")
        
        _spcName = str(self.leSpcFilter.text())
        _tplpath = ""
        _redshift = ""
        _iextremaredshift = 0
        _diffthres = str(self.leDiffThres.text())
        if _diffthres=="":
            _diffthres = -1
        _failureindex = "0"
        aview.plotRes(_resDir, _spcName, _tplpath, _redshift, _iextremaredshift, _diffthres, _failureindex)


    def bt_showResult(self):
        _resDir = str(self.leResDir.text())
        
        _spcName = str(self.leResultSpcName.text())
        _tplpath = ""
        
        #current_index = int(self.leResultIndex.text() )-1
        _redshift = ""#self.resList.list[current_index].zcalc

        _iextremaredshift = 0
        _diffthres = str(self.leDiffThres.text())
        if _diffthres=="":
            _diffthres = -1
        _failureindex = "0"
        aview.plotRes(_resDir, _spcName, _tplpath, _redshift, _iextremaredshift, _diffthres, _failureindex)
        
    
    def bt_setResultDir(self):
        _resDirDefault = os.path.abspath(str(self.leResDir.text()))
        _resDirDefault = _resDirDefault[:_resDirDefault.index(os.sep)] if os.sep in _resDirDefault else _resDirDefault
        
        _resDir = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Directory",_resDirDefault))
        if os.path.exists(_resDir):
            #check the diff file is present, if not, do something...
            print("_resDir = {}".format(_resDir))
            print("previous from settings = {}".format(self.settings.value("resDir").toString()))
            
            if not _resDir == self.settings.value("resDir").toString():
                _resDirPrevious = self.settings.value("resDir").toString()
                self.settings.setValue("resDirPrevious", _resDirPrevious);
            self.leResDir.setText(_resDir)
            self.settings.setValue("resDir", _resDir);
            self.enableCtrls(True)
        
    def bt_setResultDirPrevious(self):
        _resDirPrevious = self.settings.value("resDirPrevious").toString()
        self.leResDir.setText(_resDirPrevious)
        self.settings.setValue("resDir", _resDirPrevious);
        self.enableCtrls(True)
        
    
    def enableCtrls(self, val):
        self.leDiffThres.setEnabled(val)
        #self.leSpcFilter.setEnabled(val)
        self.btn.setEnabled(val)
    
    def setCurrentDir(self):
        _savedCurDir = os.path.realpath(os.path.abspath(os.path.curdir))
        _resDir = str(self.leResDir.text())
        _resParser = resparser.ResParser(_resDir)

        #fixed
        #_amazedDir = "/home/aschmitt/data/pfs/pfs2_simu20151118_jenny/amazed"    
        #os.chdir(_amazedDir)
    
        #try to find the amazed working directory by identifying the linecatalog relative path
        _abscurpath = os.path.realpath(os.path.abspath(_resDir))
        _relCatalogDir = os.path.normpath( _resParser.getConfigVal("linecatalog"))
        if len(glob.glob(_relCatalogDir))>0:
            return
        #print("_relCatalogDir = {}".format(_relCatalogDir))
        continueDown = True;
        while continueDown:
            os.chdir(_abscurpath)
            #print("\ncurdir = {}".format(_abscurpath))
            subfolder = "../"      
            
            _possible_root = glob.glob(_relCatalogDir) 
            _possible_level1 = glob.glob("*/" + _relCatalogDir)
            _possible_level2 = glob.glob("*/" + "*/" + _relCatalogDir)            
            #print("_possible_root = {}".format(_possible_root))     
            #print("_possible_level1 = {}".format(_possible_level1))
            
            if _abscurpath==os.path.sep:
                continueDown = False
                os.chdir(_savedCurDir)
            elif len(_possible_root)==1:
                continueDown = False
                print("SETTING CURRENT DIRECTORY: _amazedDir = {}".format(_abscurpath))
            elif len(_possible_level1)==1:
                subdir = _possible_level1[0].split(os.path.sep)[0]
                _amazedDir = os.path.join(_abscurpath, subdir)
                os.chdir(_amazedDir)
                continueDown = False
                print("SETTING CURRENT DIRECTORY: _amazedDir = {}".format(_amazedDir))
            elif len(_possible_level2)==1:
                subdir0 = _possible_level2[0].split(os.path.sep)[0]
                subdir1 = _possible_level2[0].split(os.path.sep)[1]
                _amazedDir = os.path.join(_abscurpath, subdir0, subdir1)
                os.chdir(_amazedDir)
                continueDown = False
                print("SETTING CURRENT DIRECTORY: _amazedDir = {}".format(_amazedDir))
            else:
                _abscurpath = os.path.realpath(os.path.abspath(os.path.join(_abscurpath,subfolder)))  
 
def main():
    app = QtGui.QApplication(sys.argv)
    ins = ObjectAViewGui()
    ex = AViewGui(ins)
    sys.exit(app.exec_())
 
if __name__ == '__main__':
    main()