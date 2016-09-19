# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 09:06:12 2016

@author: aschmitt
"""

 
from PyQt5 import QtGui, QtCore, QtWidgets
import sys
import os
import inspect
import glob
import time
import threading

import shutil

import aview
import resultstat
import resparser
import aviewgui


class ObjectAProcessGui(object):
    def __init__(self, parent=None):
        self.value = []
        self.parent = parent
        


def file_len(fname):
    if os.path.exists(fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    else:
        return 0
        
class AProcessGui(QtWidgets.QWidget):
    def __init__(self, obj=None, initReprocess=None):
        super(AProcessGui, self).__init__()
        self.obj = obj
        wdg = self#QtGui.QWidget()
        layout = QtWidgets.QGridLayout(wdg)    
        
        layoutRow = 0
        separator_height = 20        
        separator_ncols = 5
        
        #Add the configuration separator
        self.lblInputSection = QtWidgets.QLabel('Input', wdg) 
        self.lblInputSection.setAlignment(QtCore.Qt.AlignCenter)           
        self.lblInputSection.setFixedHeight(separator_height)
        self.lblInputSection.setProperty("coloredcell", True)
        layout.addWidget(self.lblInputSection, layoutRow, 0, 1, 1) 
        for i in range(1, separator_ncols):
            lbl = QtWidgets.QLabel('', wdg)
            if i==separator_ncols-1:
                lbl.setFixedWidth(100)
            lbl.setProperty("coloredcell", True)
            layout.addWidget(lbl, layoutRow, i, 1, 1)         
 
        
        #Add the spectrumdir setup ctrls
        layoutRow += 1
        self.lblapcdirPath = QtWidgets.QLabel('Spectrum directory path', wdg)
        layout.addWidget(self.lblapcdirPath, layoutRow, 0, 1, 1)
        
        self.leSpcdir = QtWidgets.QLineEdit(wdg)
        self.leSpcdir.setFixedWidth(500) 
        layout.addWidget(self.leSpcdir, layoutRow, 1, 1, 10)
        
        self.btnBrowseSpcdir = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseSpcdir.setToolTip('Browse to select the spectrum directory ...')
        self.btnBrowseSpcdir.clicked.connect(self.bt_setSpcdir)
        layout.addWidget(self.btnBrowseSpcdir, layoutRow, 2, 1, 1)
        
        self.btnBrowseSpcdirPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseSpcdirPrevious.setToolTip('Back to the previous directory path...')
        self.btnBrowseSpcdirPrevious.clicked.connect(self.bt_setSpcdirPrevious)
        layout.addWidget(self.btnBrowseSpcdirPrevious, layoutRow, 3, 1, 1)
        
        #Add the spectrumlist setup ctrls
        layoutRow += 1
        self.ckSpclist = QtWidgets.QCheckBox('Spectrum LIST path', wdg)
        self.ckSpclist.stateChanged.connect(self.ck_Spclist)
        layout.addWidget(self.ckSpclist, layoutRow, 0, 1, 1)
        
        self.leSpclist = QtWidgets.QLineEdit(wdg)
        self.leSpclist.setFixedWidth(500) 
        layout.addWidget(self.leSpclist, layoutRow, 1, 1, 10)
        
        self.btnBrowseSpclist = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseSpclist.setToolTip('Browse to select the spectrumlist file...')
        self.btnBrowseSpclist.clicked.connect(self.bt_setSpclist)
        layout.addWidget(self.btnBrowseSpclist, layoutRow, 2, 1, 1)
        
        self.btnBrowseSpclistPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseSpclistPrevious.setToolTip('Back to the previous file path...')
        self.btnBrowseSpclistPrevious.clicked.connect(self.bt_setSpclistPrevious)
        layout.addWidget(self.btnBrowseSpclistPrevious, layoutRow, 3, 1, 1) 
        
                
        #Add the spectrum fits setup ctrls
        layoutRow += 1
        self.ckSpcFits= QtWidgets.QCheckBox('Spectrum .FITS name', wdg)
        self.ckSpcFits.stateChanged.connect(self.ck_SpcFits)
        layout.addWidget(self.ckSpcFits, layoutRow, 0, 1, 1)
        
        self.leSpcFits = QtWidgets.QLineEdit(wdg)
        self.leSpcFits.setFixedWidth(400) 
        layout.addWidget(self.leSpcFits, layoutRow, 1, 1, 10)
        
        self.btnBrowseSpcFits = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseSpcFits.setToolTip('Browse to select the spectrum file...')
        self.btnBrowseSpcFits.clicked.connect(self.bt_setSpcFits)
        layout.addWidget(self.btnBrowseSpcFits, layoutRow, 2, 1, 1)
        
        self.btnBrowseSpcFitsPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseSpcFitsPrevious.setToolTip('Back to the previous file path...')
        self.btnBrowseSpcFitsPrevious.clicked.connect(self.bt_setSpcFitsPrevious)
        layout.addWidget(self.btnBrowseSpcFitsPrevious, layoutRow, 3, 1, 1) 
                        
        #Add the noise fits setup ctrls
        layoutRow += 1
        self.lblNoiseFits= QtWidgets.QLabel('     + Noise .FITS name', wdg)
        layout.addWidget(self.lblNoiseFits, layoutRow, 0, 1, 1)
        
        self.leNoiseFits = QtWidgets.QLineEdit(wdg)
        self.leNoiseFits.setFixedWidth(400) 
        layout.addWidget(self.leNoiseFits, layoutRow, 1, 1, 10)
        
        self.btnBrowseNoiseFits = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseNoiseFits.setToolTip('Browse to select the spectrum file...')
        self.btnBrowseNoiseFits.clicked.connect(self.bt_setNoiseFits)
        layout.addWidget(self.btnBrowseNoiseFits, layoutRow, 2, 1, 1)
        
        self.btnBrowseNoiseFitsPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseNoiseFitsPrevious.setToolTip('Back to the previous file path...')
        self.btnBrowseNoiseFitsPrevious.clicked.connect(self.bt_setNoiseFitsPrevious)
        layout.addWidget(self.btnBrowseNoiseFitsPrevious, layoutRow, 3, 1, 1) 
        

        #Add the configuration separator
        layoutRow += 1        
        self.lblConfigurationSection = QtWidgets.QLabel('Configuration', wdg)
        self.lblConfigurationSection.setAlignment(QtCore.Qt.AlignCenter)             
        self.lblConfigurationSection.setFixedHeight(separator_height)  
        self.lblConfigurationSection.setProperty("coloredcell", True)  
        layout.addWidget(self.lblConfigurationSection, layoutRow, 0, 1, 1)
        for i in range(1, separator_ncols):
            lbl = QtWidgets.QLabel('', wdg)
            lbl.setProperty("coloredcell", True)
            layout.addWidget(lbl, layoutRow, i, 1, 1)         
 

        #Add the import config button        
        layoutRow += 1
        self.btnImportFromConfigFile = QtWidgets.QPushButton(' Import ', wdg)
        self.btnImportFromConfigFile.setToolTip('Import config/parameters from existing config.txt file')
        self.btnImportFromConfigFile.clicked.connect(self.bt_importConfigFile)
        layout.addWidget(self.btnImportFromConfigFile, layoutRow, 0, 1, 1)  
         

        #Add the method selection        
        layoutRow += 1
        self.lblMethodCombo = QtWidgets.QLabel('Method', wdg)
        layout.addWidget(self.lblMethodCombo, layoutRow, 0, 1, 1)
        self.cbMethod = QtWidgets.QComboBox()
        self.cbMethod.addItem("linemodel")
        self.cbMethod.addItem("chisquare2solve")
        #self.cbMethod.addItem("amazed0_2")
        self.cbMethod.currentIndexChanged.connect(self.cb_method_selectionchange)
        layout.addWidget(self.cbMethod, layoutRow, 1, 1, 1)  
        
                
        #Add the linecatalog file ctrls
        layoutRow += 1
        self.lblLinecatalogPath = QtWidgets.QLabel('Lines catalog path', wdg)
        layout.addWidget(self.lblLinecatalogPath, layoutRow, 0, 1, 1)
        
        self.leLinecatalog = QtWidgets.QLineEdit(wdg)
        self.leLinecatalog.setFixedWidth(500) 
        layout.addWidget(self.leLinecatalog, layoutRow, 1, 1, 10)
        
        self.btnBrowseLinecatalog = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseLinecatalog.setToolTip('Browse to select the linecatalog file...')
        self.btnBrowseLinecatalog.clicked.connect(self.bt_setLinecatalogFilePath)
        layout.addWidget(self.btnBrowseLinecatalog, layoutRow, 2, 1, 1)
        
        self.btnBrowseLinecatalogPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseLinecatalogPrevious.setToolTip('Back to the previous file path...')
        self.btnBrowseLinecatalogPrevious.clicked.connect(self.bt_setLinecatalogFilePathPrevious)
        layout.addWidget(self.btnBrowseLinecatalogPrevious, layoutRow, 3, 1, 1)  

        #Add the lince catalog convert to air toggle
        layoutRow += 1
        self.lblConvertVac2Air = QtWidgets.QLabel('Convert VACUUM->AIR', wdg)
        layout.addWidget(self.lblConvertVac2Air, layoutRow, 0, 1, 1)
        
        self.ckConvertVac2Air = QtWidgets.QCheckBox('Enable conversion', wdg)
        self.ckConvertVac2Air.setFixedWidth(150) 
        self.ckConvertVac2Air.stateChanged.connect(self.ck_convertVac2Air)
        layout.addWidget(self.ckConvertVac2Air, layoutRow, 1, 1, 10)
        
        
        #Add the template directory setup ctrls
        layoutRow += 1
        self.lblTpldirPath = QtWidgets.QLabel('Templates directory path', wdg)
        layout.addWidget(self.lblTpldirPath, layoutRow, 0, 1, 1)
        
        self.leTpldir = QtWidgets.QLineEdit(wdg)
        self.leTpldir.setFixedWidth(500) 
        layout.addWidget(self.leTpldir, layoutRow, 1, 1, 10)
        
        self.btnBrowseTpldir = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseTpldir.setToolTip('Browse to select the spectrum directory ...')
        self.btnBrowseTpldir.clicked.connect(self.bt_setTpldir)
        layout.addWidget(self.btnBrowseTpldir, layoutRow, 2, 1, 1)
        
        self.btnBrowseTpldirPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseTpldirPrevious.setToolTip('Back to the previous directory path...')
        self.btnBrowseTpldirPrevious.clicked.connect(self.bt_setTpldirPrevious)
        layout.addWidget(self.btnBrowseTpldirPrevious, layoutRow, 3, 1, 1)

        
        #Add the parameters file ctrls
        layoutRow += 1
        self.lblMethodParametersPath = QtWidgets.QLabel('Parameters file (.json) path', wdg)
        layout.addWidget(self.lblMethodParametersPath, layoutRow, 0, 1, 1)
        
        self.leMethodParametersPath = QtWidgets.QLineEdit(wdg)
        self.leMethodParametersPath.setFixedWidth(500) 
        layout.addWidget(self.leMethodParametersPath, layoutRow, 1, 1, 10)
        
        self.btnBrowseMethodParameters = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseMethodParameters.setToolTip('Browse to select the parameters file...')
        self.btnBrowseMethodParameters.clicked.connect(self.bt_setParametersFilePath)
        layout.addWidget(self.btnBrowseMethodParameters, layoutRow, 2, 1, 1)
        
        self.btnBrowseMethodParametersPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseMethodParametersPrevious.setToolTip('Back to the previous file path...')
        self.btnBrowseMethodParametersPrevious.clicked.connect(self.bt_setParametersFilePathPrevious)
        layout.addWidget(self.btnBrowseMethodParametersPrevious, layoutRow, 3, 1, 1)        
        
        self.ckParametersAllDefaults = QtWidgets.QCheckBox('all defaults', wdg)
        self.ckParametersAllDefaults.setFixedWidth(150) 
        self.ckParametersAllDefaults.stateChanged.connect(self.ck_parametersAllDefaults)
        layout.addWidget(self.ckParametersAllDefaults, layoutRow, 4, 1, 10)
        

                
        #Add the thread count setup ctrls
        layoutRow += 1
        self.lblProcThreadCount = QtWidgets.QLabel('Processing thread count', wdg)
        layout.addWidget(self.lblProcThreadCount, layoutRow, 0, 1, 1)
        
        self.leProcThreadCount = QtWidgets.QLineEdit(wdg)
        self.leProcThreadCount.setFixedWidth(50) 
        layout.addWidget(self.leProcThreadCount, layoutRow, 1, 1, 10)
        
        #Add the processing section separator
        layoutRow += 1
        self.lblProcessingCtrlsSection = QtWidgets.QLabel('Processing', wdg)    
        self.lblProcessingCtrlsSection.setAlignment(QtCore.Qt.AlignCenter)         
        self.lblProcessingCtrlsSection.setFixedHeight(separator_height)  
        self.lblProcessingCtrlsSection.setProperty("coloredcell", True)  
        layout.addWidget(self.lblProcessingCtrlsSection, layoutRow, 0, 1, 1)
        for i in range(1, separator_ncols):
            lbl = QtWidgets.QLabel('', wdg)
            lbl.setProperty("coloredcell", True)
            layout.addWidget(lbl, layoutRow, i, 1, 1)         
 
        #Add the processing settings buttons : amazed bin selection 
        layoutRow += 1
        self.layoutProcessLeftCol= QtWidgets.QVBoxLayout()
        self.btnSetWorkDirPath = QtWidgets.QPushButton(' Workspace ', wdg)
        self.btnSetWorkDirPath.setToolTip('Set the working directory (ex: /usr/tmp/amazed)')
        self.btnSetWorkDirPath.clicked.connect(self.bt_setWorkspace)
        self.layoutProcessLeftCol.addWidget(self.btnSetWorkDirPath) 
        
        self.btnSetAmazedBinPath = QtWidgets.QPushButton(' Amazed Bin Path ', wdg)
        self.btnSetAmazedBinPath.setToolTip('Set AMAZED bin path (ex: /home/user/amazed/amazed-0.2.5)')
        self.btnSetAmazedBinPath.clicked.connect(self.bt_setAmazedBinPath)
        self.layoutProcessLeftCol.addWidget(self.btnSetAmazedBinPath) 
        
        self.btnPrintAmazedUsage = QtWidgets.QPushButton(' Amazed usage ', wdg)
        self.btnPrintAmazedUsage.setToolTip('Print AMAZED usage in the console')
        self.btnPrintAmazedUsage.clicked.connect(self.bt_printAmazedBinUsage)
        self.layoutProcessLeftCol.addWidget(self.btnPrintAmazedUsage) 
        layout.addLayout(self.layoutProcessLeftCol, layoutRow, 0, 1, 1)
        
        
        #Add the process button
        self.btnProcess = QtWidgets.QPushButton('Start Processing', wdg)
        self.btnProcess.setFixedWidth(500)
        self.btnProcess.setFixedHeight(100)
        self.btnProcess.clicked.connect(self.bt_process)
        self.btnProcess.setToolTip('Process with the selected parameters/configuration...')
        layout.addWidget(self.btnProcess, layoutRow, 1, 1, 1)

        self.btnShow = QtWidgets.QPushButton('AView', wdg)
        self.btnShow.setFixedHeight(100)
        self.btnShow.clicked.connect(self.bt_show)
        self.btnShow.setToolTip('Start aview to display the results')
        layout.addWidget(self.btnShow, layoutRow, 2, 1, 1)


        #Progress bar        
        layoutRow += 1
        self.barProgress = QtWidgets.QProgressBar(self)
        layout.addWidget(self.barProgress, layoutRow, 1, 1, 1)
        #self.progress.setGeometry(200, 80, 250, 20)

 
        #self.setCentralWidget(wdg)
        self.setLayout(layout)        
        self.setWindowTitle('AProcess')        

        
        ### INITIAL POPULATE       
        self.settings = QtCore.QSettings("amazed", "aprocessgui")  
        #self.settings.clear() 

        #some directories settings
        _workspace = self.settings.value("workspacePath", defaultValue = "/var/tmp/amazed", type=str)
        self.setWorkspace(_workspace)    
         
        _amazed_bin_path = self.settings.value("amazedBinPath", defaultValue = "/home/aschmitt/gitlab/amazed/bin/amazed-0.0.0", type=str)
        self.setAmazedBinPath(_amazed_bin_path)
        print("\n")
             
        #Set spclist path from settings
        _spclistPath = self.settings.value("spclistPath", defaultValue="", type=str)
        print("Settings: _spclistPath = {}".format(_spclistPath))
        self.setSpclist(_spclistPath)
            
        _spcDirPath = self.settings.value("spcDirPath", defaultValue="", type=str)
        print("Settings: _spcDirPath = {}".format(_spcDirPath))
        self.setSpcdir(_spcDirPath)
        
        _ck_useSpcList = self.settings.value("enableSpectrumlist")
        _ck_useSpcListBool = _ck_useSpcList=="True"
        if _ck_useSpcListBool:        
            self.ckSpclist.toggle()   
            self.ckSpclist.toggle()   
            self.ckSpclist.toggle()
        else:   
            self.ckSpcFits.toggle()
            self.ckSpcFits.toggle()
            self.ckSpcFits.toggle()
        _spcFitsName = self.settings.value("spcFitsName", defaultValue="", type=str)
        print("Settings: _spcFitsName = {}".format(_spcFitsName))
        self.setSpcFits(_spcFitsName)
        
        _noiseFitsPath = self.settings.value("noiseFitsPath", defaultValue="", type=str)
        self.leNoiseFits.setText(_noiseFitsPath) 
        
            
        #templates dir path        
        _tplDirPath = self.settings.value("tplDirPath")
        self.leTpldir.setText(_tplDirPath) 
        #Set parameters file path from settings
        _ParametersFilePath = self.settings.value("ParametersFilePath")
        self.leMethodParametersPath.setText(_ParametersFilePath)
        #set checbox all default parameters from settings
        _cb_parametersAllDefault = self.settings.value("parametersAllDefault")
        _cb_parametersAllDefaultBool = _cb_parametersAllDefault=="True"
        #print("_cb_parametersAllDefault = {}".format(_cb_parametersAllDefault))
        #self.ck_parametersAllDefaults(_cb_parametersAllDefaultBool)           
        if _cb_parametersAllDefaultBool:        
            self.ckParametersAllDefaults.toggle()
        
        #Set linecatalog file path from settings
        _linecatalogFilePath = self.settings.value("linecatalogPath")
        self.leLinecatalog.setText(_linecatalogFilePath)  
        _cb_convertVac2Air = self.settings.value("linecatalogConvertToAir")
        _cb_convertVac2AirBool = _cb_convertVac2Air=="True"
        if _cb_convertVac2AirBool:        
            self.ckConvertVac2Air.toggle()
        #set the thread count        
        self.leProcThreadCount.setText("1") 
        self.btnShow.setEnabled(False)
        
        #OPTIONALLY populate from external Result Directory and Spectrum Tag
        if not initReprocess==None:
            _resDir = initReprocess[0]
            _spcTag = initReprocess[1]
            #_resDir = "/home/aschmitt/data/keck/amazed/res_20160914_cosm4n1_linemodel_wavelets8_nosuperstrong"
            #_spcTag = "spec1d.m4-n1.005.COSMOS-1450404_F"
            self.initForReprocess( _resDir, _spcTag)
            
        wdg.setStyleSheet("*[coloredcell=\"true\"] {background-color:rgb(215,215,215);}")
        self.show()
        
    def enableCtrls(self, val):
        return
   
    def bt_setSpclist(self):
        self.setSpclist()
           
    def setSpclist(self, newPath=None):
        tag = "setSpclist"
        if newPath==None:
            _spclistDefault = os.path.abspath(str(self.leSpclist.text()))
            #_spclistDefault = _spclistDefault[:_spclistDefault.index(os.sep)] if os.sep in _spclistDefault else _spclistDefault
            
            _spclistPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select File",_spclistDefault)
        else:
            _spclistPath = newPath
            
        if os.path.exists(_spclistPath):
            #check the diff file is present, if not, do something...
            print("{}: _spclistPath = {}".format(tag, _spclistPath))
            print("{}: previous from settings = {}".format(tag, self.settings.value("spclistPath")))
            
            if not _spclistPath == self.settings.value("spclistPath"):
                _spclistPathPrevious = self.settings.value("spclistPath")
                self.settings.setValue("spclistPathPrevious", _spclistPathPrevious);
            self.leSpclist.setText(_spclistPath)
            self.settings.setValue("spclistPath", _spclistPath);
            self.enableCtrls(True)

    def bt_setSpclistPrevious(self):
        _spclistPathPrevious = self.settings.value("spclistPathPrevious")
        self.leSpclist.setText(_spclistPathPrevious)
        self.settings.setValue("spclistPath", _spclistPathPrevious);
        self.enableCtrls(True) 
        
    
    def ck_Spclist(self, state):
        if state == QtCore.Qt.Checked:                
            self.settings.setValue("enableSpectrumlist", "True");            
            self.leSpclist.setEnabled(True)   
            self.btnBrowseSpclist.setEnabled(True) 
            self.btnBrowseSpclistPrevious.setEnabled(True) 
                      
            if self.ckSpcFits.isChecked():
                self.ckSpcFits.toggle()
        else:                
            self.settings.setValue("enableSpectrumlist", "False");            
            self.leSpclist.setEnabled(False)
            self.btnBrowseSpclist.setEnabled(False) 
            self.btnBrowseSpclistPrevious.setEnabled(False)
            if not self.ckSpcFits.isChecked():
                self.ckSpcFits.toggle()
        
            
    def bt_setSpcFits(self):
        self.setSpcFits()
        
    def setSpcFits(self, newName=None):
        tag = "setSpcFits"
        
        _spcDirPath = str(self.leSpcdir.text())
            
        if newName==None:
            _SpcFitsDefault = _spcDirPath
            
            _SpcFitsPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select File",_SpcFitsDefault)
            _SpcFitsName = os.path.split(_SpcFitsPath)[1]
        else:
            _SpcFitsName = newName
            
        #print("_SpcFitsName = {}".format(_SpcFitsName))
        _SpcFitsPath = os.path.join(_spcDirPath , _SpcFitsName)
        if os.path.exists(_SpcFitsPath):
            current_settings_value = self.settings.value("spcFitsName")
            #check the diff file is present, if not, do something...
            #print("{}: _SpcFitsName = {}".format(tag, _SpcFitsName))
            #print("{}: previous from settings = {}".format(tag, current_settings_value))

            if not _SpcFitsName == current_settings_value:
                _SpcFitsNamePrevious = current_settings_value
                self.settings.setValue("spcFitsNamePrevious", _SpcFitsNamePrevious)
            self.leSpcFits.setText(_SpcFitsName)
            self.settings.setValue("spcFitsName", _SpcFitsName) 
            #print("{}: settings SpcFitsName set to = {}".format(tag, _SpcFitsName))
            
            self.enableCtrls(True)
        else:
            print("{}: Spectrum fits file not located in spectrum directory. Please update spectrum directory path first.".format(tag))
            #msg = QtGui.QMessageBox()
            #msg.setText("Spectrum fits file not located in spectrum directory. Please update spectrum directory path first.")   
            #msg.setWindowTitle("Error")	
            #msg.setWindowTitle("MessageBox demo")
            #retval = msg.exec_()


    def bt_setSpcFitsPrevious(self):
        return
        _SpcFitsPathPrevious = self.settings.value("SpcFitsPathPrevious")
        self.leSpcFits.setText(_SpcFitsPathPrevious)
        self.settings.setValue("SpcFitsPath", _SpcFitsPathPrevious);
        self.enableCtrls(True) 
        
    
    def ck_SpcFits(self, state):        
        if state == QtCore.Qt.Checked:                
            self.settings.setValue("enableSpectrumlist", "False");
            self.leSpcFits.setEnabled(True) 
            self.btnBrowseSpcFits.setEnabled(True) 
            self.btnBrowseSpcFitsPrevious.setEnabled(True) 
            self.leNoiseFits.setEnabled(True) 
            self.btnBrowseNoiseFits.setEnabled(True) 
            self.btnBrowseNoiseFitsPrevious.setEnabled(True) 
            
            if self.ckSpclist.isChecked():
                self.ckSpclist.toggle()
        else:                
            self.settings.setValue("enableSpectrumlist", "True");            
            self.leSpcFits.setEnabled(False)
            self.btnBrowseSpcFits.setEnabled(False) 
            self.btnBrowseSpcFitsPrevious.setEnabled(False)            
            self.leNoiseFits.setEnabled(False)
            self.btnBrowseNoiseFits.setEnabled(False) 
            self.btnBrowseNoiseFitsPrevious.setEnabled(False) 
            
            if not self.ckSpclist.isChecked():
                self.ckSpclist.toggle()
            
   
    def bt_setNoiseFits(self):
        self.setNoiseFits()
           
    def setNoiseFits(self, newName=None):
        tag = "setNoiseFits"
        
        _spcDirPath = str(self.leSpcdir.text())
        
        if newName==None:
            _NoiseFitsDefault = _spcDirPath
            
            _NoiseFitsPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select File",_NoiseFitsDefault)
            _NoiseFitsName = os.path.split(_NoiseFitsPath)[1]
        else:
            _NoiseFitsName = newName
            
        _NoiseFitsPath = os.path.join(_spcDirPath , _NoiseFitsName)
        if os.path.exists(_NoiseFitsPath):
            current_settings_value = self.settings.value("noiseFitsPath")
            #check the diff file is present, if not, do something...
            print("{}: _NoiseFitsName = {}".format(tag, _NoiseFitsName))
            print("{}: previous from settings = {}".format(tag, current_settings_value))
            
            
            if not _NoiseFitsPath == current_settings_value:
                _NoiseFitsPathPrevious = current_settings_value
                self.settings.setValue("noiseFitsPathPrevious", _NoiseFitsPathPrevious);
            self.leNoiseFits.setText(_NoiseFitsName)
            self.settings.setValue("noiseFitsPath", _NoiseFitsName);
            self.enableCtrls(True)
        else:
            print("{}: Noise fits file not located in spectrum directory. Please update spectrum directory path first.".format(tag))
            msg = QtGui.QMessageBox()
            msg.setText("Noise fits file not located in spectrum directory. Please update spectrum directory path first.")   
            msg.setWindowTitle("Error")	
            msg.setWindowTitle("MessageBox demo")
            retval = msg.exec_()

    def bt_setNoiseFitsPrevious(self):
        _NoiseFitsPathPrevious = self.settings.value("NoiseFitsPathPrevious")
        self.leNoiseFits.setText(_NoiseFitsPathPrevious)
        self.settings.setValue("NoiseFitsPath", _NoiseFitsPathPrevious);
        self.enableCtrls(True) 
        
   
    def bt_setSpcdir(self):
        self.setSpcdir()
        
    def setSpcdir(self, newDir=None):    
        tag = "setSpcDir"       
        
        if newDir==None:
            _spcdirDefault = str(self.leSpcdir.text())
            if not os.path.exists(_spcdirDefault):
                _spcdirDefault = os.path.abspath(str(self.leSpcdir.text()))
                #_spcdirDefault = _spcdirDefault[:_spcdirDefault.index(os.sep)] if os.sep in _spcdirDefault else _spcdirDefault
                
            _spcdirPath = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Spectrum Directory",_spcdirDefault)
        else:
            _spcdirPath = newDir
            
        if os.path.exists(_spcdirPath):
            #check the diff file is present, if not, do something...
            print("{}: _spcdirPath = {}".format(tag, _spcdirPath))
            print("{}: previous from settings = {}".format(tag, self.settings.value("spcdirPath")))
            
            if not _spcdirPath == self.settings.value("spcdirPath"):
                _spcdirPathPrevious = self.settings.value("spcdirPath")
                self.settings.setValue("spcdirPathPrevious", _spcdirPathPrevious);
            self.leSpcdir.setText(_spcdirPath)
            self.settings.setValue("spcDirPath", _spcdirPath);
            self.enableCtrls(True)

    def bt_setSpcdirPrevious(self):
        _spcDirPathPrevious = self.settings.value("spcDirPathPrevious")
        self.leSpcdir.setText(_spcDirPathPrevious)
        self.settings.setValue("spcDirPath", _spcDirPathPrevious);
        self.enableCtrls(True) 

        
    def initForReprocess(self, resDir, spcTag):
        """
        Initialize the dialog/ctrls with custom parameters coming from an existing batch.
        """
        #self.setInputSpcToReprocess()
        _resParser = resparser.ResParser(resDir)

        #retrieve config params to ui controls
        _configPath = _resParser.configpath        
        self.importConfigFile(_configPath)
        #override parameters .json from the resdir
        newParamsPath = os.path.join(self.m_workspace, "params_reprocess.json")
        shutil.copyfile(_resParser.parameterspath, newParamsPath)
        self.setParametersFilePath(newParamsPath)

        #retrieve spectrum directory to ui controls
        _spcDir = _resParser.getConfigVal('spectrumdir')
        self.setSpcdir(_spcDir)
        
        #retrieve/set spectrum/noise fits path
        splListLine = _resParser.getSpectrumlistline(spcTag)
        _SpcFits = splListLine[0]
        self.setSpcFits(_SpcFits)
        _SpcNoise = splListLine[1]
        self.setNoiseFits(_SpcNoise)
           
       
    def bt_importConfigFile(self):
                
        _configFilePathDefault = "/home/aschmitt/Documents/amazed/amazed_v0.2_implementation/data-0.2.4/config_example.txt"
        #_TpldirDefault = _TpldirDefault[:_TpldirDefault.index(os.sep)] if os.sep in _TpldirDefault else _TpldirDefault
        _configFilePath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select config.txt file",_configFilePathDefault)
        
        self.importConfigFile(_configFilePath)
               
    def importConfigFile(self, configPath):
        """
        Import the config settings from a file to the UI controls. (except for inputs (spectrum dir, spectrum list) that are not imported)
        """
        tag = "bt_importConfigFile"

        _configFilePath = configPath
        _configFileDir = os.path.split(_configFilePath)[0]
        print("config dir = {}".format(_configFileDir))
        
        if not os.path.exists(_configFilePath):
            print("{}: ERROR, config file does not exist!".format(tag))
            return
            
        print("{}: importing from {}".format(tag, _configFilePath))
        convertToAirFound = False
        f = open(_configFilePath)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                valStr = ""
                #print lineStr
                keyword = "method="
                if keyword in lineStr:
                    valStr = lineStr.replace(keyword, "")
                    print("{}: Found {} {}".format(tag, keyword, valStr))
                    methodSet = False
                    for count in range(self.cbMethod.count()):
                        if valStr.lower() == str(self.cbMethod.itemText(count)).lower():
                            self.cbMethod.setCurrentIndex(count)
                            methodSet = True
                            break
                    if not methodSet:
                        print("{}: ERROR unable to set the method ({}) from config file".format(tag, valStr))
                            
                keyword = "linecatalog="
                if keyword in lineStr:
                    valStr = lineStr.replace(keyword, "")
                    print("{}: Found {} {}".format(tag, keyword, valStr))
                    lcatPath = os.path.abspath(os.path.join(_configFileDir, valStr))
                    self.bt_setLinecatalogFilePath(lcatPath)
                keyword = "linecatalog-convert="
                if keyword in lineStr:
                    valStr = lineStr.replace(keyword, "")
                    print("{}: Found {} {}".format(tag, keyword, valStr))
                    convertToAirFound = True
                    
                keyword = "parameters="
                if keyword in lineStr:
                    valStr = lineStr.replace(keyword, "")
                    print("{}: Found {} {}".format(tag, keyword, valStr))
                    paramPath = os.path.abspath(os.path.join(_configFileDir, valStr))
                    self.setParametersFilePath(paramPath)
                keyword = "templatedir="
                if keyword in lineStr:
                    valStr = lineStr.replace(keyword, "")
                    print("{}: Found {} {}".format(tag, keyword, valStr))
                    tplPath = os.path.abspath(os.path.join(_configFileDir, valStr))
                    self.setTpldir(tplPath)
        f.close()
        
        if convertToAirFound:
            self.settings.setValue("linecatalogConvertToAir", "True");
            if not self.ckConvertVac2Air.isChecked():        
                self.ckConvertVac2Air.toggle()
        else:
            self.settings.setValue("linecatalogConvertToAir", "False");
            if self.ckConvertVac2Air.isChecked():        
                self.ckConvertVac2Air.toggle()

   
    def bt_setTpldir(self):
        self.setTpldir()
        
    def setTpldir(self, newPath=None):
        tag = "setTpldir"
        if newPath==None:
            _TpldirDefault = os.path.abspath(str(self.leTpldir.text()))
            #_TpldirDefault = _TpldirDefault[:_TpldirDefault.index(os.sep)] if os.sep in _TpldirDefault else _TpldirDefault
            
            _TpldirPath = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Templates directory",_TpldirDefault)
        else:
            _TpldirPath = newPath
            
        if os.path.exists(_TpldirPath):
            current_settings_value = self.settings.value("tplDirPath")
            #check the diff file is present, if not, do something...
            print("{}: _TpldirPath = {}".format(tag, _TpldirPath))
            print("{}: previous from settings = {}".format(tag, current_settings_value))
            
            if not _TpldirPath == current_settings_value:
                _TpldirPathPrevious = current_settings_value
                self.settings.setValue("TpldirPathPrevious", _TpldirPathPrevious);
            self.leTpldir.setText(_TpldirPath)
            self.settings.setValue("tplDirPath", _TpldirPath);
            self.enableCtrls(True)

    def bt_setTpldirPrevious(self):
        _TpldirPathPrevious = self.settings.value("TpldirPathPrevious")
        self.leTpldir.setText(_TpldirPathPrevious)
        self.settings.setValue("tplDirPath", _TpldirPathPrevious);
        self.enableCtrls(True) 
        
        
    def bt_setParametersFilePath(self, newPath=None):
        self.setParametersFilePath()
                
    def setParametersFilePath(self, newPath=None):
        tag = "setParametersFilePath"
        if newPath==None:
            _MethodParametersDefault = os.path.abspath(str(self.leMethodParametersPath.text()))
            #_MethodParametersDefault = _MethodParametersDefault[:_MethodParametersDefault.index(os.sep)] if os.sep in _MethodParametersDefault else _MethodParametersDefault
            
            _MethodParametersPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select File",_MethodParametersDefault)
        else:
            _MethodParametersPath = newPath
            
        if os.path.exists(_MethodParametersPath):
            #check the diff file is present, if not, do something...
            print("{}: _MethodParametersPath = {}".format(tag, _MethodParametersPath))
            print("{}: previous from settings = {}".format(tag, self.settings.value("ParametersFilePath")))
            
            if not _MethodParametersPath == self.settings.value("ParametersFilePath"):
                _MethodParametersPathPrevious = self.settings.value("ParametersFilePath")
                self.settings.setValue("ParametersFilePathPrevious", _MethodParametersPathPrevious);
            self.leMethodParametersPath.setText(_MethodParametersPath)
            self.settings.setValue("ParametersFilePath", _MethodParametersPath);
            self.enableCtrls(True)

    def bt_setParametersFilePathPrevious(self):
        _MethodParametersPathPrevious = self.settings.value("ParametersFilePathPrevious")
        self.leMethodParametersPath.setText(_MethodParametersPathPrevious)
        self.settings.setValue("ParametersFilePath", _MethodParametersPathPrevious);
        self.enableCtrls(True) 
        
    def ck_parametersAllDefaults(self, state):
        if state == QtCore.Qt.Checked:
            self.leMethodParametersPath.setEnabled(False)        
            self.btnBrowseMethodParameters.setEnabled(False)        
            self.btnBrowseMethodParametersPrevious.setEnabled(False)                
            self.settings.setValue("parametersAllDefault", "True");
        else:
            self.leMethodParametersPath.setEnabled(True)         
            self.btnBrowseMethodParameters.setEnabled(True)        
            self.btnBrowseMethodParametersPrevious.setEnabled(True)                 
            self.settings.setValue("parametersAllDefault", "False");
            
        
    def bt_setLinecatalogFilePath(self, newPath=None):
        tag = "setLinecatalogFilePath"
        if newPath==None:
            _linecatalogDefault = os.path.abspath(str(self.leLinecatalog.text()))
            #_linecatalogDefault = _linecatalogDefault[:_linecatalogDefault.index(os.sep)] if os.sep in _linecatalogDefault else _linecatalogDefault
            #print("default path = {}".format(_linecatalogDefault))
                    
            _linecatalogPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select File",_linecatalogDefault)
        else:
            _linecatalogPath = newPath
            
        if os.path.exists(_linecatalogPath):
            #check the diff file is present, if not, do something...
            print("{}: _linecatalogPath = {}".format(tag, _linecatalogPath))
            print("{}: previous from settings = {}".format(tag, self.settings.value("linecatalogPath")))
            
            if not _linecatalogPath == self.settings.value("linecatalogPath"):
                _linecatalogPathPrevious = self.settings.value("linecatalogPath")
                self.settings.setValue("linecatalogPathPrevious", _linecatalogPathPrevious);
            self.leLinecatalog.setText(_linecatalogPath)
            self.settings.setValue("linecatalogPath", _linecatalogPath);
            self.enableCtrls(True)

    def bt_setLinecatalogFilePathPrevious(self):
        _linecatalogPathPrevious = self.settings.value("linecatalogPathPrevious")
        self.leLinecatalog.setText(_linecatalogPathPrevious)
        self.settings.setValue("linecatalogPath", _linecatalogPathPrevious);
        self.enableCtrls(True)

    def ck_convertVac2Air(self, state):
        if state == QtCore.Qt.Checked:                
            self.settings.setValue("linecatalogConvertToAir", "True");
        else:                
            self.settings.setValue("linecatalogConvertToAir", "False");
      
    def cb_method_selectionchange(self,i):
        tag = "cb_method_selectionchange"
        print("{}, Items in the list are :".format(tag))
        for count in range(self.cbMethod.count()):
            print("    {}".format(self.cbMethod.itemText(count)))
        print("{}: Current index = {}, selection changed to {}".format(tag, i, self.cbMethod.currentText()))

    def bt_setWorkspace(self, newWorkspace=None):
        self.setWorkspace()
        
    def setWorkspace(self, newWorkspace=None):
        tag = "bt_setWorkspace"
        if newWorkspace==None:
            _wokspaceDefault = self.m_workspace
            
            _workspace = QtWidgets.QFileDialog.getExistingDirectory(self, "Select working directory",_wokspaceDefault)
        else:
            _workspace = newWorkspace
        if os.path.exists(_workspace):
            self.m_workspace = _workspace
            self.settings.setValue("workspacePath", self.m_workspace);
        else:
            print("{}: ERROR: Unable to set workspace {}".format(tag, _workspace))
        
        print("{}: Using workspace : {}".format(tag, self.m_workspace))       
    
    def bt_setAmazedBinPath(self, newBinPath=None):
        self.setAmazedBinPath()
        
    def setAmazedBinPath(self, newBinPath=None):
        tag = "bt_setAmazedBinPath"
        if newBinPath==None:
            _amazedBinPathDefault = self.amazed_bin_path
            #_TpldirDefault = _TpldirDefault[:_TpldirDefault.index(os.sep)] if os.sep in _TpldirDefault else _TpldirDefault
            
            _amazedBinPath, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select AMAZED binary (amazed-x-x-x)",_amazedBinPathDefault)
        else:
            _amazedBinPath = newBinPath
            
        if os.path.exists(_amazedBinPath):
            self.amazed_bin_path = _amazedBinPath
            self.settings.setValue("amazedBinPath", self.amazed_bin_path);
        else:
            print("{}: ERROR: Unable to set amazed_bin_path {}".format(tag, _amazedBinPath))
        
        print("{}: Using amazed_bin_path : {}".format(tag, self.amazed_bin_path))
        
    def bt_printAmazedBinUsage(self):
        tag = "bt_printAmazedBinUsage"
        bashCommand = "{} -h".format(self.amazed_bin_path)
        os.system(bashCommand) 

    def bt_process(self):
        tag = ""
        self.completed = 0  

        #init workspace  
        enable_reinit_workspace= False     
        if enable_reinit_workspace:
            if os.path.exists(self.m_workspace):
                shutil.rmtree(self.m_workspace) 
            os.mkdir(self.m_workspace) 
        else:
            if not os.path.exists(self.m_workspace):
                os.mkdir(self.m_workspace)
        
        #init output path
        self.m_outputPath = os.path.join(self.m_workspace, "output") 
                
        if self.ckSpclist.isChecked():
            _spclistPath = str(self.leSpclist.text())
        else:
            _spclistPath = self.createSpcListFromSpcNoiseFits()
        self.saveEmptyRefFile(_spclistPath)

            
        _configFilePath = self.saveConfigFile()        
        
        while self.completed <= 1.0: #run the first x percent initially
            self.completed += 0.0001
            self.barProgress.setValue(self.completed)
        time.sleep(0.1)
        
        self.execAmazedThread(_configFilePath)    
    
        while self.completed<99.999:
            sleep_time_seconds = 10.0
            #print("Sleeping for {} seconds...".format(sleep_time_seconds))
            print(".")
            time.sleep(sleep_time_seconds)          
            self.completed = self.getProgress()
            self.barProgress.setValue(self.completed)
            
        
        self.btnShow.setEnabled(True)
        print("Process finished.")
        return

    def msgbtn(i):
        print("Button pressed is:".format(i.text()))        
    
    def getProgress(self):

        config_output_fname = os.path.join(self.m_workspace, "output/config.txt".format())  
        if os.path.exists(config_output_fname):
            return 100.0
        
        if self.ckSpclist.isChecked():
            _spclistPath = str(self.leSpclist.text())
        else:
            _spclistPath = self.createSpcListFromSpcNoiseFits()
            
        calc_fname = os.path.join(self.m_workspace, "output/redshift.csv".format())
        nProcessedSpectra = file_len(calc_fname)
        nToProcessSpectra = file_len(_spclistPath)
        progressPercentage = float(nProcessedSpectra)/float(nToProcessSpectra)*100.0
        #print("nToProcessSpectra = {}".format(nToProcessSpectra))
        #print("nProcessedSpectra = {}".format(nProcessedSpectra))
        print("Progress is : {}/{}".format(nProcessedSpectra, nToProcessSpectra))
        return progressPercentage
        
    def createSpcListFromSpcNoiseFits(self):
        tag = "createSpcListFromSpcNoiseFits"

        _spcPath = str(self.leSpcFits.text())
        _noisePath = str(self.leNoiseFits.text())


        spclistPath = os.path.join(self.m_workspace, "spclist.spectrumlist")
        f = open(spclistPath, 'w')
        f.write("{}\t{}".format(_spcPath, _noisePath))
        f.close()
        
        return spclistPath
        
    def prepareConfigString(self):
        tag = "prepareConfigString"
        _str = ""
        
        #add the input line
        if self.ckSpclist.isChecked():
            _spclistPath = str(self.leSpclist.text())
        else:
            _spclistPath = self.createSpcListFromSpcNoiseFits()
        if os.path.exists(_spclistPath):
            _str += "input={}\n".format(os.path.abspath(_spclistPath))
        else:
            print("{}: ERROR: Spectrumlist file path is not valid, aborting...".format(tag))
        #add the input directory
        _spcDirPath = str(self.leSpcdir.text())
        if os.path.exists(_spcDirPath):
            _str += "spectrumdir={}\n".format(os.path.abspath(_spcDirPath))
        else:
            print("{}: ERROR: Spectrum directory path is not valid, aborting...".format(tag))
            
        
        #add the method line
        _method = str(self.cbMethod.currentText())
        _str += "method={}\n".format(_method)
        
        #add the linecatalog line
        _linecatalogPath = str(self.leLinecatalog.text())
        if os.path.exists(_linecatalogPath):
            _str += "linecatalog={}\n".format(os.path.abspath(_linecatalogPath))
        else:
            print("{}: ERROR: linecatalog file path is not valid, aborting...".format(tag))
        #add the line catalog conversion to AIR
        if self.ckConvertVac2Air.isChecked():
            _str += "linecatalog-convert=\n"
            
            
        #add the templates directory
        _tplDirPath = str(self.leTpldir.text())
        if os.path.exists(_tplDirPath):
            _str += "templatedir={}\n".format(os.path.abspath(_tplDirPath))
        else:
            print("{}: ERROR: Templates directory path is not valid, aborting...".format(tag))
            
        #add the parameters file line
        if not self.ckParametersAllDefaults.isChecked():
            _parametersPath = str(self.leMethodParametersPath.text())
            if os.path.exists(_parametersPath):
                _str += "parameters={}\n".format(os.path.abspath(_parametersPath))
            else:
                print("{}: ERROR: linecatalog file path is not valid, aborting...".format(tag))
        
        #add the ourput dir 
        _str += "output={}\n".format(self.m_outputPath)       
        
        
        #add the thread count line
        _thread_count = str(self.leProcThreadCount.text())
        _str += "thread-count={}\n".format(_thread_count)       

        
        return _str
        
    def saveConfigFile(self):
        tag = "saveConfigFile"
        _str = self.prepareConfigString()
        configFilePath = os.path.join(self.m_workspace, "config_tmp.txt")
        f = open(configFilePath, 'w')
        f.write(_str)
        f.close()
        return configFilePath
        
    def saveEmptyRefFile(self, spcListFilePath):
        tag = "saveEmptyRefFile"
        _strRef = ""
        
        dataNames = []
        f = open(spcListFilePath)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                data = lineStr.split("\t")
                if len(data)<1:
                    data = lineStr.split(" ")
                data = [r for r in data if r != '']
                #print len(data)
                dataname = str(data[0]).replace(".fits", "").replace(".FITS", "")
                dataNames.append(dataname)
                _strRef += "{}\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n".format(dataname)
        f.close()
        
        fPath = os.path.join(self.m_workspace, "ref_tmp.txt")
        fref = open(fPath, 'w')
        fref.write(_strRef)
        fref.close()
        return fPath
        
    def execAmazedThread(self, config_fname):
        threads = []
        t = threading.Thread(target=self.execAmazedFcn, args=(config_fname,))
        threads.append(t)
        t.start()
    
    def execAmazedFcn(self, config_fname):
        bashCommand = "{} -c {}".format(self.amazed_bin_path, config_fname)#self.amazed_bin_path = "~/gitlab/amazed/bin/amazed-0.0.0"
        os.system(bashCommand) 
        
    def bt_show(self):
        frefpath = os.path.join(self.m_workspace, "ref_tmp.txt")
        self.aviewWindow = aviewgui.AViewGui(parent=None, init_resdir=self.m_outputPath, init_refpath=frefpath)
        self.aviewWindow.show()
        return
        
def main():
    print("\n")
    app = QtWidgets.QApplication(sys.argv)
    ins = ObjectAProcessGui(parent=app)
    ex = AProcessGui(ins)
    sys.exit(app.exec_())
 
if __name__ == '__main__':
    main()
        
        