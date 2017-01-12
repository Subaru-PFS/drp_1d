# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 11:21:31 2015

@author: aschmitt
"""
 

#import mpl only to force Qt5Agg
import matplotlib as mpl
mpl.use('Qt5Agg')

from PyQt5 import QtGui, QtCore, QtWidgets
import sys
import os
import inspect
import glob
import time
import traceback

import aview
import aviewwidget
import aprocessgui
import resultstat
import resparser

subfolder = "../stats"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],subfolder)))
if cmd_subfolder not in sys.path:
    #print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
    sys.path.insert(0, cmd_subfolder)
     
import processAmazedOutputStats
 
class ObjectAViewGui(object):
    def __init__(self, parent=None):
        self.value = []
        self.parent = parent
 
 
class AViewGui(QtWidgets.QWidget):
    def __init__(self, parent=None, init_resdir="", init_refpath=""):
        super(AViewGui, self).__init__(parent)
        self.parent = parent
        wdg = self#QtGui.QWidget()
        layout = QtWidgets.QGridLayout(wdg)     
        
        layoutRow = 0
        separator_height = 20        
        separator_ncols = 4
        
        #Add the configuration separator
        self.lblInputSection = QtWidgets.QLabel('Load results', wdg) 
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
 
 
        #Add the result dir. setup ctrls
        layoutRow += 1
        self.lblResdir = QtWidgets.QLabel(' AMAZED Output Result Dir. ', wdg)
        layout.addWidget(self.lblResdir, layoutRow, 0, 1, 1)
        
        self.leResDir = QtWidgets.QLineEdit(wdg)
        self.leResDir.setFixedWidth(500) 
        layout.addWidget(self.leResDir, layoutRow, 1, 1, 10)
        
        self.btnBrowseResdir = QtWidgets.QPushButton(' Browse ', wdg)
        self.btnBrowseResdir.setToolTip('Browse to select the AMAZED output directory...')
        self.btnBrowseResdir.clicked.connect(self.bt_setResultDir)
        layout.addWidget(self.btnBrowseResdir, layoutRow, 2, 1, 1)
        
        self.btnBrowseResdirPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnBrowseResdirPrevious.setToolTip('Back to the previous ResDir...')
        self.btnBrowseResdirPrevious.clicked.connect(self.bt_setResultDirPrevious)
        layout.addWidget(self.btnBrowseResdirPrevious, layoutRow, 3, 1, 1)
        
        #Add the favorite results ctrls
        layoutRow += 1
        layoutFavorites = QtWidgets.QHBoxLayout()        
        self.btnSetResdirFavorite = QtWidgets.QPushButton(' Set as Favorite ', wdg)
        self.btnSetResdirFavorite.setToolTip('Set the favorite ResDir...')
        self.btnSetResdirFavorite.clicked.connect(self.bt_setResultDirFavorite)
        layoutFavorites.addWidget(self.btnSetResdirFavorite)
        
        self.btnLoadResdirFavorite = QtWidgets.QPushButton(' Retrieve Favorite ', wdg)
        self.btnLoadResdirFavorite.setToolTip('Retrieve favorite ResDir...')
        self.btnLoadResdirFavorite.clicked.connect(self.bt_loadResultDirFavorite)
        layoutFavorites.addWidget(self.btnLoadResdirFavorite)
        layout.addLayout(layoutFavorites, layoutRow, 1, 1, 1)
        
        #Add the failure threshold filter ctrls
        layoutRow += 1
        self.lblDiffThres = QtWidgets.QLabel(' Filter by z diff. threshold: ', wdg)
        layout.addWidget(self.lblDiffThres, layoutRow, 0, 1, 1)
        
        self.leDiffThres = QtWidgets.QLineEdit(wdg)
        self.leDiffThres.setToolTip('Enter threshold value to filter the dataset by the z error value (ex. 0.01)...')
        self.leDiffThres.setFixedWidth(200) 
        layout.addWidget(self.leDiffThres, layoutRow, 1, 1, 10)  
        
        #Add the zref range threshold filter ctrls
        layoutRow += 1
        self.lblzrefrange = QtWidgets.QLabel(' Filter by zref range: ', wdg)
        layout.addWidget(self.lblzrefrange, layoutRow, 0, 1, 1)
        
        layoutZrange = QtWidgets.QHBoxLayout() 
        self.lezrefrangemin = QtWidgets.QLineEdit(wdg)
        self.lezrefrangemin.setToolTip('Enter zref min value...')
        self.lezrefrangemin.setFixedWidth(100) 
        layoutZrange.addWidget(self.lezrefrangemin)

        self.lezrefrangemax = QtWidgets.QLineEdit(wdg)
        self.lezrefrangemax.setToolTip('Enter zref max value...')
        self.lezrefrangemax.setFixedWidth(100) 
        layoutZrange.addWidget(self.lezrefrangemax)
        #add empty label to fill the hbox
        lbl = QtWidgets.QLabel(' ', wdg)
        layoutZrange.addWidget(lbl)
        layout.addLayout(layoutZrange, layoutRow, 1, 1, 10) 
        
        #Add the sfrref range threshold filter ctrls
        layoutRow += 1
        self.lblsfrrange = QtWidgets.QLabel(' Filter by sfr range: ', wdg)
        layout.addWidget(self.lblsfrrange, layoutRow, 0, 1, 1)
        
        layoutSFRrange = QtWidgets.QHBoxLayout() 
        self.lesfrrangemin = QtWidgets.QLineEdit(wdg)
        self.lesfrrangemin.setToolTip('Enter sfr min value...')
        self.lesfrrangemin.setFixedWidth(100) 
        layoutSFRrange.addWidget(self.lesfrrangemin)

        self.lesfrrangemax = QtWidgets.QLineEdit(wdg)
        self.lesfrrangemax.setToolTip('Enter sfr max value...')
        self.lesfrrangemax.setFixedWidth(100) 
        layoutSFRrange.addWidget(self.lesfrrangemax)
        #add empty label to fill the hbox
        lbl = QtWidgets.QLabel(' ', wdg)
        layoutSFRrange.addWidget(lbl)
        layout.addLayout(layoutSFRrange, layoutRow, 1, 1, 10) 
        
        #Add the spectrum filter ctrls
        layoutRow += 1
        self.lblSpcFilter = QtWidgets.QLabel(' Filter by spectrum name: ', wdg)
        layout.addWidget(self.lblSpcFilter, layoutRow, 0, 1, 1)
        
        self.leSpcFilter = QtWidgets.QLineEdit(wdg)
        self.leSpcFilter.setToolTip('Enter a tag to filter the dataset by the spectrum name...')
        self.leSpcFilter.setFixedWidth(300) 
        layout.addWidget(self.leSpcFilter, layoutRow, 1, 1, 10)
        
        #Add the show button
        layoutRow += 1
        self.btnLoad = QtWidgets.QPushButton('Load Result List', wdg)
        self.btnLoad.setFixedWidth(500)
        self.btnLoad.setFixedHeight(50)
        self.btnLoad.clicked.connect(self.bt_loadresults)
        layout.addWidget(self.btnLoad, layoutRow, 1, 1, 1)
        
        #Add the result list N
        layoutRow += 1
        self.lblResultsN = QtWidgets.QLabel('N Results loaded:', wdg)
        layout.addWidget(self.lblResultsN, layoutRow, 0, 1, 1) 
        self.leResultsN = QtWidgets.QLineEdit(wdg)
        self.leResultsN.setFixedWidth(100) 
        layout.addWidget(self.leResultsN, layoutRow, 1, 1, 10) 
        self.leResultsN.setEnabled(False)
        
        ######### Result Details Section ###############################################
        #Add the result list section separator
        layoutRow += 1
        self.lblResultListSection = QtWidgets.QLabel('Browse results', wdg) 
        self.lblResultListSection.setAlignment(QtCore.Qt.AlignCenter)           
        self.lblResultListSection.setFixedHeight(separator_height)
        self.lblResultListSection.setProperty("coloredcell", True)
        layout.addWidget(self.lblResultListSection, layoutRow, 0, 1, 1) 
        for i in range(1, separator_ncols):
            lbl = QtWidgets.QLabel('', wdg)
            if i==separator_ncols-1:
                lbl.setFixedWidth(100)
            lbl.setProperty("coloredcell", True)
            layout.addWidget(lbl, layoutRow, i, 1, 1)   
        
        #Add the result index
        layoutRow += 1
        self.lblResultIndex = QtWidgets.QLabel('Result #', wdg)
        layout.addWidget(self.lblResultIndex, layoutRow, 0, 1, 1) 
        self.leResultIndex = QtWidgets.QLineEdit(wdg)
        self.leResultIndex.editingFinished.connect(self.leResultIndex_handleEditingFinished)
        self.leResultIndex.setFixedWidth(100) 
        layout.addWidget(self.leResultIndex, layoutRow, 1, 1, 2) 

        self.btnSetResultPrevious = QtWidgets.QPushButton(' < ', wdg)
        self.btnSetResultPrevious.setToolTip('Got to the previous Result in the list...')
        self.btnSetResultPrevious.clicked.connect(self.bt_previousResultIndex)
        layout.addWidget(self.btnSetResultPrevious, layoutRow, 2, 1, 1)
        self.btnSetResultNext = QtWidgets.QPushButton(' > ', wdg)
        self.btnSetResultNext.setToolTip('Go to the next Result in the list...')
        self.btnSetResultNext.clicked.connect(self.bt_nextResultIndex)
        layout.addWidget(self.btnSetResultNext, layoutRow, 3, 1, 1)
        #Add the result name 
        layoutRow += 1
        self.lblResultName = QtWidgets.QLabel('        - name', wdg)
        layout.addWidget(self.lblResultName, layoutRow, 0, 1, 1) 
        self.leResultSpcName = QtWidgets.QLineEdit(wdg)
        self.leResultSpcName.setFixedWidth(500) 
        layout.addWidget(self.leResultSpcName, layoutRow, 1, 1, 2) 
        self.leResultSpcName.setEnabled(False)
        #Add the result zdiff
        layoutRow += 1
        self.lblResultZdiff = QtWidgets.QLabel('        - zdiff', wdg)
        layout.addWidget(self.lblResultZdiff, layoutRow, 0, 1, 1) 
        self.leResultSpcZdiff = QtWidgets.QLineEdit(wdg)
        self.leResultSpcZdiff.setFixedWidth(300) 
        layout.addWidget(self.leResultSpcZdiff, layoutRow, 1, 1, 2) 
        self.leResultSpcZdiff.setEnabled(False)
        #Add the result zcalc 
        layoutRow += 1
        self.lblResultZcalc = QtWidgets.QLabel('        - zcalc', wdg)
        layout.addWidget(self.lblResultZcalc, layoutRow, 0, 1, 1) 
        self.leResultSpcZcalc = QtWidgets.QLineEdit(wdg)
        self.leResultSpcZcalc.setFixedWidth(100) 
        layout.addWidget(self.leResultSpcZcalc, layoutRow, 1, 1, 2) 
        self.leResultSpcZcalc.setEnabled(False)
        #Add the result zref
        layoutRow += 1
        self.lblResultZref = QtWidgets.QLabel('        - zref', wdg)
        layout.addWidget(self.lblResultZref, layoutRow, 0, 1, 1) 
        self.leResultSpcZref = QtWidgets.QLineEdit(wdg)
        self.leResultSpcZref.setFixedWidth(100) 
        layout.addWidget(self.leResultSpcZref, layoutRow, 1, 1, 2) 
        self.leResultSpcZref.setEnabled(False) 
        #Mag, sfr, ebmv info        
        layoutRow += 1
        layout_magsfrebmv = QtWidgets.QHBoxLayout()   
        #Add the result mag
        self.lblResultMagref = QtWidgets.QLabel(' mag =', wdg)
        layout_magsfrebmv.addWidget(self.lblResultMagref)
        self.leResultSpcMagref = QtWidgets.QLineEdit(wdg)
        self.leResultSpcMagref.setFixedWidth(100) 
        layout_magsfrebmv.addWidget(self.leResultSpcMagref) 
        self.leResultSpcMagref.setEnabled(False)  
        #Add the result sfr
        self.lblResultSfrref = QtWidgets.QLabel(' sfr =', wdg)
        layout_magsfrebmv.addWidget(self.lblResultSfrref)
        self.leResultSpcSfrref = QtWidgets.QLineEdit(wdg)
        self.leResultSpcSfrref.setFixedWidth(100) 
        layout_magsfrebmv.addWidget(self.leResultSpcSfrref)
        self.leResultSpcSfrref.setEnabled(False)  
        #Add the result ebmv
        self.lblResultEbmvref = QtWidgets.QLabel(' ebmv =', wdg)
        layout_magsfrebmv.addWidget(self.lblResultEbmvref)
        self.leResultSpcEbmvref = QtWidgets.QLineEdit(wdg)
        self.leResultSpcEbmvref.setFixedWidth(100) 
        layout_magsfrebmv.addWidget(self.leResultSpcEbmvref)
        self.leResultSpcEbmvref.setEnabled(False)  
        #add empty label to fill the hbox
        lbl = QtWidgets.QLabel(' ', wdg)
        layout_magsfrebmv.addWidget(lbl)
        layout.addLayout(layout_magsfrebmv, layoutRow, 1, 1, 1)

        ######### SHOW section ###############################################    
        layoutRow += 1    
        #Add the show button
        layoutRow += 1
        self.btn = QtWidgets.QPushButton('Show', wdg)
        self.btn.setFixedWidth(500)
        self.btn.setFixedHeight(50)
        self.btn.clicked.connect(self.bt_showAViewWidget)
        self.btn.setToolTip('Display the Chi2/spectrum/fitted template/linemodel results successively...')
        layout.addWidget(self.btn, layoutRow, 1, 1, 1)
        self.btnReprocess = QtWidgets.QPushButton('Reprocess', wdg)
        self.btnReprocess.setFixedWidth(100)
        self.btnReprocess.setFixedHeight(50)
        self.btnReprocess.clicked.connect(self.bt_reprocess)
        self.btnReprocess.setToolTip('Reprocess this source...')
        layout.addWidget(self.btnReprocess, layoutRow, 2, 1, 1)

        #Add the extremum choice ctrls
        layoutRow += 1
        layoutRow += 1
        self.lblExtremumChoice = QtWidgets.QLabel(' Extremum: ', wdg)
        layout.addWidget(self.lblExtremumChoice, layoutRow, 0, 1, 1)
        
        self.leExtremumChoice = QtWidgets.QLineEdit(wdg)
        self.leExtremumChoice.setToolTip('Enter the extremum num to be selected for display...')
        self.leExtremumChoice.setFixedWidth(100) 
        self.leExtremumChoice.setEnabled(False)
        layout.addWidget(self.leExtremumChoice, layoutRow, 1, 1, 10)
        
        self.ckExtremumChoiceOverride = QtWidgets.QCheckBox('Extremum Override', wdg)
        self.ckExtremumChoiceOverride.stateChanged.connect(self.ck_extremumOverride)
        layout.addWidget(self.ckExtremumChoiceOverride, layoutRow, 2, 1, 10)
        #self.ckExtremumChoiceOverride.toggle()
        
        #self.setCentralWidget(wdg)
        self.setLayout(layout)        
        self.setWindowTitle('Aview')
        

        #Set path from settings
        self.settings = QtCore.QSettings("amazed", "aviewgui")
        _resDir = self.settings.value("resDir", defaultValue = "", type=str)
        self.leResDir.setText(_resDir)      
        if _resDir=="":
            self.enableCtrls(False)
            
        #Set spc filter
        self.leSpcFilter.setText("")
        #self.leSpcFilter.setEnabled(False)
            
        #Set zdiff Threshold
        self.leDiffThres.setText("")
            

        #auto load from inputs args
        if not init_resdir=="" and not init_refpath=="":
            self.leResDir.setText(init_resdir) 
            self.loadresults(init_refpath)
            
        wdg.setStyleSheet("*[coloredcell=\"true\"] {background-color:rgb(215,215,215);}")
        self.show()
 
    def bt_loadresults(self):
        self.loadresults()

    def loadresults(self, init_refPath=""):
        self.setCurrentDir()
    
        _resDir = str(self.leResDir.text())
        print("using _resDir: {}".format(_resDir))
        
        statsPath = os.path.join(_resDir, "stats")
        diffPath = os.path.join(statsPath, "diff.txt")
        if not os.path.exists(diffPath):
            print("diff file not found... computing the diff file is necessary: please select the reference redshift file list in order continue !")
            calcFile = os.path.join(_resDir, "redshift.csv") 
            
            if not init_refPath=="":
                print("-using init_refPath={}".format(init_refPath))
                _refzfilepath = init_refPath
            else:
                print("-browse path: select ref file path...".format())
                _refzfilepathDefault = self.settings.value("refzfilepath", defaultValue = "", type=str)
                _refzfilepath = str(QtWidgets.QFileDialog.getOpenFileName(self, "Select Reference Redshift list file", _refzfilepathDefault))
                
            if os.path.exists(_refzfilepath) :
                print("Selected zref file found... processing with pfs type !")
                
                if os.path.isdir(statsPath)==False:
                    os.mkdir( statsPath, 0o755 );
                self.settings.setValue("refzfilepath", _refzfilepath)
                processAmazedOutputStats.setPFSRefFileType()
                processAmazedOutputStats.ProcessDiff( _refzfilepath, calcFile, diffPath , reftype='pfs')
            else:
                print("Selected zref file does not exist... processing without refFile !")
                
                if os.path.isdir(statsPath)==False:
                    os.mkdir( statsPath, 0o755 );
                processAmazedOutputStats.setNoRefFileType()
                processAmazedOutputStats.ProcessDiff( _refzfilepath, calcFile, diffPath , reftype='no')
                
        
        _spcName = str(self.leSpcFilter.text())
        print("INFO: using spcName filter = {}".format(_spcName))
        _diffthres = str(self.leDiffThres.text())
        if _diffthres=="":
            _diffthres = -1
        
        #filter by zref
        try:
            _zrefmin = float(str(self.lezrefrangemin.text()))
            if not (_zrefmin >0.0 and _zrefmin < 20.0):
                _zrefmin = -1.0
        except:
            print("ERROR with zref min value, please check your input ! (using zrefmin=-1)")
            _zrefmin = -1
        try:
            _zrefmax = float(str(self.lezrefrangemax.text()))
            if not (_zrefmax >0.0 and _zrefmax < 20.0):
                _zrefmax = 20.0
        except:
            print("ERROR with zref max value, please check your input ! (using zrefmax=20)")
            _zrefmax = 20.0   
            
        #filter by zref
        try:
            _sfrrefmin = float(str(self.lesfrrangemin.text()))
            if not (_sfrrefmin >0.0 and _sfrrefmin < 1000.0):
                _sfrrefmin = -1.0
        except:
            print("ERROR with sfr min value, please check your input ! (using sfrrefmin=-1)")
            _sfrrefmin = -1
        try:
            _sfrrefmax = float(str(self.lesfrrangemax.text()))
            if not (_sfrrefmax >0.0 and _sfrrefmax < 20.0):
                _sfrrefmax = 20.0
        except:
            print("ERROR with sfr max value, please check your input ! (using sfrrefmax=20)")
            _sfrrefmax = 20.0
            
        self.resList = resultstat.ResultList(_resDir, spcName=_spcName, diffthreshold=float(_diffthres), opt='brief',  
                                             zrefmin=_zrefmin, zrefmax=_zrefmax, 
                                             sfrrefmin=_sfrrefmin, sfrrefmax=_sfrrefmax) 
        self.leResultsN.setText(str(self.resList.n))
        #init result index to 1 (first index)
        self.leResultIndex.setText(str(1))
        self.refreshResultDetails()

    def refreshResultDetails(self, reinit=False):
        if not reinit:
            current_index = int(self.leResultIndex.text() )-1
            self.leResultSpcName.setText(self.resList.list[current_index].name)
            self.leResultSpcZdiff.setText(str(self.resList.list[current_index].zdiff))
            self.leResultSpcZcalc.setText("{:.5}".format(self.resList.list[current_index].zcalc))
            self.leResultSpcZref.setText("{:.5}".format(self.resList.list[current_index].zref))
            self.leResultSpcMagref.setText("{:.5}".format(self.resList.list[current_index].refValues['mag']))
            self.leResultSpcSfrref.setText("{:.5}".format(self.resList.list[current_index].refValues['sfr']))
            self.leResultSpcEbmvref.setText("{:.5}".format(self.resList.list[current_index].refValues['ebmv']))
        else:
            self.leResultSpcName.setText("")
            self.leResultSpcZdiff.setText("")
            self.leResultSpcZcalc.setText("")
            self.leResultSpcZref.setText("")
            self.leResultSpcMagref.setText("")
            self.leResultSpcSfrref.setText("")
            self.leResultSpcEbmvref.setText("")
    
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
        
    def leResultIndex_handleEditingFinished(self):
        if self.leResultIndex.isModified():
            iManual = float(self.leResultIndex.text())
            
            self.refreshResultDetails()
            print('Editing Finished: using imanual = {}'.format(iManual))
        self.leResultIndex.setModified(False)
        
        
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
            
            _refzfilepathDefault = self.settings.value("refzfilepath", defaultValue = "", type=str)
            _refzfilepath = str(QtWidgets.QFileDialog.getOpenFileName(self, "Select Reference Redshift list file", _refzfilepathDefault))
            if os.path.exists(_refzfilepath) :
                if os.path.isdir(statsPath)==False:
                    os.mkdir( statsPath, 0o755 );
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

    #not used, this just calls the aview command line interface
    def bt_showResult(self):
        _resDir = str(self.leResDir.text())
        
        _spcName = str(self.leResultSpcName.text())
        _tplpath = ""
        
        #current_index = int(self.leResultIndex.text() )-1
        _redshift = ""#self.resList.list[current_index].zcalc

                
        if self.ckExtremumChoiceOverride.isChecked():
            iextrFromCtrl = str(self.leExtremumChoice.text())
            if iextrFromCtrl=="":
                _iextremaredshift = 0
            else:
                _iextremaredshift = float(iextrFromCtrl)
        else:
            _iextremaredshift = 0
            
        _diffthres = str(self.leDiffThres.text())
        if _diffthres=="":
            _diffthres = -1
        _failureindex = "0"
        aview.plotRes(_resDir, _spcName, _tplpath, _redshift, _iextremaredshift, _diffthres, _failureindex)
        
        
    def bt_showAViewWidget(self):
        tag = "bt_showAViewWidget"
        _resDir = str(self.leResDir.text())
        
        _spcName = str(self.leResultSpcName.text())
        _spcIdx = -1
        for idx in range(self.resList.n):
            if _spcName == self.resList.list[idx].name:
                print("INFO: found index for this spectrum and diffthreshold : i={}\n".format(idx))
                _spcIdx = idx
                break
            
        if self.ckExtremumChoiceOverride.isChecked():
            iextrFromCtrl = str(self.leExtremumChoice.text())
            if iextrFromCtrl=="":
                _iextremaredshift = 0
            else:
                _iextremaredshift = float(iextrFromCtrl)
        else:
            _iextremaredshift = 0
            
        _resParser = resparser.ResParser(_resDir)

        #self.AViewWidget = aviewwidget.AViewWidget(parent=self, resParser=_resParser, resList=self.resList, resIdx=_spcIdx, iextremaredshift=_iextremaredshift)
        try:
            self.AViewWidget = aviewwidget.AViewWidget(parent=self, resParser=_resParser, resList=self.resList, resIdx=_spcIdx, iextremaredshift=_iextremaredshift)
            print("AviewWidget created.")            
            self.AViewWidget.show()
        except Exception as e:
            if 1:
                print( e )
                traceback.print_exc()
            print("\n".format())
            print("ERROR: Unable to show this result... (i={})".format(_spcIdx))
            print("NB: Maybe this source hasn't been processed successfully by amazed...(cf. zcalc=-1)".format())
            print("NB: Maybe the intermediate results are missing for this spectrum...".format())
            print("NB: Could be another problem too...".format())
            
    def bt_reprocess(self):
        self.aviewwidget_reprocess()
        
    def bt_setResultDir(self):
        self.refreshResultDetails(reinit=True) 
        self.setResultDir()
    
    def setResultDir(self, newPath=None):
        if newPath==None:
            _resDirDefault = os.path.abspath(str(self.leResDir.text()))
            #_resDirDefault = _resDirDefault[:_resDirDefault.index(os.sep)] if os.sep in _resDirDefault else _resDirDefault
            
            _resDir = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory",_resDirDefault))
        else:
            _resDir = newPath
            
        if os.path.exists(_resDir):
            #check the diff file is present, if not, do something...
            print("_resDir = {}".format(_resDir))
            print("previous from settings = {}".format(self.settings.value("resDir", defaultValue = "", type=str)))
            
            if not _resDir == self.settings.value("resDir", defaultValue = "", type=str):
                _resDirPrevious = self.settings.value("resDir", defaultValue = "", type=str)
                self.settings.setValue("resDirPrevious", _resDirPrevious);
            self.leResDir.setText(_resDir)
            self.settings.setValue("resDir", _resDir);
            self.enableCtrls(True)
        
    def bt_setResultDirPrevious(self):
        _resDirPrevious = self.settings.value("resDirPrevious", defaultValue = "", type=str)
        self.leResDir.setText(_resDirPrevious)
        self.settings.setValue("resDir", _resDirPrevious);
        self.enableCtrls(True)
        
    def bt_setResultDirFavorite(self):
        _path = str(self.leResDir.text())
        if os.path.exists(_path):
            self.settings.setValue("resDirFavorite", _path)
    
    def bt_loadResultDirFavorite(self):
        _path = self.settings.value("resDirFavorite", defaultValue = "", type=str)
        self.setResultDir(_path)
        
    def ck_extremumOverride(self, state):
        if state == QtCore.Qt.Checked:
            self.leExtremumChoice.setEnabled(True)
            iextrFromCtrl = str(self.leExtremumChoice.text())
            if iextrFromCtrl=="":
                self.leExtremumChoice.setText("0")
        else:
            self.leExtremumChoice.setEnabled(False)
            
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
    
    def aviewwidget_nextResult(self):
        self.AViewWidget.close()
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        self.bt_nextResultIndex()
        self.repaint()
        time.sleep(1.0)
        self.bt_showAViewWidget()
        QtWidgets.QApplication.restoreOverrideCursor()
        
        
    def aviewwidget_reprocess(self):
        try:
            self.AViewWidget.close()
        except:
            pass
        
        _resDir = str(self.leResDir.text())
        _spcName = str(self.leResultSpcName.text())
        _initReprocessData = [_resDir, _spcName]
        print("AVIEWGUI: reprocess with : _resDir={}".format(_resDir))
        print("AVIEWGUI: reprocess with : _spcName={}".format(_spcName))
        
        self.aprocessWindow = aprocessgui.AProcessGui(obj=None, initReprocess=_initReprocessData)
        self.aprocessWindow.show()
    
 
def main():
    app = QtWidgets.QApplication(sys.argv)
    ins = ObjectAViewGui()
    ex = AViewGui()
    sys.exit(app.exec_())
 
if __name__ == '__main__':
    main()