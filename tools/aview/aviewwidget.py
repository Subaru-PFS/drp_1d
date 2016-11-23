# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 12:26:11 2016

@author: aschmitt

This handles the display of amazed results 
"""

 
import sys
import os
import inspect
import glob
import time
import threading
import shutil

import numpy as np


#import mpl only to force Qt5Agg
import matplotlib as mpl
mpl.use('Qt5Agg')

#from matplotlib.backends import qt_compat
from PyQt5 import QtGui, QtCore, QtWidgets

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas    
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

#print(plt.style.available)

plt.style.use('ggplot')
#plt.style.use('seaborn-white')
#plt.style.use('seaborn-darkgrid')
#plt.style.use('seaborn-dark-palette')
#[u'seaborn-darkgrid', u'seaborn-notebook', u'classic', u'seaborn-ticks', u'grayscale', u'bmh', 
#u'seaborn-talk', u'dark_background', u'ggplot', u'fivethirtyeight', u'seaborn-colorblind', 
#u'seaborn-deep', u'seaborn-whitegrid', u'seaborn-bright', u'seaborn-poster', u'seaborn-muted', 
#u'seaborn-paper', u'seaborn-white', u'seaborn-pastel', u'seaborn-dark', u'seaborn-dark-palette']

#plt.rcParams['axes.facecolor']='w'
#plt.style.use('seaborn-white')
#import seaborn as sns
#sns.set_context("poster")
#sns.set_style("whitegrid")

import resparser
import aviewplot
import chisquare as chisq

class Graph(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, fig=None, graphType="merit"):
        if not fig == None:
            self.fig = fig
        else:
            self.fig = Figure()
        FigureCanvas.__init__(self, self.fig)
        self.ax = self.fig.add_subplot(111)
        self.parent = parent
        self.type = graphType
        
    def mouseDoubleClickEvent(self, event):
        if self.type=='merit':
            inv = self.ax.transData.inverted()
            _x, _y = inv.transform((event.x(), self.frameSize().height() - event.y()))
            print("Double Click x={} - y={}".format(_x, _y))
            self.parent.plotCandidate(_x)


class candidateListWidget(QtWidgets.QListWidget):
    def __init__(self, parent=None):
        self.parent = parent
        QtWidgets.QListWidget.__init__(self)
        
    def Clicked(self, item):
        #k = item.data(QtCore.Qt.UserRole).toPyObject()
        k = item.data(QtCore.Qt.UserRole)
        print("candidateListWidget: you clicked = {}".format(k))        
        self.parent.plotCandidate(zClick=-1, iextremaredshift=k)
      
class ObjectAViewWidget(object):
    def __init__(self, parent=None):
        self.value = []
        self.parent = parent
        
class AViewWidget(QtWidgets.QWidget):
    def __init__(self, parent=None, resParser=None, resList=[], resIdx=-1, iextremaredshift=-1):
        super(AViewWidget, self).__init__()
        self.resParser = resParser
        self.resIdx = resIdx
        self.resList = resList
        self.iextremaredshift = iextremaredshift
        self.initCandidatesList()
        
        self.parent = parent
        wdg = self#QtGui.QWidget()
        
        pal=QtGui.QPalette()
        role = QtGui.QPalette.Background
        pal.setColor(role, QtGui.QColor(255, 255, 255))
        self.setPalette(pal)
      
        layout = QtWidgets.QGridLayout(wdg) 
        layoutRow = 0
        
        #Add the redshift manual override
        layoutRedshift = QtWidgets.QHBoxLayout()
        self.lblRedshiftChoice = QtWidgets.QLabel('Redshift: ', wdg)
        layoutRedshift.addWidget(self.lblRedshiftChoice)
        
        self.leRedshiftChoice = QtWidgets.QLineEdit(wdg)
        self.leRedshiftChoice.setToolTip('Enter the manual redshift to be used for display...')
        self.leRedshiftChoice.editingFinished.connect(self.le_handleEditingFinished)
        #self.leRedshiftChoice.setFixedHeight(30) 
        self.leRedshiftChoice.setEnabled(False)
        layoutRedshift.addWidget(self.leRedshiftChoice)
        
        self.ckRedshiftChoiceOverride = QtWidgets.QCheckBox('Redshift Override', wdg)        
        #self.ckRedshiftChoiceOverride.setFixedHeight(30) 
        self.ckRedshiftChoiceOverride.stateChanged.connect(self.ck_redshiftOverride)
        layoutRedshift.addWidget(self.ckRedshiftChoiceOverride)
        
        
        #Add the navigation shortcuts
        layoutResNavigation = QtWidgets.QHBoxLayout()
        self.lblNavigate = QtWidgets.QLabel('Navigate: ', wdg)
        layoutResNavigation.addWidget(self.lblNavigate)
        self.btnNavReprocess = QtWidgets.QPushButton(' Go to re-process ', wdg)
        self.btnNavReprocess.clicked.connect(self.bt_setNavReprocess)
        layoutResNavigation.addWidget(self.btnNavReprocess)
        self.btnNavNext = QtWidgets.QPushButton(' Next ', wdg)
        self.btnNavNext.clicked.connect(self.bt_setNavNextResult)
        layoutResNavigation.addWidget(self.btnNavNext)
        

        #set the layout    
        layoutLeftColumn = QtWidgets.QVBoxLayout()
        #canddidates list
        self.lst_Candidates = candidateListWidget(self)
        self.lst_Candidates.setWindowTitle('Redshift Candidates')
        self.lst_Candidates.setToolTip("Left-click on a candidate (row) to display")
        self.lst_Candidates.itemClicked.connect(self.lst_Candidates.Clicked)
   
        layoutLeftColumn.addWidget(self.lst_Candidates) 
        self.populateCandidatesList()
        #chi2plot
        self.layoutChi2Plot = QtWidgets.QVBoxLayout()
        #add the chi2 operator merit plot combo selection
        self.cbChi2Selection = QtWidgets.QComboBox()
        self.cbChi2Selection.currentIndexChanged.connect(self.cb_chi2Selection_selectionchange)
        self.cbChi2SelectedIndex = -1
        self.layoutChi2Plot.addWidget(self.cbChi2Selection)  
        #add the tpl selcetion override for the chi2plot        
        self.cbChi2ModelSelection = QtWidgets.QComboBox()
        self.cbChi2ModelSelection.currentIndexChanged.connect(self.cb_chi2ModelSelection_selectionchange)
        self.cbChi2ModelSelectedIndex = -1
        self.layoutChi2Plot.addWidget(self.cbChi2ModelSelection)  
        self.populateChi2ModelCombo()
        if '.chisquare' in self.cbChi2Selection.currentText():
            self.cbChi2ModelSelection.setEnabled(True)
        else:
            self.cbChi2ModelSelection.setEnabled(False)
        
        #add chisquare.py figure
        self.loadMeritPlot(self.iextremaredshift)
        #self.layoutChi2Plot.addWidget(self.toolbarChi2)
        #self.layoutChi2Plot.addWidget(self.canvasChi2)      
        layoutLeftColumn.addLayout(self.layoutChi2Plot)
        
        #adding separators
        lbl = QtWidgets.QLabel('', wdg)
        lbl.setProperty("coloredcell", True)
        lbl.setFixedHeight(5)
        layoutLeftColumn.addWidget(lbl)         
        #lbl = QtWidgets.QLabel('', wdg)
        #layoutLeftColumn.addWidget(lbl) 

        layoutLeftColumn.addLayout(layoutRedshift)
        #add seperator
        lbl = QtWidgets.QLabel('', wdg)
        lbl.setProperty("coloredcell", True)
        lbl.setFixedHeight(5)
        layoutLeftColumn.addWidget(lbl) 
        
        layoutLeftColumn.addLayout(layoutResNavigation)
        
        v_widget = QtWidgets.QWidget()
        v_widget.setLayout(layoutLeftColumn)
        v_widget.setFixedWidth(500)
        
        layout.addWidget(v_widget, layoutRow, 1, 1, 1)
        
        lbl = QtWidgets.QLabel('', wdg)
        lbl.setFixedWidth(5)
        lbl.setProperty("coloredcell", True)
        layout.addWidget(lbl, layoutRow, 2, 1, 1) 

        # set the layout
        self.layoutAviewPlot = QtWidgets.QVBoxLayout()
        self.cbAviewPlotSelection = QtWidgets.QComboBox()
        self.cbAviewPlotSelection.currentIndexChanged.connect(self.cb_aviewplotselection_selectionchange)
        for a in self.candid_displayParamsBundle:
            op_name = a['operator']
            self.cbAviewPlotSelection.addItem(op_name)
            
        self.layoutAviewPlot.addWidget(self.cbAviewPlotSelection)
        self.loadAViewPlot(self.layoutAviewPlot)
        #self.layoutAviewPlot.addWidget(self.canvasAView)
        layout.addLayout(self.layoutAviewPlot, layoutRow, 3, 1, 1)        
        
        #self.setCentralWidget(wdg)
        #wdg.setLayout(layoutAviewPlot)     
        #self.setLayout(layoutAviewPlot)        
        self.setWindowTitle('AViewWidget')   

        screen = QtWidgets.QDesktopWidget().screenGeometry()
        wdg.setGeometry(screen.width()*0.2, screen.height()*0.1, screen.width()*0.75, screen.height()*0.75) 
                                               
        wdg.setStyleSheet("*[coloredcell=\"true\"] {background-color:rgb(215,215,215);}")
        self.show()
        
    def cb_aviewplotselection_selectionchange(self, i):
        tag = "cb_aviewplotselection_selectionchange"
        print("{}, Items in the list are :".format(tag))
        for count in range(self.cbAviewPlotSelection.count()):
            print("    {}".format(self.cbAviewPlotSelection.itemText(count)))
        print("{}: Current index = {}, selection changed to {}".format(tag, i, self.cbAviewPlotSelection.currentText()))
        if not i==self.candid_displayParamsBundleIdx:
            self.candid_displayParamsBundleIdx = i
            self.plotCandidate()
            
    def cb_chi2Selection_selectionchange(self,i):
        tag = "cb_chi2Selection_selectionchange"
        print("{}, Items in the list are :".format(tag))
        for count in range(self.cbChi2Selection.count()):
            print("    {}".format(self.cbChi2Selection.itemText(count)))
        print("{}: Current index = {}, selection changed to {}".format(tag, i, self.cbChi2Selection.currentText()))
        if not i==self.cbChi2SelectedIndex and not self.cbChi2SelectedIndex==-1:
            #if the selection corresponds to a chi2:
            if '.chisquare' in self.cbChi2Selection.currentText():
                self.cbChi2ModelSelection.setEnabled(True)
            else:
                self.cbChi2ModelSelection.setEnabled(False)
            self.cbChi2SelectedIndex = i
            self.plotCandidate()
            
    def cb_chi2ModelSelection_selectionchange(self,i):
        tag = "cb_chi2modelSelection_selectionchange"
        print("{}, Items in the list are :".format(tag))
        for count in range(self.cbChi2ModelSelection.count()):
            print("    {}".format(self.cbChi2ModelSelection.itemText(count)))
        print("{}: Current index = {}, selection changed to {}".format(tag, i, self.cbChi2ModelSelection.currentText()))
        if not i==self.cbChi2ModelSelectedIndex and not self.cbChi2ModelSelectedIndex==-1:
            self.cbChi2ModelSelectedIndex = i
            self.plotCandidate( zClick=-1, iextremaredshift=-1, tplNameOverride=self.cbChi2ModelSelection.currentText())    
            
    def cb_chi2ModelSelection_setselection(self, tplname):
        tag = "cb_chi2ModelSelection_setselection"
        #print("{}, Items in the list are :".format(tag))
        for count, name in enumerate(range(self.cbChi2ModelSelection.count())):
            if name==tplname:
                self.cbChi2ModelSelectedIndex = count
                self.cbChi2ModelSelection.SetSelection()
            print("    {}".format(self.cbChi2ModelSelection.itemText(count)))
        print("{}: Current index = {}, selection set to {}".format(tag, i, self.cbChi2ModelSelection.currentText()))
        
        
    def ck_redshiftOverride(self, state):
        if state == QtCore.Qt.Checked:
            self.leRedshiftChoice.setEnabled(True)
            zextrFromCtrl = str(self.leRedshiftChoice.text())
            if zextrFromCtrl=="":
                self.leRedshiftChoice.setText("0")
        else:
            self.leRedshiftChoice.setEnabled(False)
            
    def le_handleEditingFinished(self):
        if self.leRedshiftChoice.isModified():
            zManual = float(self.leRedshiftChoice.text())
            
            self.loadAViewPlot(self.layoutAviewPlot, zManual)
            print 'Editing Finished: using zmanual = {}'.format(zManual)
        self.leRedshiftChoice.setModified(False)
        
    def getClosestCandidateIdx(self, zClick):
        zArray = np.array(self.candid_zvalCandidates)
        zArrayDiff = np.abs(zArray-zClick)
        ind = np.argmin(zArrayDiff)
        
        return ind
        
    def plotCandidate(self, zClick=-1, iextremaredshift=-1, tplNameOverride=""):
        QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        
        if self.ckRedshiftChoiceOverride.isChecked():
            self.ckRedshiftChoiceOverride.toggle()        

        if iextremaredshift==-1 and not zClick==-1:
            #get the closest candidate for the clicked redshift
            self.iextremaredshift = self.getClosestCandidateIdx(zClick)
            #self.iextremaredshift += 1
        
        if not iextremaredshift==-1 and zClick==-1:
            #get the closest candidate for the clicked redshift
            self.iextremaredshift = iextremaredshift               
          
        #refresh meritplot
        self.loadMeritPlot(self.iextremaredshift, tplNameOverride=tplNameOverride)               
                
        #refresh aviewplot
        #self.layoutAviewPlot.removeWidget(self.canvasAView)
        #self.canvasAView.setParent(None)
        #remove all widgtes... 
        #for i in reversed(range(self.layoutAviewPlot.count())): 
        #    self.layoutAviewPlot.itemAt(i).widget().deleteLater()
        #self.show()
        self.loadAViewPlot(self.layoutAviewPlot)


    def populateChi2ModelCombo(self):
        self.cbChi2ModelSelection.clear()

        tplPathList = self.resParser.getTplFullPathList(component='raw')
            
        for k, tplPath in enumerate(tplPathList):
            tplname = os.path.split(tplPath)[1]
            self.cbChi2ModelSelection.addItem(tplname)
        self.cbChi2ModelSelectedIndex=0

    def populateCandidatesList(self):
        for k, zcand in enumerate(self.candid_zvalCandidates):
            tplPath = self.candid_displayParamsBundle[self.candid_displayParamsBundleIdx]['tplPaths'][k]
            tplname = os.path.split(tplPath)[1]
            item = QtWidgets.QListWidgetItem("C{:<8}{:<15}{}".format(k, zcand, tplname))
            data = k
            item.setData(QtCore.Qt.UserRole, data)

            self.lst_Candidates.addItem(item);

    def initCandidatesList(self):
        tag = "initCandidatesList"
        # 
        self.candid_zvalCandidates = []
        self.candid_displayParamsBundle = []        
        self.candid_displayParamsBundleIdx = 0
    
        spcnametag = self.resList.list[self.resIdx].name
        self.candid_zvalCandidates, self.candid_displayParamsBundle = self.resParser.getAutoCandidatesList(spcnametag)
        print("{}: {} candidates loaded".format(tag, len(self.candid_zvalCandidates)))
                
    def loadMeritPlot(self, idxCandidate=-1, tplNameOverride=""):
        tag = "loadMeritPlot"
        print("{}: idxCandidate={}".format(tag, idxCandidate))
        #reset the merit plot canvas if already existing
        try:
            self.canvasChi2.ax.cla()
            self.layoutChi2Plot.removeWidget(self.canvasChi2)
            self.canvasChi2.setParent(None)
        except:
            pass
        
        #retrieve the chisquare data PATH
        _sname = self.resList.list[self.resIdx].name
        if idxCandidate==-1:
            tplnametag = ""
        else:
            #find the first 'chi2' displayParamsBundleIdx in order to extract the template name: 
            displayParamsBundleIdx_chi2 = 0
            for a in range(len(self.candid_displayParamsBundle)):
                if 'chi2' in self.candid_displayParamsBundle[a]['operator']:
                    displayParamsBundleIdx_chi2 = a
                    break
            
            tplPath = self.candid_displayParamsBundle[displayParamsBundleIdx_chi2]['tplPaths'][idxCandidate]
            tplnametag = os.path.split(tplPath)[1]
            if not tplNameOverride=="":
                tplnametag = tplNameOverride
                print("INFO: overriding the tpl for merit plot with tpl named : {}".format(tplnametag))
        print("looking for the chi2 path using the tplnameteg = {}".format(tplnametag))
        [chipathlist, chinamelist] = self.resParser.getAutoChi2FullPath(_sname, tplnametag)

        #select the merit curve w. regard to the combobox
        if self.cbChi2SelectedIndex==-1:      
            self.cbChi2Selection.clear()
            for name in chinamelist:
                self.cbChi2Selection.addItem(name)
            self.cbChi2SelectedIndex=0
            
        #actually load the chi2 data
        print("using cbChi2SelectedIndex={}".format(self.cbChi2SelectedIndex))
        chi2path = chipathlist[self.cbChi2SelectedIndex]
        chi2 = chisq.ResultChisquare(chi2path)
        
        #plot the merit curve in the canvas
        try:
            self.figureChi2 = chi2.plot(showContinuumEstimate=False, showExtrema=True, showAmbiguities=False, enablePlot=False, exportPath="", enableReturnFig=True)       
        except:
            print("unable to load the merit curve plot...")
            self.figureChi2 = plt.plot()
        self.canvasChi2 = Graph(self, self.figureChi2, graphType="merit")
        #self.canvasChi2.setFixedHeight(550.0)
        #self.canvasChi2.setFixedWidth(500.0)
        self.canvasChi2.setToolTip("Double left-click on the extrema (red circles) to display")
        
        #refresh the toolbar canvas
        try:
            self.toolbarChi2
            self.layoutChi2Plot.removeWidget(self.toolbarChi2)
            self.toolbarChi2.setParent(None)
        except AttributeError:
            pass
        self.toolbarChi2 = NavigationToolbar(self.canvasChi2, self)
        
        self.layoutChi2Plot.addWidget(self.toolbarChi2)
        self.layoutChi2Plot.addWidget(self.canvasChi2) 
        
        
    def loadAViewPlot(self, layout, zManual=-1):
        tag = "loadAViewPlot"
        
        #reset the aview plot canvas if already existing
        try:
            self.layoutAviewPlot.removeWidget(self.canvasAView)
            self.canvasAView.setParent(None)
        except:
            pass
        
        
        _sname = self.resList.list[self.resIdx].name
        spath = self.resParser.getSpcFullPath(_sname)
        try:
            npath = self.resParser.getNoiseFullPath(_sname)
        except:
            npath = ""
            
        if not self.iextremaredshift == -1: 
            idxExtrema = int(self.iextremaredshift)
        else:
            idxExtrema = int(0)
        print('{}: retrieve tpl-model path - using idxExtrema: {}'.format(tag, idxExtrema))
            
        # todo, when idxExtrema is -1, z should be found from the "zval = s.getRedshiftVal(spcName)"
        zvalCandidate = self.candid_zvalCandidates[idxExtrema]
        tplpathCandidate =  self.candid_displayParamsBundle[self.candid_displayParamsBundleIdx]['tplPaths'][idxExtrema]

        forceTplAmplitudeCandidate = self.candid_displayParamsBundle[self.candid_displayParamsBundleIdx]['forceTplAmplitudes'][idxExtrema]
        forceTplDustCoeffCandidate = self.candid_displayParamsBundle[self.candid_displayParamsBundleIdx]['forceTplDustCoeff'][idxExtrema]
        
        forceTplDoNotRedShiftCandidate = self.candid_displayParamsBundle[self.candid_displayParamsBundleIdx]['forceTplDoNotRedShifts'][idxExtrema]
        print("{}: full tpl-model path: {}".format(tag, tplpathCandidate))
        if not os.path.exists(tplpathCandidate):
            print("{}: tplPath does not exist : reset tpl-model path = {}".format(tag, tplpathCandidate))
            tplpathCandidate = ""
            forceTplAmplitudeCandidate = -1
            forceTplDoNotRedShiftCandidate = 0        

        #tpath = "/home/aschmitt/data/vuds/VUDS_flag3_4/amazed/templates/ExtendedTemplatesMarch2016_v2/emission/NEW-Sbc-extended.txt"
        #tpath = "/home/aschmitt/tmp/output/sc_530002397_F53P002_join_A_125_1_atm_clean/linemodelsolve.linemodel_spc_extrema_0.csv"
        
        cpath = self.resParser.getCatalogFullPath()
        if not os.path.exists(cpath):
            print("{}: ERROR: linecatalog file not found, please check : {}".format(tag, cpath))
            
        if zManual==-1:
            z = zvalCandidate
        else:
            z = zManual
        
        spectrumContinuumPath = self.resParser.getContinuumPath(_sname, methodForced='auto')
        print("using continuum path = {}".format(spectrumContinuumPath))
        #print("Z for aviewplot = {}".format(z))
        avp = aviewplot.AViewPlot(spath, 
                                               npath, 
                                               tplpathCandidate, 
                                               cpath, 
                                               z, 
                                               forceTplAmplitude=forceTplAmplitudeCandidate, 
                                               forceTplDustCoeff=forceTplDustCoeffCandidate,
                                               forceTplDoNotRedShift=forceTplDoNotRedShiftCandidate, 
                                               enablePlot=False, 
                                               scontinuumpath = spectrumContinuumPath) 
                                               
        figureAView = avp.plot(enableReturnFig=True)
        self.canvasAView = FigureCanvas(figureAView)  
        try:
            self.toolbarAView
            self.layoutAviewPlot.removeWidget(self.toolbarAView)
            self.toolbarAView.setParent(None)
        except AttributeError:
            pass
        self.toolbarAView = NavigationToolbar(self.canvasAView, self)
        self.layoutAviewPlot.addWidget(self.toolbarAView)
        layout.addWidget(self.canvasAView)  
        
        #set the manual redshift value
        self.leRedshiftChoice.setText(str(z))
        
        #set the list candidates current row
        self.lst_Candidates.setCurrentRow( self.iextremaredshift );
        
        QtWidgets.QApplication.restoreOverrideCursor()
        print('{}: finished\n\n'.format(tag))
    
        
    def bt_setNavNextResult(self):
        #self.close()
        self.parent.aviewwidget_nextResult()
        
    def bt_setNavReprocess(self):
        self.parent.aviewwidget_reprocess()
        
    def closeEvent(self, event):
        print "AviewWidget: Closing"
        self.canvasChi2.ax.cla()
        #super(MainWidget, self).closeEvent(event)
        
def main():
    print("\nAViewWidget\n")
    app = QtWidgets.QApplication(sys.argv)
    ins = ObjectAViewWidget(parent=app)
    ex = AViewWidget(ins)
    
    
    #w = 1280; h = 720
#    w = 1920; h = 1200
#    frm = QtGui.QFrame ()
#    frm.sizeHint = lambda: QtCore.QSize (w, h)
#    ex.sizeHint = lambda: QtCore.QSize (w, h)

    sys.exit(app.exec_())
 
if __name__ == '__main__':
    main()
