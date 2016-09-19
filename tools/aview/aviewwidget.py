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

#from matplotlib.backends import qt_compat
from PyQt5 import QtGui, QtCore, QtWidgets
    
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas    
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

plt.style.use('ggplot')
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
        self.loadAViewPlot(self.layoutAviewPlot)
        #self.layoutAviewPlot.addWidget(self.canvasAView)
        layout.addLayout(self.layoutAviewPlot, layoutRow, 3, 1, 1)        
        
        #self.setCentralWidget(wdg)
        #wdg.setLayout(layoutAviewPlot)     
        #self.setLayout(layoutAviewPlot)        
        self.setWindowTitle('AViewWidget')   
                    
        wdg.setStyleSheet("*[coloredcell=\"true\"] {background-color:rgb(215,215,215);}")
        self.show()
        
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
        zArray = np.array(self.zvalCandidates)
        zArrayDiff = np.abs(zArray-zClick)
        ind = np.argmin(zArrayDiff)
        
        return ind
        
    def plotCandidate(self, zClick=-1, iextremaredshift=-1):
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
        self.loadMeritPlot(self.iextremaredshift)               
                
        #refresh aviewplot
        #self.layoutAviewPlot.removeWidget(self.canvasAView)
        #self.canvasAView.setParent(None)
        #remove all widgtes... 
        #for i in reversed(range(self.layoutAviewPlot.count())): 
        #    self.layoutAviewPlot.itemAt(i).widget().deleteLater()
        #self.show()
        self.loadAViewPlot(self.layoutAviewPlot)
        
    def populateCandidatesList(self):
        for k, zcand in enumerate(self.zvalCandidates):
            tplname = os.path.split(self.tplpathCandidates[k])[1]
            item = QtWidgets.QListWidgetItem("C{:<5}{:<12}{}".format(k, zcand, tplname))
            data = k
            item.setData(QtCore.Qt.UserRole, data)

            self.lst_Candidates.addItem(item);

    def initCandidatesList(self):
        tag = "initCandidatesList"
        # 
        self.zvalCandidates = []
        self.tplpathCandidates = []
        self.forceTplAmplitudeCandidates = []
        self.forceTplDoNotRedShiftCandidate = []
    
        _sname = self.resList.list[self.resIdx].name
        self.zvalCandidates, self.tplpathCandidates, self.forceTplAmplitudeCandidates, self.forceTplDoNotRedShiftCandidates = self.resParser.getAutoCandidatesList(_sname)
        print("{}: {} candidates loaded".format(tag, len(self.zvalCandidates)))
        
 
        
    def loadMeritPlot(self, idxCandidate=-1):
        tag = "loadMeritPlot"
                        
        try:
            self.canvasChi2.ax.cla()
            self.layoutChi2Plot.removeWidget(self.canvasChi2)
            self.canvasChi2.setParent(None)
        except:
            pass
        
        _sname = self.resList.list[self.resIdx].name
        #chi2path = "/home/aschmitt/tmp/output/sc_530002397_F53P002_join_A_125_1_atm_clean/linemodelsolve.linemodel.csv"
        if idxCandidate==-1:
            tplnametag = ""
        else:
            tplnametag = os.path.split(self.tplpathCandidates[idxCandidate])[1]
        [chipathlist, chinamelist] = self.resParser.getAutoChi2FullPath(_sname, tplnametag)
        chi2path = chipathlist[0]
        
        chi2 = chisq.ResultChisquare(chi2path)
        self.figureChi2 = chi2.plot(showContinuumEstimate=False, showExtrema=True, showAmbiguities=False, enablePlot=False, exportPath="", enableReturnFig=True)       
        self.canvasChi2 = Graph(self, self.figureChi2, graphType="merit")
        self.canvasChi2.setFixedHeight(475.0)
        self.canvasChi2.setFixedWidth(475.0)
        self.canvasChi2.setToolTip("Double left-click on the extrema (red circles) to display")
        
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
        
        try:
            self.layoutAviewPlot.removeWidget(self.canvasAView)
            self.canvasAView.setParent(None)
        except:
            pass
        
        #QObjectCleanupHandler().add(self.canvasAView())
        #self.layoutAviewPlot.addWidget(self.canvasAView)
        _sname = self.resList.list[self.resIdx].name
        #spath = "/home/aschmitt/data/vuds/VUDS_flag3_4/VUDS_flag3_4/ECDFS/spec1d_spec1dnoise/sc_530002397_F53P002_join_A_125_1_atm_clean.fits"
        spath = self.resParser.getSpcFullPath(_sname)
        #npath = "/home/aschmitt/data/vuds/VUDS_flag3_4/VUDS_flag3_4/ECDFS/spec1d_spec1dnoise/sc_530002397_F53P002_join_A_125_1_noise.fits"
        npath = self.resParser.getNoiseFullPath(_sname)
        

        if not self.iextremaredshift == -1: 
            idxExtrema = int(self.iextremaredshift)
        else:
            idxExtrema = int(0)
        print('{}: retrieve tpl-model path - using idxExtrema: {}'.format(tag, idxExtrema))
            
        # todo, when idxExtrema is -1, z should be found from the "zval = s.getRedshiftVal(spcName)"
        zvalCandidate = self.zvalCandidates[idxExtrema]
        tplpathCandidate = self.tplpathCandidates[idxExtrema]
        forceTplAmplitudeCandidate = self.forceTplAmplitudeCandidates[idxExtrema]
        forceTplDoNotRedShiftCandidate = self.forceTplDoNotRedShiftCandidates[idxExtrema]
        print("{}: full tpl-model path: {}".format(tag, tplpathCandidate))
        if not os.path.exists(tplpathCandidate):
            print("{}: tplPath does not exist : reset tpl-model path = {}".format(tag, tplpathCandidate))
            tplpathCandidate = ""
            forceTplAmplitudeCandidate = -1
            forceTplDoNotRedShiftCandidate = 0        

        #tpath = "/home/aschmitt/data/vuds/VUDS_flag3_4/amazed/templates/ExtendedTemplatesMarch2016_v2/emission/NEW-Sbc-extended.txt"
        #tpath = "/home/aschmitt/tmp/output/sc_530002397_F53P002_join_A_125_1_atm_clean/linemodelsolve.linemodel_spc_extrema_0.csv"
        tpath = tplpathCandidate
        
        cpath = self.resParser.getCatalogFullPath()
        if not os.path.exists(cpath):
            print("{}: ERROR: linecatalog file not found, please check : {}".format(tag, cpath))
            
        if zManual==-1:
            z = zvalCandidate
        else:
            z = zManual
        forceTplAmplitude = forceTplAmplitudeCandidate
        forceTplDoNotRedShift = forceTplDoNotRedShiftCandidate
        
        avp = aviewplot.AViewPlot(spath, 
                                               npath, 
                                               tpath, 
                                               cpath, 
                                               z, 
                                               forceTplAmplitude=forceTplAmplitude, 
                                               forceTplDoNotRedShift=forceTplDoNotRedShift, 
                                               enablePlot=False, 
                                               scontinuumpath = "") 
                                               
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
