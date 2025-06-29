﻿# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "13 Jun 2025"

import sys
import os
import re
import platform
from functools import partial
from collections import OrderedDict
import webbrowser
import numpy as np
import matplotlib as mpl
from matplotlib.figure import Figure
from pyparsing import ParseBaseException
sys.path.append(os.path.join('..'))  # analysis:ignore
from XAFSmass import XAFSmassCalc as xc
from XAFSmass.__init__ import (__version__, __author__, __license__)
from XAFSmass.materials_simple import elemental, compounds

# os.environ["QT_API"] = 'pyside6'
import qtpy
from qtpy import QtGui, QtCore
import qtpy.QtWidgets as QtWidgets

from matplotlib.backends import qt_compat
import matplotlib.backends.backend_qtagg as mpl_qt

Canvas = mpl_qt.FigureCanvasQTAgg
ToolBar = mpl_qt.NavigationToolbar2QT

# ==this can make the formulas look nicer:=====================================
mpl.rcParams['mathtext.fontset'] = 'custom'

# mpl.rcParams['mathtext.rm'] = 'serif'
# mpl.rcParams['mathtext.it'] = 'serif:italic'
# mpl.rcParams['mathtext.bf'] = 'serif:bold'
# mpl.rcParams['mathtext.cal'] = 'cursive'
# mpl.rcParams['mathtext.rm'] = 'serif'
# mpl.rcParams['mathtext.tt'] = 'monospace'

mpl.rcParams['mathtext.rm'] = 'cmr10'
mpl.rcParams['mathtext.it'] = 'cmmi10'
mpl.rcParams['mathtext.bf'] = 'cmb10'
mpl.rcParams['mathtext.cal'] = 'cmsy10'
mpl.rcParams['mathtext.bf'] = 'cmtt10'
mpl.rcParams['mathtext.bf'] = 'cmss10'
# =============================================================================

# Set the default color cycle:
# mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'k']

MAC = "qt_mac_set_native_menubar" in dir()

POWDER, FOIL, GAS, XCONTENT = range(4)
whats = 'powder', 'foil, film, glass etc.', 'gas', 'has unknown concentration'
formulas = (
    r"$\nu = (\mu_T d)\cdot S\cdot\left(N_A 2r_0 \lambda \sum_i{N_i f_i''}"
    r"\right)^{-1}; \quad m = M\cdot\nu$",
    r"$d = (\mu_T d)\cdot M\cdot\left(\rho N_A 2r_0 \lambda \sum_i{N_i f_i''}"
    r"\right)^{-1}$",
    r"$p = -\ln(1-{\rm att})\cdot kT \cdot"
    r"\left(d 2r_0 \lambda \sum_i{N_i f_i''}\right)^{-1}$",
    r"$\Delta\mu/\mu_T = N_x\Delta f_x'' \cdot"
    r"\left(\sum_{i\ne x}{N_if_i''} + N_xf_x''\right)^{-1}$")
examples = (
    "Cu(NO3)2  <i>or</i>  Cu%1Zn%1((Al2O3)%10SiO2)", "Cu%25Zn",
    "Ar <i>or</i> N2 <i>or</i> (Ar)(N2) to explore mixtures",
    "FexSiO2")
tables = ("Henke", "Brennan&Cowan", "Chantler (NIST)", "Chantler total (NIST)")
tablesF = ("Henke", "BrCo", "Chantler", "Chantler total")
edges = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "N1", "N2", "N3")

gasMixerN = 5
gasMixerPressureMin = 0
gasMixerPressureMax = 3000
gasMixerPressurePageStep = 50
gasMixerPressureDefault = 1000


class ComboBoxWithPlaceholder(QtWidgets.QComboBox):
    def paintEvent(self, event):
        painter = QtWidgets.QStylePainter(self)

        opt = QtWidgets.QStyleOptionComboBox()
        self.initStyleOption(opt)
        painter.drawComplexControl(QtWidgets.QStyle.CC_ComboBox, opt)

        if self.currentIndex() < 0:
            painter.setPen(QtGui.QColor('#888888'))
            try:
                if self.placeholderText():
                    opt.currentText = self.placeholderText()
            except AttributeError:
                pass
        painter.drawControl(QtWidgets.QStyle.CE_ComboBoxLabel, opt)


class MyFormulaMplCanvas(Canvas):
    def __init__(self, parent=None, width=5, height=0.4):
        fig = Figure(figsize=(width, height), dpi=96)
        self.fig = fig
        Canvas.__init__(self, fig)
        fig.patch.set_visible(False)
        self.setParent(parent)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Fixed)
        self.updateGeometry()
        fm = QtGui.QFontMetrics(self.font())
        self.fontsize = int(fm.height()) / 1.33
        locos = platform.platform(terse=True)
        if 'Linux' in locos:
            self.fontsize = int(fm.height()) / 1.5
        else:
            self.fontsize = int(fm.height()) / 1.33
        self.setStyleSheet("background-color:transparent;")

    def update_formula(self, formula=None):
        self.fig.clf()
        self.fig.suptitle(formula, x=0.5, y=0.48, ha='center', va='center',
                          fontsize=self.fontsize)
        self.draw()


class MyMplCanvas(Canvas):
    def __init__(self, parent=None, width=8, height=5):
        self.showf1 = False
        self.plotData = None
        fig = Figure(figsize=(width, height))
        self.fig = fig
        self.axes = fig.add_subplot(111)
        self.axesr = self.axes.twinx()
        self.fig.subplots_adjust(left=0.12, right=0.88, bottom=0.15, top=0.97)
#        self.axes.hold(False)  # clear axes every time plot() is called
        Canvas.__init__(self, fig)
        self.setParent(parent)
        self.updateGeometry()
        fm = QtGui.QFontMetrics(self.font())
        self.fontsize = int(fm.height()) / 1.33

        self.axes.set_xlabel('energy (keV)', fontsize=self.fontsize)
        self.axes.set_ylabel(r"$f''$", fontsize=self.fontsize)
        self.axesr.set_ylabel(r"$f'$", fontsize=self.fontsize)
        self.axesr.set_visible(self.showf1)

    def plot(self, compound, E, table):
        self.plotData = compound, E, table
        for artist in self.axes.lines + self.axes.collections:
            artist.remove()
        if self.showf1:
            for artist in self.axesr.lines + self.axesr.collections:
                artist.remove()
        self.axesr.set_visible(self.showf1)
        if compound is None:
            ll = self.axes.legend([], title='no given elements',
                                  loc='upper right', fontsize=self.fontsize)
        else:
            ced = xc.calculate_element_dict(compound.asList(), E, table)
            if isinstance(ced, str):
                QtWidgets.QMessageBox.critical(self, "Error", ced)
                return
            elementsDict = ced[0]
            for iel, (elName, el) in enumerate(elementsDict.items()):
                color = 'C{0}'.format(iel % 10)
                label = elName if el[5] == 0 else\
                    r"{0}, $f''$={1}, $\Delta f''$={2}".format(
                        elName, xc.round_to_n(el[4]), xc.round_to_n(el[5]))
                self.axes.plot(el[0].E*1e-3, el[0].f2, '-', marker='.',
                               label=label, color=color)
                if self.showf1:
                    self.axesr.plot(el[0].E*1e-3, el[0].f1, '--', marker='.',
                                    color=color)
            # self.axes.set_xlim(el[0].E[0]*1e-3, el[0].E[-1]*1e-3)
            ll = self.axes.legend(loc='upper right', fontsize=self.fontsize)
            ylim = self.axes.get_ylim()
            dy = ylim[1] - ylim[0]
            ylim = ylim[0] + 0.2*dy, ylim[1] - 0.2*dy
            self.axes.plot([E*1e-3, E*1e-3], ylim, '--', color='gray',
                           label=None)
        mpl.artist.setp(ll.get_title(), fontsize=self.fontsize)
        self.draw()


class PlotDlg(QtWidgets.QDialog):
    def __init__(self, parent, compound, E, table):
        super(PlotDlg, self).__init__(parent)
        bl = QtWidgets.QVBoxLayout(self)
        self.plotCanvas = MyMplCanvas(self)
        self.plotCanvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                      QtWidgets.QSizePolicy.Expanding)
        self.toolbar = ToolBar(self.plotCanvas, self)
        self.f1CB = QtWidgets.QCheckBox("show f'")
        self.f1CB.clicked.connect(self.showf1)
        self.toolbar.insertWidget(self.toolbar.actions()[-1], self.f1CB)

        bl.addWidget(self.toolbar)
        bl.addWidget(self.plotCanvas)
        pg = parent.frameGeometry()
        self.move(parent.x()+pg.width(), parent.y())
        pg = parent.geometry()
        self.resize(int(pg.width()*1.8), int(pg.height()))
        self.setWindowTitle("plots of f' and f''")
        self.setWindowFlags(QtCore.Qt.Window)
        self.show()
        self.plotCanvas.plot(compound, E, table)

#    def closeEvent(self, event):  # is not invoked by esc. press
    def done(self, event):
        self.parent().plotDlg = None
        super(PlotDlg, self).done(event)

    def showf1(self, checked):
        self.plotCanvas.showf1 = checked
        if self.plotCanvas.plotData is None:
            return
        self.plotCanvas.plot(*self.plotCanvas.plotData)


class MainDlg(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(MainDlg, self).__init__(parent)

        self.whatLabel = QtWidgets.QLabel(r"sample:")
        self.whatCB = QtWidgets.QComboBox()
        self.whatCB.addItems(whats)
        self.whatCB.currentIndexChanged.connect(self.updateUi)
        self.what = 0

        self.formula = MyFormulaMplCanvas(self, width=3.2, height=0.5)

        self.compoundLabel = QtWidgets.QLabel("compound:")
        self.compoundExLabel = QtWidgets.QLabel("")
        self.compoundExLabel.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlags(QtCore.Qt.TextSelectableByMouse))
        self.compoundEdit = QtWidgets.QLineEdit()
        self.compoundEdit.setPlaceholderText(
            'type here or select from the lists above')

        self.compoundList1 = ComboBoxWithPlaceholder()
        self.compoundList1.addItems(
            [k for k in OrderedDict(
                sorted(compounds.items(), key=lambda it: it[1]['rho']))]),
        try:
            self.compoundList1.setPlaceholderText(
                    "compounds (sorted by density)")
        except AttributeError:
            pass
        self.compoundList1.setToolTip("sorted by density")
        self.compoundList1.setMaxVisibleItems(18)
        self.compoundList1.setCurrentIndex(-1)
        self.compoundList1.activated.connect(self.compoundActivated)
        self.compoundList2 = ComboBoxWithPlaceholder()
        self.compoundList2.addItems(
            [k for k in OrderedDict(  # have to sort it in Py2
                sorted(elemental.items(), key=lambda it: it[1]['Z']))]),
        try:
            self.compoundList2.setPlaceholderText("elements")
        except AttributeError:
            pass
        self.compoundList2.setMaxVisibleItems(18)
        self.compoundList2.setCurrentIndex(-1)
        self.compoundList2.activated.connect(self.elementActivated)

        self.compoundMassLabel = QtWidgets.QLabel("M (g/mol) = ")
        self.compoundMass = QtWidgets.QLabel("")

        self.muTdLabel = QtWidgets.QLabel("")
        self.muTdEdit = QtWidgets.QLineEdit()
#        font = QtGui.QFont(self.font())
#        font.setPointSize(font.pointSize()+2)
#        fm = QtGui.QFontMetrics(font)
#        self.muTdEdit.setMinimumSize(fm.width("8.88"), fm.height())
        self.areaLabel = QtWidgets.QLabel("")
        self.areaEdit = QtWidgets.QLineEdit()
#        self.areaEdit.setMinimumSize(fm.width("8.88"), fm.height())
        self.dmudLabel = QtWidgets.QLabel(u"δµd = ")
        self.dmudEdit = QtWidgets.QLineEdit()
#        self.dmudEdit.setMinimumSize(fm.width("8.888"), fm.height())

        self.energyLabel = QtWidgets.QLabel("E (eV) =")
        self.energyCB = QtWidgets.QComboBox()
        self.read_energies()
        self.energyCB.addItems(self.energies)
        self.energyCB.setEditable(True)
        self.energyCB.currentIndexChanged.connect(self.energySelected)
        self.energy = 9029
        self.energyCB.lineEdit().setText(str(self.energy))
        self.energyCB.setMaxVisibleItems(25)

        self.tableLabel = QtWidgets.QLabel("f data table:")
        self.tableCB = QtWidgets.QComboBox()
        self.tableCB.addItems(tables)
        self.tableCB.setCurrentIndex(3)  # Chantler total
        self.tableCB.currentIndexChanged.connect(self.calculate)
        self.tablePlotButton = QtWidgets.QPushButton('Plot f"')
        if not MAC:
            self.tablePlotButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.tablePlotButton.clicked.connect(self.plotf)

        self.resNuLabel = QtWidgets.QLabel("")
        self.resNu = QtWidgets.QLabel("")
        self.resMassLabel = QtWidgets.QLabel("")
        self.resMass = QtWidgets.QLabel("")

        self.stepLabel = QtWidgets.QLabel("absorptance step = ")
        self.stepCB = QtWidgets.QComboBox()
        self.stepCB.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                  QtWidgets.QSizePolicy.Fixed)
#        self.stepCB.setMinimumWidth(200)

        self.extraRhoLabel = QtWidgets.QLabel(u"ρ (g/cm<sup>3</sup>) = ")
        self.extraRhoEdit = QtWidgets.QLineEdit()
#        self.extraRhoEdit.setMinimumSize(fm.width("8.88"), fm.height())
        self.extraRhoEdit.setSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,
                                        QtWidgets.QSizePolicy.Fixed)
        self.extraDLabel = QtWidgets.QLabel(u"d (µm) = ")
        self.extraD = QtWidgets.QLabel(u"")

        self.buttonCalculate = QtWidgets.QPushButton("Calculate")
        self.buttonCalculate.clicked.connect(self.calculate)
        self.buttonAbout = QtWidgets.QPushButton("About...")
        self.buttonAbout.clicked.connect(self.about)
        self.buttonHelp = QtWidgets.QPushButton("Help...")
        self.buttonHelp.clicked.connect(self.myhelp)

        whatLayout = QtWidgets.QHBoxLayout()
        whatLayout.addWidget(self.whatLabel)
        whatLayout.addWidget(self.whatCB)
        whatLayout.addStretch()

        compoundLayout = QtWidgets.QHBoxLayout()
        compoundLayout.addWidget(self.compoundLabel)
        compoundLayout.addWidget(self.compoundExLabel)
        compoundLayout.addStretch()

        listsLayout = QtWidgets.QHBoxLayout()
        # listsLayout.addWidget(self.compoundListsLabel)
        listsLayout.addWidget(self.compoundList1)
        listsLayout.addWidget(self.compoundList2)

        compoundMassLayout = QtWidgets.QHBoxLayout()
        compoundMassLayout.addWidget(self.compoundMassLabel)
        compoundMassLayout.addWidget(self.compoundMass)
        compoundMassLayout.addStretch()

        muTdLayout = QtWidgets.QHBoxLayout()
        muTdLayout.addWidget(self.muTdLabel)
        muTdLayout.addWidget(self.muTdEdit)
        muTdLayout.addStretch()
        muTdLayout.addWidget(self.areaLabel)
        muTdLayout.addWidget(self.areaEdit)
        muTdLayout.addStretch()
        muTdLayout.addWidget(self.dmudLabel)
        muTdLayout.addWidget(self.dmudEdit)
        muTdLayout.addStretch()

        energyLayout = QtWidgets.QHBoxLayout()
        energyLayout.addWidget(self.energyLabel)
        energyLayout.addWidget(self.energyCB)
        energyLayout.addStretch()

        tableLayout = QtWidgets.QHBoxLayout()
        tableLayout.addWidget(self.tableLabel)
        tableLayout.addWidget(self.tableCB)
        tableLayout.addStretch()
        tableLayout.addWidget(self.tablePlotButton)
        tableLayout.addStretch()

        resLayout = QtWidgets.QHBoxLayout()
        resLayout.addWidget(self.resNuLabel)
        resLayout.addWidget(self.resNu)
        resLayout.addStretch()
        resLayout.addWidget(self.resMassLabel)
        resLayout.addWidget(self.resMass)
        resLayout.addStretch()

        stepLayout = QtWidgets.QHBoxLayout()
        stepLayout.addWidget(self.stepLabel)
        stepLayout.addWidget(self.stepCB)
#        stepLayout.addStretch()

        extraLayout = QtWidgets.QHBoxLayout()
        extraLayout.addWidget(self.extraRhoLabel)
        extraLayout.addWidget(self.extraRhoEdit)
        extraLayout.addStretch()
        extraLayout.addWidget(self.extraDLabel)
        extraLayout.addWidget(self.extraD)
        extraLayout.addStretch()

        self.gasesInMixture = []
        self.gasesInMixtureN = 0
        if gasMixerN > 0:
            mixerLayout = QtWidgets.QVBoxLayout()
            mixerLayout.setContentsMargins(0, 0, 0, 0)
            for iGas in range(gasMixerN):
                gasInMixture = QtWidgets.QWidget()
                self.gasesInMixture.append(gasInMixture)
                gasLayout = QtWidgets.QHBoxLayout(gasInMixture)
                gasLayout.setContentsMargins(0, 0, 0, 0)

                gasLabel = QtWidgets.QLabel()

                gasSlider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
                gasSlider.setMinimum(gasMixerPressureMin)
                gasSlider.setMaximum(gasMixerPressureMax)
                gasSlider.setTickPosition(QtWidgets.QSlider.TicksAbove)
                gasSlider.setTickInterval(gasMixerPressurePageStep)
                gasSlider.setPageStep(gasMixerPressurePageStep)
                gasSlider.setValue(gasMixerPressureDefault)
                gasSlider.valueChanged.connect(
                    partial(self.gasSliderChanged, iGas))

                gasValue = QtWidgets.QLabel(
                    "{0:.0f} mbar".format(gasMixerPressureDefault))
                gasInMixture.pressure = gasMixerPressureDefault

                gasLayout.addWidget(gasLabel)
                gasLayout.addWidget(gasSlider)
                gasLayout.addWidget(gasValue)
                mixerLayout.addWidget(gasInMixture)
            self.calculateGasMixture()

        buttonLayout = QtWidgets.QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.buttonCalculate)  # 1st button gets Enter
        buttonLayout.addWidget(self.buttonAbout)
        buttonLayout.addWidget(self.buttonHelp)
        buttonLayout.addStretch()

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(whatLayout)
        layout.addWidget(self.formula)
        hline1 = QtWidgets.QFrame(self)
        hline1.setFrameStyle(QtWidgets.QFrame.HLine)
        layout.addWidget(hline1)
        layout.addLayout(compoundLayout)
        layout.addLayout(listsLayout)
        layout.addWidget(self.compoundEdit)
        layout.addLayout(compoundMassLayout)
        hline2 = QtWidgets.QFrame(self)
        hline2.setFrameStyle(QtWidgets.QFrame.HLine)
        layout.addWidget(hline2)
        layout.addLayout(muTdLayout)
        layout.addLayout(energyLayout)
        hline3 = QtWidgets.QFrame(self)
        hline3.setFrameStyle(QtWidgets.QFrame.HLine)
        layout.addWidget(hline3)
        layout.addLayout(tableLayout)
        hline4 = QtWidgets.QFrame(self)
        hline4.setFrameStyle(QtWidgets.QFrame.HLine)
        layout.addWidget(hline4)
        layout.addLayout(resLayout)
        layout.addLayout(stepLayout)
        layout.addLayout(extraLayout)
        if gasMixerN > 0:
            layout.addLayout(mixerLayout)
        layout.addStretch()
        layout.addLayout(buttonLayout)
        self.setLayout(layout)

        self.setWindowTitle("XAFSmassQt")
        selfDir = os.path.dirname(__file__)
        icon = QtGui.QIcon(os.path.abspath(os.path.join(
            selfDir, 'help', '_static', 'XAFSmassQt.ico')))
        self.setWindowIcon(icon)
        self.setWindowFlags(QtCore.Qt.Window)
        try:
            xy = QtWidgets.QApplication.primaryScreen().availableGeometry().center() - \
                self.rect().center()
        except AttributeError:
            xy = QtWidgets.QApplication.desktop().screen().rect().center() - \
                self.rect().center()
        self.move(xy)
        self.updateUi()
        self.plotDlg = None
#        pg = self.geometry()
#        self.formula.setFixedHeight(pg.width()*0.5/4)

    def sizeHint(self):
        return QtCore.QSize(0, 0)  # set to minimum possible

    def updateUi(self):
        self.what = self.whatCB.currentIndex()
        self.formula.update_formula(formulas[self.what])
        self.compoundExLabel.setText("<i>e.g.</i> " + examples[self.what])

        if self.what in [POWDER, FOIL, XCONTENT]:
            self.muTdLabel.setText(u"µ<sub>T</sub>d = ")
        elif self.what == GAS:
            self.muTdLabel.setText(u"attenuation = ")

        if self.what in [POWDER, FOIL, GAS]:
            self.dmudLabel.hide()
            self.dmudEdit.hide()
        elif self.what == XCONTENT:
            self.dmudLabel.show()
            self.dmudEdit.show()

        if self.what == POWDER:
            self.areaLabel.setText(u"S (cm<sup>2</sup>) = ")
            self.resNuLabel.setText(u"ν (mmol) = ")
            self.resMassLabel.setText("m (mg) = ")
            self.muTdEdit.setPlaceholderText('typ. 2.6')
            self.areaEdit.setPlaceholderText('1.33 for \u230013mm ')
            self.resNuLabel.show()
            self.resNu.show()
        elif self.what == FOIL:
            self.areaLabel.setText(u"ρ (g/cm<sup>3</sup>) = ")
            self.resMassLabel.setText(u"d (µm) = ")
            self.muTdEdit.setPlaceholderText('typ. 2.6')
            self.areaEdit.setPlaceholderText('')
            self.resNuLabel.hide()
            self.resNu.hide()
        elif self.what == GAS:
            self.areaLabel.setText(u"d (cm) = ")
            self.resMassLabel.setText(u"p (mbar) = ")
            self.muTdEdit.setPlaceholderText('(0, 1)')
            self.areaEdit.setPlaceholderText('')
            self.resNuLabel.hide()
            self.resNu.hide()
        elif self.what == XCONTENT:
            self.areaLabel.setText(u"Δµd = ")
            self.resNuLabel.setText("N<sub>x</sub> = ")
            self.resMassLabel.setText(r"wt%<sub>x</sub> = ")
            self.muTdEdit.setPlaceholderText('')
            self.areaEdit.setPlaceholderText('')
            self.resNuLabel.show()
            self.resNu.show()

        if self.what in [POWDER, FOIL]:
            self.stepLabel.show()
            self.stepCB.show()
        else:
            self.stepLabel.hide()
            self.stepCB.hide()

        if self.what == POWDER:
            self.extraRhoLabel.show()
            self.extraRhoEdit.show()
            self.extraDLabel.show()
            self.extraD.show()
        else:
            self.extraRhoLabel.hide()
            self.extraRhoEdit.hide()
            self.extraDLabel.hide()
            self.extraD.hide()

        for gas in self.gasesInMixture:
            gas.hide()

        self.resize(0, 0)
        self.calculate()

    def setupGasMixtures(self, action='parse'):
        if action == 'clear':
            for gas in self.gasesInMixture:
                gas.hide()
            return
        if len(self.gasesInMixture) == 0:
            return
        compound = str(self.compoundEdit.text())
        groups = re.findall(r'\((.*?)\)', compound)
        self.gasesInMixtureN = len(groups)
        if self.gasesInMixtureN == 0:
            self.setupGasMixtures('clear')
            return
        for g, gas in zip(groups, self.gasesInMixture):
            gasLabel = gas.layout().itemAt(0).widget()
            gasLabel.setText(g)
            gas.show()
        for gas in self.gasesInMixture[len(groups):]:
            gas.hide()

    def gasSliderChanged(self, gasNo, value):
        gasValue = self.gasesInMixture[gasNo].layout().itemAt(2).widget()
        self.gasesInMixture[gasNo].pressure = value
        gasValue.setText("{0:.0f} mbar".format(value))
        self.calculateGasMixture()

    def calculateGasMixture(self):
        if self.gasesInMixtureN == 0:
            return
        try:
            E = float(self.energyCB.lineEdit().text())
            length = float(self.areaEdit.text())
        except ValueError:
            return
        table = tablesF[self.tableCB.currentIndex()]
        nus2 = 0.
        sumPressure = 0.
        for gasNo in range(self.gasesInMixtureN):
            gasLabel = self.gasesInMixture[gasNo].layout().itemAt(0).widget()
            compound = str(gasLabel.text())
            parsedCompound = xc.formula.parseString(compound, parseAll=True)
            formulaList = parsedCompound.asList()
            if isinstance(formulaList[0], dict):
                sumSigma2 = formulaList[1]
            else:
                res = xc.calculate_element_dict(formulaList, E, table=table)
                if isinstance(res, str):
                    return res
                sumSigma2 = res[1]
            nu = self.gasesInMixture[gasNo].pressure / (xc.R*xc.T)*length*1e-4
            nus2 += nu * sumSigma2
            sumPressure += self.gasesInMixture[gasNo].pressure
        attenuation = 1 - np.exp(-nus2)
        self.muTdEdit.setText("{0:.3f}".format(attenuation))
        self.resMass.setText("<b>{0:.0f}</b>".format(sumPressure))
        compound = ""
        for gasNo in range(self.gasesInMixtureN):
            gasLabel = self.gasesInMixture[gasNo].layout().itemAt(0).widget()
            normPressure = self.gasesInMixture[gasNo].pressure/sumPressure
            if abs(normPressure - 1) < 1e-3:
                textNormPressure = ''
            else:
                textNormPressure = xc.round_to_n(normPressure, 3)
            compound += "({0}){1}".format(gasLabel.text(), textNormPressure)
        self.compoundEdit.setText(compound)

    def read_energies(self):
        selfDir = os.path.dirname(__file__)
        efname = os.path.join(selfDir, 'data', 'Energies.txt')
        with open(efname, 'r') as f:
            f.readline()
            f.readline()
            self.energies = []
            for line in f.readlines():
                cs = line.strip().split()
                if len(cs[0]) == 1:
                    cs[0] = '0' + cs[0]
                pre = cs[0] + ' ' + cs[1] + ' '
                for ic, c in enumerate(cs[2:]):
                    self.energies.append(pre + edges[ic] + ' ' + c + ' + 50')

    def compoundActivated(self):
        txt = str(self.compoundList1.currentText())  # otherwise u'' in Py2
        self.compoundEdit.setText(compounds[txt]['formula'])
        if self.what in [POWDER, FOIL]:
            rho = compounds[txt]['rho']
            if self.what == POWDER:
                self.extraRhoEdit.setText(str(rho))
            elif self.what == FOIL:
                self.areaEdit.setText(str(rho))
        self.calculate()

    def elementActivated(self):
        txt = str(self.compoundList2.currentText())  # otherwise u'' in Py2
        self.compoundEdit.setText(txt)
        if self.what in [POWDER, FOIL]:
            rho = elemental[txt]['rho']
            if self.what == POWDER:
                self.extraRhoEdit.setText(str(rho))
            elif self.what == FOIL:
                self.areaEdit.setText(str(rho))
        self.calculate()

    def energySelected(self):
        txt = self.energyCB.lineEdit().text()
        try:
            self.energy = float(txt)
        except ValueError:
            st = str(txt).strip().split()
            if len(st) == 6:
                e = float(st[3]) + float(st[5])
                e = int(np.floor(e)) if e == np.floor(e) else e
                self.energy = e
                self.energyCB.lineEdit().setText(str(self.energy))
        self.calculate()

    def clear_results(self):
        self.resNu.setText("")
        self.resMass.setText("")
        self.extraD.setText("")

    def parse_compound(self):
        compound = str(self.compoundEdit.text())
        self.parsedCompound = None
        if len(compound) == 0:
            self.compoundMass.setText("")
            return
        try:
            parsedCompound = xc.formula.parseString(compound, parseAll=True)
            withX = False
            cMass = 0.
            for e in parsedCompound.asList():
                if e[1] < 0:
                    withX = True
                else:
                    cMass += xc.rm.read_atomic_data(e[0])*e[1]
            if withX:
                self.compoundMass.setText("(of matrix) {0}".format(
                    xc.round_to_n(cMass, 5)))
            elif '%' in compound:
                self.compoundEdit.setText(xc.reconstruct(parsedCompound))
                self.calculate()
                return
            else:
                self.compoundMass.setText(
                    '{0}'.format(xc.round_to_n(cMass, 5)))
        except (ParseBaseException, ValueError, ZeroDivisionError):
            self.compoundMass.setText("wrong compound formula")
            return

        if (self.what == XCONTENT) ^ withX:
            if (self.what == XCONTENT):
                outStr = ", give 'x'"
            else:
                outStr = ", remove 'x'"
            self.compoundMass.setText("wrong compound formula"+outStr)
            self.clear_results()
            return

        self.parsedCompound = parsedCompound
        return True

    def calculate(self):
        if not self.parse_compound():
            if self.what == GAS:
                self.setupGasMixtures('clear')
            return
        if self.what == GAS:
            self.setupGasMixtures()
        try:
            E = float(self.energyCB.lineEdit().text())
        except ValueError:
            return
        table = tablesF[self.tableCB.currentIndex()]
        if self.plotDlg is not None:
            self.plotDlg.plotCanvas.plot(self.parsedCompound, E, table)

        self.stepCB.clear()
        try:
            xxx = float(self.muTdEdit.text())
            if xxx < 0:
                xxx = abs(xxx)
                self.muTdEdit.setText(repr(xxx))
            yyy = float(self.areaEdit.text())
            if yyy < 0:
                yyy = abs(yyy)
                self.areaEdit.setText(repr(yyy))
            if self.what == POWDER:
                muTd = xxx
                area = yyy
                rho = self.extraRhoEdit.text()
                if len(rho) == 0:
                    rho = None
                    self.extraD.setText("")
                else:
                    rho = float(rho)
                    if rho < 0:
                        rho = abs(rho)
                        self.extraRhoEdit.setText(repr(rho))
            elif self.what == FOIL:
                muTd = xxx
                rho = yyy
            elif self.what == GAS:
                attenuation = xxx
                if not (0 < attenuation < 1):
                    self.muTdEdit.setText("")
                    self.muTdEdit.setFocus()
                    raise ValueError
                thickness = yyy
            elif self.what == XCONTENT:
                muTd = xxx
                Dmud = yyy
        except ValueError:
            self.clear_results()
            return

        if self.what == XCONTENT:
            try:
                dmud = float(self.dmudEdit.text())
            except ValueError:
                dmud = 0

        if self.what == POWDER:
            res = xc.calculate_powder(
                self.parsedCompound.asList(), E, muTd, area, rho, table=table)
            if isinstance(res, str):
                self.clear_results()
                QtWidgets.QMessageBox.critical(self, "Error", res)
                return
            nu, m, th, eDict = res
            self.resNu.setText('{0}'.format(xc.round_to_n(nu, 3)))
            self.resMass.setText("<b>{0}</b>".format(xc.round_to_n(m, 3)))
            if rho is not None:
                self.extraD.setText('{0}'.format(xc.round_to_n(th, 4)))
            iSelect = 0
            for i, (k, v) in enumerate(eDict.items()):
                self.stepCB.addItem("{0}({1}mg={2}wt%): {3}".format(
                    k, xc.round_to_n(v[-2], 3), xc.round_to_n(v[-2]/m*100., 2),
                    xc.round_to_n(v[-1], 3)))
                if v[2] > 0:
                    iSelect = i
            self.stepCB.setCurrentIndex(iSelect)
        elif self.what == FOIL:
            res = xc.calculate_foil(
                self.parsedCompound.asList(), E, muTd, rho, table=table)
            if isinstance(res, str):
                self.clear_results()
                QtWidgets.QMessageBox.critical(self, "Error", res)
                return
            th, eDict = res
            self.resMass.setText("<b>{0}</b>".format(xc.round_to_n(th, 3)))
            iSelect = 0
            for i, (k, v) in enumerate(eDict.items()):
                self.stepCB.addItem("{0}: {1}".format(
                    k, xc.round_to_n(v[-1], 3)))
                if v[1] > 0:
                    iSelect = i
            self.stepCB.setCurrentIndex(iSelect)
        elif self.what == GAS:
            res = xc.calculate_gas(
                self.parsedCompound.asList(), E, attenuation, thickness,
                table=table)
            if isinstance(res, str):
                self.clear_results()
                QtWidgets.QMessageBox.critical(self, "Error", res)
                return
            P = res
            self.resMass.setText("<b>{0}</b>".format(xc.round_to_n(P, 3)))
        elif self.what == XCONTENT:
            res = xc.calculate_x(
                self.parsedCompound.asList(), E, muTd, Dmud, dmud, table=table)
            if isinstance(res, str):
                self.clear_results()
                QtWidgets.QMessageBox.critical(self, "Error", res)
                return
            nu, wt, wt1 = res
            self.resNu.setText('{0}'.format(xc.round_to_n(nu, 3)))
            outStr = '{0}'.format(xc.round_to_n(wt, 5))
            if wt1 > 0:
                outStr += u'±{0}'.format(xc.round_to_n(wt1, 3))
            self.resMass.setText(u"<b>{0}</b>".format(outStr))

    def about(self):
        # https://stackoverflow.com/a/69325836/2696065
        def isWin11():
            return True if sys.getwindowsversion().build > 22000 else False

        locos = platform.platform(terse=True)
        if 'Linux' in locos:
            try:
                locos = " ".join(platform.linux_distribution())
            except AttributeError:  # no platform.linux_distribution in py3.8
                try:
                    import distro
                    locos = " ".join([distro.name(), distro.version()])
                except ImportError:
                    print("do 'pip install distro' for a better view of Linux"
                          " distro string")
        elif 'Windows' in locos:
            if isWin11():
                locos = 'Winows 11'

        title = "About XAFSmassQt"
        txt = \
            """<b>XAFSmass(Qt)</b> v {0}
            <ul>
            <li>{1[0]}
            <li>{1[1]}
            </ul>
            <p>Open source, {2}. Available on PyPI and GitHub<p>
            <p>Your system:
            <ul>
            <li>{3}
            <li>Python {4}
            <li>Qt: {5} {6}
            </ul>""".format(
                __version__, __author__.split(','), __license__,
                locos, platform.python_version(),
                qtpy.API_NAME, qtpy.QT_VERSION)

        msg = QtWidgets.QMessageBox()
        msg.setStyleSheet("#qt_msgbox_label{min-width: 400px;}")
        msg.setWindowIcon(self.windowIcon())
        msg.setIconPixmap(QtGui.QPixmap(os.path.join(
             'help', '_static', 'XAFSmassQt.ico')))
        msg.setText(txt)
        msg.setWindowTitle(title)
        msg.exec_()

    def myhelp(self):
        hname = os.path.join(os.path.dirname(__file__), 'help', 'index.html')
        webbrowser.open(hname)

    def plotf(self):
        self.parse_compound()
        E = float(self.energyCB.lineEdit().text())
        table = tablesF[self.tableCB.currentIndex()]
        if self.plotDlg is None:
            self.plotDlg = PlotDlg(self, self.parsedCompound, E, table)
            self.plotDlg.setMinimumWidth(600)
        else:
            self.plotDlg.plotCanvas.plot(self.parsedCompound, E, table)


def run():
    app = QtWidgets.QApplication(sys.argv)
    form = MainDlg()
    form.show()
    app.exec_()


if __name__ == "__main__":
    run()
