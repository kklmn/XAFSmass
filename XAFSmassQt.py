# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev, Roman Chernikov"
__date__ = "15 May 2017"

import sys
import os
import webbrowser
import numpy as np
from pyparsing import ParseBaseException

#=this may help to resolve conflict betwen Qt4 and Qt5=========================
# import matplotlib as mpl
# mpl.use("Qt4Agg")
#==============================================================================

try:
    from matplotlib.backends import qt_compat
except ImportError:
    from matplotlib.backends import qt4_compat
    qt_compat = qt4_compat
use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE
if use_pyside:
    QtName = "PySide"
    import PySide
    from PySide import QtGui, QtCore
    import PySide.QtGui as myQtGUI
    import matplotlib.backends.backend_qt4agg as mpl_qt
else:
    try:
        QtName = "PyQt4"
        from PyQt4 import QtGui, QtCore
        import PyQt4.QtGui as myQtGUI
        import matplotlib.backends.backend_qt4agg as mpl_qt
    except ImportError:
        QtName = "PyQt5"
        from PyQt5 import QtGui, QtCore
        import PyQt5.QtWidgets as myQtGUI
        import matplotlib.backends.backend_qt5agg as mpl_qt

QDialog, QApplication, QLabel, QComboBox, QLineEdit, QPushButton,\
    QSizePolicy, QHBoxLayout, QVBoxLayout, QFrame, QMessageBox = \
    myQtGUI.QDialog, myQtGUI.QApplication, myQtGUI.QLabel,\
    myQtGUI.QComboBox, myQtGUI.QLineEdit, myQtGUI.QPushButton,\
    myQtGUI.QSizePolicy, myQtGUI.QHBoxLayout, myQtGUI.QVBoxLayout,\
    myQtGUI.QFrame, myQtGUI.QMessageBox
Canvas = mpl_qt.FigureCanvasQTAgg
ToolBar = mpl_qt.NavigationToolbar2QT

from matplotlib.figure import Figure
import matplotlib as mpl
# Set the default color cycle
#mpl.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'm', 'c', 'k']
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',
# ...this to force siunitx to actually use your fonts:
#       r'\sisetup{detect-all}',
# set the normal font here:
#       r'\usepackage{helvet}',
# load up the sansmath so that math -> helvet:
#       r'\usepackage{sansmath}',
# <- tricky! -- gotta actually tell tex to use!:
#       r'\sansmath'
#]

import XAFSmassCalc as xc
from __init__ import (__version__, __author__, __license__)

MAC = "qt_mac_set_native_menubar" in dir()

POWDER, FOIL, GAS, XCONTENT = range(4)
whats = ('powder', 'foil, film, glass etc.', 'gas',
         'has unknown concentration')
formulas = (
    r"$\nu = (\mu_T d)\cdot S\cdot\left(\sum_i{N_AN_i2r_0\lambda f_i''}"
    r"\right)^{-1}; \quad m = M\cdot\nu$",
    r"$d = (\mu_T d)\cdot M\cdot\left(\rho\sum_i{N_AN_i2r_0\lambda f_i''}"
    r"\right)^{-1}$",
    r"$p = -\ln(1-{\rm attenuation})\cdot kT \cdot"
    r"\left(d\sum_i{N_i2r_0\lambda f_i''}\right)^{-1}$",
    r"$\Delta\mu/\mu_T = N_x\Delta f_x'' \cdot"
    r"\left(\sum_{i\ne x}{N_if_i''} + N_xf_x''\right)^{-1}$")
examples = (
    "Cu(NO3)2  <i>or</i>  Cu%1Zn%1((Al2O3)%10SiO2)",
    "Cu%25Zn", "Ar0.9Kr0.1  <i>or</i>  N2", "FexSiO2")
tables = ("Henke", "Brennan&Cowan", "Chantler (NIST)", "Chantler total (NIST)")
tablesF = ("Henke", "BrCo", "Chantler", "Chantler total")
edges = ("K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "N1", "N2", "N3")


class MyFormulaMplCanvas(Canvas):
    def __init__(self, parent=None, width=5, height=0.4):
        fig = Figure(figsize=(width, height), dpi=96)
        self.fig = fig
        Canvas.__init__(self, fig)
        bg = self.palette().window().color()
        cl = (bg.redF(), bg.greenF(), bg.blueF())
#        fig.set_edgecolor(cl)
        fig.set_facecolor(cl)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.updateGeometry()
        fm = QtGui.QFontMetrics(self.font())
        self.fontsize = int(fm.height()) / 1.25

    def update_formula(self, formula=None):
        self.fig.clf()
        self.fig.suptitle(formula, x=0.5, y=0.48, ha='center', va='center',
                          fontsize=self.fontsize)
        self.draw()


class MyMplCanvas(Canvas):
    def __init__(self, parent=None, width=6, height=5):
        fig = Figure(figsize=(width, height))
        self.fig = fig
        self.axes = fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.97)
#        self.axes.hold(False)  # clear axes every time plot() is called
        Canvas.__init__(self, fig)
        self.setParent(parent)
        self.updateGeometry()
        fm = QtGui.QFontMetrics(self.font())
        self.fontsize = int(fm.height()) / 1.25

    def plot(self, compound, E, table):
        self.axes.cla()
        self.axes.set_xlabel('energy (keV)', fontsize=self.fontsize)
        self.axes.set_ylabel(r"$f''$", fontsize=self.fontsize)
        if compound is None:
            ll = self.axes.legend([], title='no given elements',
                                  loc='upper right', fontsize=self.fontsize)
        else:
            ced = xc.calculate_element_dict(compound.asList(), E, table)
            if isinstance(ced, str):
                QMessageBox.critical(self, "Error", ced)
                return
            elementsDict = ced[0]
            for elName, el in elementsDict.items():
                label = elName if el[5] == 0 else\
                    r"{0}, $f''$={1}, $\Delta f''$={2}".format(
                        elName, xc.round_to_n(el[4]), xc.round_to_n(el[5]))
                self.axes.plot(el[0].E*1e-3, el[0].f2, '-', marker='.',
                               label=label)
            self.axes.set_xlim(el[0].E[0]*1e-3, el[0].E[-1]*1e-3)
            ll = self.axes.legend(loc='upper right', fontsize=self.fontsize)
            ylim = self.axes.get_ylim()
            self.axes.plot([E*1e-3, E*1e-3], ylim, '--', color='gray',
                           label=None)
        mpl.artist.setp(ll.get_title(), fontsize=self.fontsize)
        self.draw()


class PlotDlg(QDialog):
    def __init__(self, parent, compound, E, table):
        super(PlotDlg, self).__init__(parent)
        bl = QVBoxLayout(self)
        self.plotCanvas = MyMplCanvas(self)
        self.plotCanvas.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = ToolBar(self.plotCanvas, self)
        bl.addWidget(self.toolbar)
        bl.addWidget(self.plotCanvas)
        pg = parent.frameGeometry()
        self.move(parent.x()+pg.width(), parent.y())
        pg = parent.geometry()
        self.resize(pg.width()*1.5, pg.height())
        self.setWindowTitle("plots of f''")
        self.setWindowFlags(QtCore.Qt.Window)
        self.show()
        self.plotCanvas.plot(compound, E, table)

#    def closeEvent(self, event):  # is not invoked by esc. press
    def done(self, event):
        self.parent().plotDlg = None
        super(PlotDlg, self).done(event)


class MainDlg(QDialog):
    def __init__(self, parent=None):
        super(MainDlg, self).__init__(parent)

        self.whatLabel = QLabel(r"&sample:")
        self.whatCB = QComboBox()
        self.whatCB.addItems(whats)
        self.whatCB.currentIndexChanged.connect(self.updateUi)
        self.whatLabel.setBuddy(self.whatCB)
        self.what = 0

        self.formula = MyFormulaMplCanvas(self, width=3.2, height=0.6)

        self.compoundLabel = QLabel(r"&compound:")
        self.compoundExLabel = QLabel("")
        self.compoundExLabel.setTextInteractionFlags(
            QtCore.Qt.TextInteractionFlags(QtCore.Qt.TextSelectableByMouse))
        self.compoundEdit = QLineEdit()
        self.compoundLabel.setBuddy(self.compoundEdit)

        self.compoundMassLabel = QLabel("M (g/mol) = ")
        self.compoundMass = QLabel("")

        self.muTdLabel = QLabel("")
        self.muTdEdit = QLineEdit()
#        font = QtGui.QFont(self.font())
#        font.setPointSize(font.pointSize()+2)
#        fm = QtGui.QFontMetrics(font)
#        self.muTdEdit.setMinimumSize(fm.width("8.88"), fm.height())
        self.areaLabel = QLabel("")
        self.areaEdit = QLineEdit()
#        self.areaEdit.setMinimumSize(fm.width("8.88"), fm.height())
        self.dmudLabel = QLabel(u"δµd = ")
        self.dmudEdit = QLineEdit()
#        self.dmudEdit.setMinimumSize(fm.width("8.888"), fm.height())

        self.energyLabel = QLabel(r"&E (eV) =")
        self.energyCB = QComboBox()
        self.read_energies()
        self.energyCB.addItems(self.energies)
        self.energyCB.setEditable(True)
        self.energyCB.currentIndexChanged.connect(self.energyCB_selected)
        self.energyLabel.setBuddy(self.energyCB)
        self.energy = 9029
        self.energyCB.lineEdit().setText(str(self.energy))

        self.tableLabel = QLabel(r"&data table:")
        self.tableCB = QComboBox()
        self.tableCB.addItems(tables)
        self.tableCB.setCurrentIndex(3)  # Chantler total
        self.tableCB.currentIndexChanged.connect(self.calculate)
        self.tableLabel.setBuddy(self.tableCB)
        self.tablePlotButton = QPushButton('&Plot f"')
        if not MAC:
            self.tablePlotButton.setFocusPolicy(QtCore.Qt.NoFocus)

        self.tablePlotButton.clicked.connect(self.plotf)

        self.resNuLabel = QLabel("")
        self.resNu = QLabel("")
        self.resMassLabel = QLabel("")
        self.resMass = QLabel("")

        self.stepLabel = QLabel("absorptance step = ")
        self.stepCB = QComboBox()
        self.stepCB.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
#        self.stepCB.setMinimumWidth(200)

        self.extraRhoLabel = QLabel(u"ρ (g/cm<sup>3</sup>) = ")
        self.extraRhoEdit = QLineEdit()
#        self.extraRhoEdit.setMinimumSize(fm.width("8.88"), fm.height())
        self.extraRhoEdit.setSizePolicy(
            QSizePolicy.MinimumExpanding, QSizePolicy.Fixed)
        self.extraDLabel = QLabel(u"d (µm) = ")
        self.extraD = QLabel(u"")

        self.buttonCalculate = QPushButton("Calculate")
        self.buttonCalculate.clicked.connect(self.calculate)
        self.buttonAbout = QPushButton(r"&About...")
        self.buttonAbout.clicked.connect(self.about)
        self.buttonHelp = QPushButton(r"&Help...")
        self.buttonHelp.clicked.connect(self.myhelp)

        whatLayout = QHBoxLayout()
        whatLayout.addWidget(self.whatLabel)
        whatLayout.addWidget(self.whatCB)
        whatLayout.addStretch()

        compoundLayout = QHBoxLayout()
        compoundLayout.addWidget(self.compoundLabel)
        compoundLayout.addWidget(self.compoundExLabel)
        compoundLayout.addStretch()

        compoundMassLayout = QHBoxLayout()
        compoundMassLayout.addWidget(self.compoundMassLabel)
        compoundMassLayout.addWidget(self.compoundMass)
        compoundMassLayout.addStretch()

        muTdLayout = QHBoxLayout()
        muTdLayout.addWidget(self.muTdLabel)
        muTdLayout.addWidget(self.muTdEdit)
        muTdLayout.addStretch()
        muTdLayout.addWidget(self.areaLabel)
        muTdLayout.addWidget(self.areaEdit)
        muTdLayout.addStretch()
        muTdLayout.addWidget(self.dmudLabel)
        muTdLayout.addWidget(self.dmudEdit)
        muTdLayout.addStretch()

        energyLayout = QHBoxLayout()
        energyLayout.addWidget(self.energyLabel)
        energyLayout.addWidget(self.energyCB)
        energyLayout.addStretch()

        tableLayout = QHBoxLayout()
        tableLayout.addWidget(self.tableLabel)
        tableLayout.addWidget(self.tableCB)
        tableLayout.addStretch()
        tableLayout.addWidget(self.tablePlotButton)
        tableLayout.addStretch()

        resLayout = QHBoxLayout()
        resLayout.addWidget(self.resNuLabel)
        resLayout.addWidget(self.resNu)
        resLayout.addStretch()
        resLayout.addWidget(self.resMassLabel)
        resLayout.addWidget(self.resMass)
        resLayout.addStretch()

        stepLayout = QHBoxLayout()
        stepLayout.addWidget(self.stepLabel)
        stepLayout.addWidget(self.stepCB)
#        stepLayout.addStretch()

        extraLayout = QHBoxLayout()
        extraLayout.addWidget(self.extraRhoLabel)
        extraLayout.addWidget(self.extraRhoEdit)
        extraLayout.addStretch()
        extraLayout.addWidget(self.extraDLabel)
        extraLayout.addWidget(self.extraD)
        extraLayout.addStretch()

        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(self.buttonCalculate)
        buttonLayout.addWidget(self.buttonAbout)
        buttonLayout.addWidget(self.buttonHelp)
        buttonLayout.addStretch()

        layout = QVBoxLayout()
        layout.addLayout(whatLayout)
        layout.addWidget(self.formula)
        hline1 = QFrame(self)
        hline1.setFrameStyle(QFrame.HLine)
        layout.addWidget(hline1)
        layout.addLayout(compoundLayout)
        layout.addWidget(self.compoundEdit)
        layout.addLayout(compoundMassLayout)
        hline2 = QFrame(self)
        hline2.setFrameStyle(QFrame.HLine)
        layout.addWidget(hline2)
        layout.addLayout(muTdLayout)
        layout.addLayout(energyLayout)
        hline3 = QFrame(self)
        hline3.setFrameStyle(QFrame.HLine)
        layout.addWidget(hline3)
        layout.addLayout(tableLayout)
        hline4 = QFrame(self)
        hline4.setFrameStyle(QFrame.HLine)
        layout.addWidget(hline4)
        layout.addLayout(resLayout)
        layout.addLayout(stepLayout)
        layout.addLayout(extraLayout)
        layout.addStretch()
        layout.addLayout(buttonLayout)
        self.setLayout(layout)

        self.setWindowTitle("XAFSmassQt")
        selfDir = os.path.dirname(__file__)
        self.setWindowIcon(QtGui.QIcon(
            os.path.join(selfDir, 'help', '_static', 'XAFSmassQt.ico')))
        self.setWindowFlags(QtCore.Qt.Window)
        self.move(QApplication.desktop().screen().rect().center() -
                  self.rect().center())
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
            self.resNuLabel.show()
            self.resNu.show()
        elif self.what == FOIL:
            self.areaLabel.setText(u"ρ (g/cm<sup>3</sup>) = ")
            self.resMassLabel.setText(u"d (µm) = ")
            self.resNuLabel.hide()
            self.resNu.hide()
        elif self.what == GAS:
            self.areaLabel.setText(u"d (cm) = ")
            self.resMassLabel.setText(u"p (mbar) = ")
            self.resNuLabel.hide()
            self.resNu.hide()
        elif self.what == XCONTENT:
            self.areaLabel.setText(u"Δµd = ")
            self.resNuLabel.setText("N<sub>x</sub> = ")
            self.resMassLabel.setText(r"wt%<sub>x</sub> = ")
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

        self.calculate()

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

    def energyCB_selected(self):
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
                    cMass += xc.read_atomic_data(e[0])*e[1]
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
            return
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
                QMessageBox.critical(self, "Error", res)
                return
            nu, m, th, eDict = res
            self.resNu.setText('{0}'.format(xc.round_to_n(nu, 3)))
            self.resMass.setText("<b>{0}</b>".format(xc.round_to_n(m, 3)))
            if rho is not None:
                self.extraD.setText('{0}'.format(xc.round_to_n(th, 4)))
            iSelect = 0
            for i, (k, v) in enumerate(eDict.items()):
                self.stepCB.addItem("{0}({1}mg={2}wt%): {3}".format(
                    k, xc.round_to_n(v[-2], 3), xc.round_to_n(v[-2]/m*100., 3),
                    xc.round_to_n(v[-1], 3)))
                if v[2] > 0:
                    iSelect = i
            self.stepCB.setCurrentIndex(iSelect)
        elif self.what == FOIL:
            res = xc.calculate_foil(
                self.parsedCompound.asList(), E, muTd, rho, table=table)
            if isinstance(res, str):
                self.clear_results()
                QMessageBox.critical(self, "Error", res)
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
                QMessageBox.critical(self, "Error", res)
                return
            P = res
            self.resMass.setText("<b>{0}</b>".format(xc.round_to_n(P, 3)))
        elif self.what == XCONTENT:
            res = xc.calculate_x(
                self.parsedCompound.asList(), E, muTd, Dmud, dmud, table=table)
            if isinstance(res, str):
                self.clear_results()
                QMessageBox.critical(self, "Error", res)
                return
            nu, wt, wt1 = res
            self.resNu.setText('{0}'.format(xc.round_to_n(nu, 3)))
            outStr = '{0}'.format(xc.round_to_n(wt, 5))
            if wt1 > 0:
                outStr += u'±{0}'.format(xc.round_to_n(wt1, 3))
            self.resMass.setText(u"<b>{0}</b>".format(outStr))

    def about(self):
        import platform
        if use_pyside:
            Qt_version = QtCore.__version__
            PyQt_version = PySide.__version__
        else:
            Qt_version = QtCore.QT_VERSION_STR
            PyQt_version = QtCore.PYQT_VERSION_STR
        QMessageBox.about(
            self, "About XAFSmassQt",
            """<b>XAFSmass(Qt)</b> v {0}
            <ul>
            <li>{1[0]}
            <li>{1[1]}
            </ul>
            <p>Open source, {2}. Available at PyPI and GitHub<p>
            <p>Your system:
            <ul>
            <li>{3}
            <li>Python {4}
            <li>Qt {5}
            <li>{6} {7}
            </ul>""".format(
                __version__, __author__.split(','), __license__,
                platform.platform(terse=True), platform.python_version(),
                Qt_version, QtName, PyQt_version))

    def myhelp(self):
        hname = os.path.join(os.path.dirname(__file__), 'help', 'index.html')
        webbrowser.open(hname)

    def plotf(self):
        self.parse_compound()
        E = float(self.energyCB.lineEdit().text())
        table = tablesF[self.tableCB.currentIndex()]
        if self.plotDlg is None:
            self.plotDlg = PlotDlg(self, self.parsedCompound, E, table)
        else:
            self.plotDlg.plotCanvas.plot(self.parsedCompound, E, table)


def run():
    app = QApplication(sys.argv)
    form = MainDlg()
    form.show()
    app.exec_()


if __name__ == "__main__":
    run()
