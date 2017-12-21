#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import threading
import time
import numpy as np
from PySide import QtGui, QtCore
from . import utils
from . import input
from . import _ember

import matplotlib
matplotlib.rcParams['backend.qt4'] = 'PySide'
matplotlib.use('Qt4Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

if sys.version_info.major == 3:
    _stringTypes = (str,)
else:
    _stringTypes = (str, unicode)

class SolverThread(threading.Thread):
    def __init__(self, *args, **kwargs):
        threading.Thread.__init__(self)
        self.solver = kwargs['solver']
        self.conf = kwargs['conf'].evaluate()
        self.solver.lock = threading.Lock()
        self.solver.progress = 0.0
        self.__stop = False
        self.daemon = True
        self.stage = 0

        self.frac0 = 0.15

    def run(self):
        done = 0
        i = 0
        while not self.__stop and not done:
            i += 1
            with self.solver.lock:
                done = self.solver.step()
            time.sleep(1e-3)
            if not i % 5:
                self.updateProgress()
            if done:
                self.solver.progress = 1.0
        self.__stop = True

    def stop(self):
        self.__stop = True

    def updateProgress(self):
        TC = self.conf.terminationCondition
        errNow = self.solver.terminationCondition
        if TC.measurement is None:
            self.solver.progress = self.solver.timeseriesWriter.t[-1] / TC.tEnd
        elif self.stage == 0:
            # First part: getting past the minimum steady-state measurement period
            self.solver.progress = self.frac0 * self.solver.timeseriesWriter.t[-1] / TC.steadyPeriod
            if errNow < 1e9:
                self.refCond = errNow
                self.stage = 1
        else:
            # Second part: vaguely linearizing the approach to steady-state
            A = (np.log10(TC.tolerance + (errNow-TC.tolerance)/self.refCond) /
                 np.log10(TC.tolerance))
            P = min(self.frac0 + (1-self.frac0) * A ** 0.5 , 1.0)
            self.solver.progress = max(P, self.solver.progress) # never go backwards

class OptionWidget(QtGui.QWidget):
    def __init__(self, label, opt, *args, **kwargs):
        QtGui.QWidget.__init__(self)

        self.opt = opt
        self.label = label
        self.optName = self.label.text()
        self.setLayout(QtGui.QHBoxLayout())
        self.layout().setContentsMargins(0,0,0,0)

    def updateOpt(self):
        if self.opt.value != self.opt.default:
            self.label.setText('<i>%s</i>' %  self.optName)
        else:
            self.label.setText('%s' %  self.optName)
        self.parent().parent().updateVisibility()


class StringOptionWidget(OptionWidget):
    def __init__(self, label, opt, *args, **kwargs):
        OptionWidget.__init__(self, label, opt)
        self.text = QtGui.QLineEdit(opt.value)
        self.layout().addWidget(self.text)
        self.text.textChanged.connect(self.updateOpt)

    def updateOpt(self):
        self.opt.value = str(self.text.text())
        OptionWidget.updateOpt(self)


class NumericOptionWidget(OptionWidget):
    def __init__(self, label, opt, *args, **kwargs):
        OptionWidget.__init__(self, label, opt)
        self.text = QtGui.QLineEdit(str(opt.value))
        self.layout().addWidget(self.text)
        self.text.textChanged.connect(self.updateOpt)

    def updateOpt(self):
        try:
            self.opt.value = float(self.text.text())
        except ValueError:
            pass
        OptionWidget.updateOpt(self)


class IntegerOptionWidget(NumericOptionWidget):
    def updateOpt(self):
        try:
            self.opt.value = int(self.text.text())
        except ValueError:
            pass
        OptionWidget.updateOpt(self)


class BoolOptionWidget(OptionWidget):
    def __init__(self, label, opt, *args, **kwargs):
        OptionWidget.__init__(self, label, opt)

        self.trueWidget = QtGui.QRadioButton('True')
        self.falseWidget = QtGui.QRadioButton('False')
        self.noneWidget = QtGui.QRadioButton('None')
        if opt.value:
            self.trueWidget.toggle()
        else:
            self.falseWidget.toggle()
        self.layout().addWidget(self.trueWidget)
        self.layout().addWidget(self.falseWidget)
        self.layout().addWidget(self.noneWidget)
        self.noneWidget.hide()

        self.trueWidget.released.connect(self.updateOpt)
        self.falseWidget.released.connect(self.updateOpt)
        self.layout().addStretch(1)

        self._savedValue = None

    def updateOpt(self):
        self.opt.value = self.trueWidget.isChecked()
        OptionWidget.updateOpt(self)

    def setEnabled(self, tf):
        if tf and not self.isEnabled() and self._savedValue is not None:
            self.trueWidget.setChecked(self._savedValue)
            self.falseWidget.setChecked(not self._savedValue)
            self._savedValue = None
        elif not tf and self.isEnabled() and self._savedValue is None:
            self._savedValue = self.trueWidget.isChecked()
            self.noneWidget.setChecked(True)

        OptionWidget.setEnabled(self, tf)


class EnumOptionWidget(OptionWidget):
    def __init__(self, label, opt, *args, **kwargs):
        OptionWidget.__init__(self, label, opt)

        self.combo = QtGui.QComboBox()
        self.items = {}
        for i,choice in enumerate(opt.choices):
            if choice == opt.value:
                startIndex = i

            self.combo.addItem(str(choice))
            self.items[i] = choice

        self.combo.setCurrentIndex(startIndex)
        self.layout().addWidget(self.combo)
        self.combo.currentIndexChanged.connect(self.updateOpt)

    def updateOpt(self):
        self.opt.value = self.items[self.combo.currentIndex()]
        OptionWidget.updateOpt(self)


class OptionsWidget(QtGui.QGroupBox):
    def __init__(self, opts, *args, **kwargs):
        QtGui.QGroupBox.__init__(self)
        self.opts = opts
        self.setLayout(QtGui.QGridLayout())
        self.layout().setSpacing(0)
        self.setTitle(self.opts.__class__.__name__)
        self.optionWidgets = []

        width = 0
        for i,(name,opt) in enumerate(self.opts):
            if opt.label:
                label = QtGui.QLabel(opt.label)
                label.setToolTip('<tt>%s</tt>' % name)
            else:
                label = QtGui.QLabel(name)
            self.layout().addWidget(label, i, 0)
            width = max(label.sizeHint().width(), width)

            if opt.choices is not None:
                w = EnumOptionWidget(label, opt)
            elif isinstance(opt, input.StringOption):
                w = StringOptionWidget(label, opt)
            elif isinstance(opt, input.IntegerOption):
                w = IntegerOptionWidget(label, opt)
            elif isinstance(opt, input.FloatOption):
                w = NumericOptionWidget(label, opt)
            elif isinstance(opt, input.BoolOption):
                w = BoolOptionWidget(label, opt)
            else:
                w = QtGui.QLabel(str(opt.value))
                w.opt = opt

            self.layout().addWidget(w, i, 1)
            self.optionWidgets.append((label,w))

        self.layout().setVerticalSpacing(4)
        self.layout().setColumnMinimumWidth(0, width + 5)
        spacer = QtGui.QSpacerItem(1, 1000, QtGui.QSizePolicy.Minimum,
                                   QtGui.QSizePolicy.Maximum)
        self.layout().addItem(spacer, i+1, 0)

    def updateVisibility(self, level, conf):
        anyVisible = False
        anyEnabled = False
        for label,w in self.optionWidgets:
            if w.opt.level > level:
                w.hide()
                label.hide()
            else:
                anyVisible = True
                w.show()
                label.show()
                e = w.opt.shouldBeEnabled(conf)
                w.setEnabled(e)
                label.setEnabled(e)
                if e:
                    anyEnabled = True
        return anyVisible, anyEnabled


class MultiOptionsWidget(QtGui.QWidget):
    """ Widget used for presenting solver configuration options """
    def __init__(self, conf, *args, **kwargs):
        QtGui.QWidget.__init__(self)
        self.conf = conf
        self.setLayout(QtGui.QHBoxLayout())
        self.optionsList = QtGui.QListWidget()
        self.optionsList.setSpacing(1)
        self.layout().addWidget(self.optionsList)
        self.activeOptionWidget = None
        self.optionWidgets = []
        self.level = 0

        height = 0
        for item in self.conf:
            listitem = QtGui.QListWidgetItem(item.__class__.__name__)
            self.optionsList.addItem(listitem)
            w = OptionsWidget(item)
            self.optionWidgets.append((listitem,w))
            height = max(height, w.minimumSizeHint().height())
            listitem.widget = w
            self.layout().addWidget(w)
            w.hide()

        self.optionsList.setCurrentRow(0)
        self.setActiveWidget(self.optionsList.currentItem())
        self.setMinimumHeight(height)
        width = self.optionsList.sizeHintForColumn(0) + 10
        self.optionsList.setMinimumWidth(width)
        self.optionsList.setMaximumWidth(width)
        self.optionsList.currentItemChanged.connect(self.setActiveWidget)
        self.optionsList.setSizePolicy(QtGui.QSizePolicy.Fixed,
                                       QtGui.QSizePolicy.Preferred)

    def setActiveWidget(self, listitem):
        if self.activeOptionWidget is not None:
            self.activeOptionWidget.hide()

        self.activeOptionWidget = listitem.widget
        self.activeOptionWidget.show()

    def updateVisibility(self, level=None):
        if level is not None:
            self.level = level
        for listitem, w in self.optionWidgets:
            visible, enabled = w.updateVisibility(self.level, self.conf)
            listitem.setHidden(not visible or not enabled)


class SolverWidget(QtGui.QWidget):
    """ Widget used to run and monitor the Ember solver """
    def __init__(self, conf, *args, **kwargs):
        QtGui.QWidget.__init__(self)

        self.conf = conf
        self.setLayout(QtGui.QVBoxLayout())

        # Buttons
        self.startButton = QtGui.QPushButton('Start')
        self.stopButton = QtGui.QPushButton('Stop')
        self.resetButton = QtGui.QPushButton('Reset')
        self.buttons = QtGui.QWidget()
        self.buttons.setLayout(QtGui.QHBoxLayout())
        self.buttons.layout().addWidget(self.startButton)
        self.buttons.layout().addWidget(self.stopButton)
        self.buttons.layout().addWidget(self.resetButton)
        self.layout().addWidget(self.buttons)

        self.startButton.pressed.connect(self.run)
        self.stopButton.pressed.connect(self.stop)
        self.resetButton.pressed.connect(self.reset)

        # Progress Bar
        self.progressBar = QtGui.QProgressBar()
        self.layout().addWidget(self.progressBar)
        self.progressBar.setRange(0, 1000)
        self.progressBar.setValue(0)

        # Graphs
        self.graphContainer = QtGui.QWidget()
        self.graphContainer.setLayout(QtGui.QHBoxLayout())
        self.layout().addWidget(self.graphContainer)

        self.fig = Figure(figsize=(600,400), dpi=72)
        self.fig.subplots_adjust(0.09, 0.08, 0.93, 0.96, wspace=0.3)
        self.ax1 = self.fig.add_subplot(1,2,1)
        self.ax1.set_xlabel('time [ms]')
        self.ax1.set_ylabel('Consumption Speed, $S_c$ [cm/s]')
        self.Sc_timeseries = self.ax1.plot([0],[0], lw=2)[0]

        self.ax2a = self.fig.add_subplot(1,2,2)
        self.ax2b = self.ax2a.twinx()
        self.ax2a.set_xlabel('flame coordinate [mm]')
        self.ax2a.set_ylabel('Temperature [K]')
        self.ax2b.set_ylabel('Heat Release Rate [MW/m$^3$]')

        self.T_profile = self.ax2a.plot([0],[0], 'b', lw=2)[0]
        self.hrr_profile = self.ax2b.plot([0],[0], 'r', lw=2)[0]

        self.canvas = FigureCanvas(self.fig)
        self.graphContainer.layout().addWidget(self.canvas)
        bgcolor = self.palette().color(QtGui.QPalette.Window)
        self.fig.set_facecolor((bgcolor.redF(), bgcolor.greenF(), bgcolor.blueF()))
        #self.fig.patch.set_alpha(0.1)

        # internals
        self.solver = None
        self.solverThread = None
        self.updateTimer = QtCore.QTimer()
        self.updateTimer.setInterval(0.5)
        self.updateTimer.timeout.connect(self.updateStatus)
        self.running = False
        self.updateButtons()

    def run(self):
        if self.solverThread is not None and self.solverThread.is_alive():
            return

        if self.solver is None:
            self.conf.validate()
            self.solver = _ember.FlameSolver(self.conf)
            self.solver.initialize()

        self.solverThread = SolverThread(solver=self.solver,
                                         conf=self.conf)
        self.solverThread.start()
        self.updateTimer.start()
        self.running = True
        self.updateButtons()

    def stop(self):
        if self.solverThread:
            self.solverThread.stop()
            self.updateTimer.stop()
        self.running = False
        self.startButton.setText('Resume')
        self.updateButtons()

    def reset(self):
        self.progressBar.setValue(0)
        self.T_profile.set_data([0], [0])
        self.hrr_profile.set_data([0], [0])
        self.Sc_timeseries.set_data([0], [0])
        self.canvas.draw()
        self.solver = None
        self.startButton.setText('Start')
        self.updateButtons()

    def updateButtons(self):
        running = self.running and self.solverThread is not None and self.solverThread.is_alive()
        self.startButton.setEnabled(not running)
        self.stopButton.setEnabled(running)
        self.resetButton.setEnabled(self.solver is not None and not running)

    def updateStatus(self):
        if not self.solver:
            return

        if not self.solverThread.is_alive():
            self.running = False
            self.updateTimer.stop()
            self.updateButtons()

        if self.solver.progress > 0:
            self.progressBar.setValue(1000 * self.solver.progress)
        with self.solver.lock:
            t = np.array(self.solver.timeseriesWriter.t)
            Sc = np.array(self.solver.timeseriesWriter.Sc)

            self.T_profile.set_data(self.solver.x * 1000,
                                    self.solver.T)
            self.hrr_profile.set_data(self.solver.x * 1000,
                                      self.solver.qDot / 1e6)
        self.Sc_timeseries.set_data(1000 * t, Sc * 100)

        for ax in (self.ax1, self.ax2a, self.ax2b):
            ax.relim()
            ax.autoscale_view(False, True, True)
        self.canvas.draw()


class MainWindow(QtGui.QMainWindow):
    def __init__(self, *args, **kwargs):
        QtGui.QMainWindow.__init__(self)

        w = QtGui.QWidget()
        self.setCentralWidget(w)
        self.resize(800,600)
        self.setWindowTitle('Simple')

        fileMenu = self.menuBar().addMenu('&File')
        optMenu = self.menuBar().addMenu('&Options')

        self.addToMenu(fileMenu, '&New', lambda: self.new())
        self.addToMenu(fileMenu, '&Open...', lambda: self.openConf())
        self.addToMenu(fileMenu, '&Save', lambda: self.saveConf(True))
        self.addToMenu(fileMenu, 'Save &as...', lambda: self.saveConf(False))
        self.addToMenu(fileMenu, '&Quit', self.close)

        optLevelGroup = QtGui.QActionGroup(optMenu)
        a = self.addToMenu(optMenu, '&Basic',
                           lambda: self.setLevel(0), optLevelGroup)
        self.addToMenu(optMenu, '&Advanced',
                       lambda: self.setLevel(1), optLevelGroup)
        self.addToMenu(optMenu, '&Expert',
                       lambda: self.setLevel(2), optLevelGroup)
        a.setChecked(True)

        self.level = 0
        self.confFileName = None
        if len(args) == 2:
            self.new(args[1])
        else:
            self.new()

    def addToMenu(self, menu, name, triggerFunc, group=None):
        a = QtGui.QAction(name, self)
        a.triggered.connect(triggerFunc)
        menu.addAction(a)
        if group:
            a.setCheckable(True)
            group.addAction(a)
        return a

    def setLevel(self, level):
        self.level = level
        self.confWidget.updateVisibility(level)

    def new(self, conf=None):
        if conf is None:
            self.conf = input.Config(input.Paths(logFile='gui-runlog.txt'))
        elif isinstance(conf, input.Config):
            self.conf = conf
        elif isinstance(conf, _stringTypes):
            if not os.path.exists(conf):
                print("Can't find input file '%s'" % conf)
                return
            localenv = {}
            execstatements = ['from numpy import *',
                              'import numpy as np',
                              'from ember.input import *']

            execstatements.extend(open(conf).readlines())
            exec('\n'.join(execstatements), localenv)
            self.conf = localenv['conf']
            self.confFileName = conf

        self.tabWidget = QtGui.QTabWidget()
        self.setCentralWidget(self.tabWidget)

        self.confWidget = MultiOptionsWidget(self.conf)
        self.tabWidget.addTab(self.confWidget, 'Configure')
        self.setLevel(self.level)

        self.runWidget = SolverWidget(self.conf)
        self.tabWidget.addTab(self.runWidget, 'Run')

        self.tabWidget.addTab(QtGui.QWidget(), 'Analyze') #TODO: unimplemented

    def openConf(self):
        fileinfo = QtGui.QFileDialog.getOpenFileName(
            self, 'Select Configuration', '.', 'Flame Configurations (*.py *.conf)')

        # Dealing with an incompatibility between PySide and PyQt
        filename = str(fileinfo[0] if isinstance(fileinfo, tuple) else fileinfo)

        if os.path.exists(filename):
            self.new(filename)

    def saveConf(self, useExisting):
        if not useExisting or self.confFileName is None:
            fileinfo = QtGui.QFileDialog.getSaveFileName(
                        self, 'Select Configuration', '.', 'Flame Configurations (*.py *.conf)')

            # Dealing with an incompatibility between PySide and PyQt
            filename = str(fileinfo[0] if isinstance(fileinfo, tuple) else fileinfo)

            if not filename:
                return

            if not filename.endswith('.py') and not filename.endswith('.conf'):
                filename += '.py'

            # Confirm before overwriting an existing file
            if os.path.exists(filename) and not useExisting:
                dlg = QtGui.QMessageBox(self.parent())
                dlg.setText("A file named '%s' already exists." % filename)
                dlg.setInformativeText("Do you wish to overwrite it?")
                dlg.setStandardButtons(dlg.Yes | dlg.No)
                dlg.setDefaultButton(dlg.Yes)

                ret = dlg.exec_()
                if ret == dlg.No:
                    self.saveConf(False)
                    return
                elif ret != dlg.Yes:
                    print('unknown return value:', ret)

            self.confFileName = filename
        else:
            filename = self.confFileName

        outFile = open(filename, 'w')
        outFile.write(self.conf.stringify())


def main():
    app = QtGui.QApplication(sys.argv)
    app.setStyle("Plastique")
    window = MainWindow(*sys.argv)
    window.show()

    sys.exit(app.exec_())
