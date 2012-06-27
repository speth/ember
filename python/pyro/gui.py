#!/usr/bin/env python

import sys
from PySide import QtGui
import utils
import input

class OptionWidget(QtGui.QWidget):
    def __init__(self, label, opt, *args, **kwargs):
        QtGui.QWidget.__init__(self)

        self.opt = opt
        self.label = label
        self.optName = self.label.text()
        self.setLayout(QtGui.QHBoxLayout())
        self.layout().setContentsMargins(0,0,0,0)

    def checkDefault(self):
        if self.opt.value != self.opt.default:
            self.label.setText('<i>%s</i>' %  self.optName)
        else:
            self.label.setText('%s' %  self.optName)


class StringOptionWidget(OptionWidget):
    def __init__(self, label, opt, *args, **kwargs):
        OptionWidget.__init__(self, label, opt)
        self.text = QtGui.QLineEdit(opt.value)
        self.layout().addWidget(self.text)
        self.text.textChanged.connect(self.updateOpt)

    def updateOpt(self):
        self.opt.value = self.text.text()
        self.checkDefault()


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
        self.checkDefault()


class BoolOptionWidget(OptionWidget):
    def __init__(self, label, opt, *args, **kwargs):
        OptionWidget.__init__(self, label, opt)

        self.trueWidget = QtGui.QRadioButton('True')
        self.falseWidget = QtGui.QRadioButton('False')
        if opt.value:
            self.trueWidget.toggle()
        else:
            self.falseWidget.toggle()
        self.layout().addWidget(self.trueWidget)
        self.layout().addWidget(self.falseWidget)

        self.trueWidget.toggled.connect(self.updateOpt)
        self.layout().addStretch(1)

    def updateOpt(self):
        self.opt.value = self.trueWidget.isChecked()
        self.checkDefault()


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
        self.checkDefault()


class OptionsWidget(QtGui.QGroupBox):
    def __init__(self, opts, *args, **kwargs):
        QtGui.QGroupBox.__init__(self)
        self.opts = opts
        self.setLayout(QtGui.QGridLayout())
        self.layout().setSpacing(0)
        self.setTitle(self.opts.__class__.__name__)

        width = 0
        for i,(name,opt) in enumerate(self.opts):
            label = QtGui.QLabel(name)
            self.layout().addWidget(label, i, 0)
            width = max(label.sizeHint().width(), width)

            if opt.choices is not None:
                w = EnumOptionWidget(label, opt)
            elif isinstance(opt, input.StringOption):
                w = StringOptionWidget(label, opt)
            elif isinstance(opt, (input.FloatOption, input.IntegerOption)):
                w = NumericOptionWidget(label, opt)
            elif isinstance(opt, input.BoolOption):
                w = BoolOptionWidget(label, opt)
            else:
                w = QtGui.QLabel(str(opt.value))

            self.layout().addWidget(w, i, 1)

        self.layout().setVerticalSpacing(4)
        self.layout().setColumnMinimumWidth(0, width + 5)
        spacer = QtGui.QSpacerItem(1, 1000, QtGui.QSizePolicy.Minimum,
                                   QtGui.QSizePolicy.Maximum)
        self.layout().addItem(spacer, i+1, 0)


class MultiOptionsWidget(QtGui.QWidget):
    def __init__(self, conf, *args, **kwargs):
        QtGui.QWidget.__init__(self)
        self.conf = conf
        self.setLayout(QtGui.QHBoxLayout())
        self.optionsList = QtGui.QListWidget()
        self.optionsList.setSpacing(1)
        self.layout().addWidget(self.optionsList)
        self.activeOptionWidget = None

        height = 0
        for item in self.conf:
            listitem = QtGui.QListWidgetItem(item.__class__.__name__)
            self.optionsList.addItem(listitem)
            w = OptionsWidget(item)
            height = max(height, w.minimumSizeHint().height())
            listitem.widget = w
            self.layout().addWidget(w)
            w.hide()

        self.optionsList.setCurrentRow(0)
        self.setActiveWidget(self.optionsList.currentItem())
        self.setMinimumHeight(height)
        width = self.optionsList.sizeHintForColumn(0) + 5
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


class MainWindow(QtGui.QMainWindow):
    def __init__(self, *args, **kwargs):
        QtGui.QMainWindow.__init__(self)

        w = QtGui.QWidget()
        self.setCentralWidget(w)
        self.resize(800,600)
        self.setWindowTitle('Simple')

        fileMenu = self.menuBar().addMenu('&File')
        a = QtGui.QAction('&New', self)
        a.triggered.connect(lambda: self.newConf())
        fileMenu.addAction(a)

        a = QtGui.QAction('&Quit', self)
        a.triggered.connect(self.close)
        fileMenu.addAction(a)

        self.newConf()

    def newConf(self):
        self.conf = input.Config()
        w = MultiOptionsWidget(self.conf)
        self.setCentralWidget(w)


def main():
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
