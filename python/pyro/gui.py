#!/usr/bin/env python

import sys
from PySide import QtGui
import utils
import input


class OptionWidget(QtGui.QWidget):
    def __init__(self, name, opt, *args, **kwargs):
        QtGui.QWidget.__init__(self)

        self.optName = name
        self.opt = opt

        self.setLayout(QtGui.QHBoxLayout())
        self.label = QtGui.QLabel(name)
        self.layout().addWidget(self.label)

    def checkDefault(self):
        if self.opt.value != self.opt.default:
            self.label.setText('<i>%s</i>' %  self.optName)
        else:
            self.label.setText('%s' %  self.optName)


class StringOptionWidget(OptionWidget):
    def __init__(self, name, opt, *args, **kwargs):
        OptionWidget.__init__(self, name, opt)
        self.text = QtGui.QLineEdit(opt.value)
        self.layout().addWidget(self.text)
        self.text.textChanged.connect(self.updateOpt)

    def updateOpt(self):
        self.opt.value = self.text.text()
        self.checkDefault()


class BoolOptionWidget(OptionWidget):
    def __init__(self, name, opt, *args, **kwargs):
        OptionWidget.__init__(self, name, opt)

        self.trueWidget = QtGui.QRadioButton('True')
        self.falseWidget = QtGui.QRadioButton('False')
        if opt.value:
            self.trueWidget.toggle()
        else:
            self.falseWidget.toggle()
        self.layout().addWidget(self.trueWidget)
        self.layout().addWidget(self.falseWidget)

        self.trueWidget.toggled.connect(self.updateOpt)

    def updateOpt(self):
        self.opt.value = self.trueWidget.isChecked()
        self.checkDefault()


class EnumOptionWidget(OptionWidget):
    def __init__(self, name, opt, *args, **kwargs):
        OptionWidget.__init__(self, name, opt)

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
        self.setLayout(QtGui.QVBoxLayout())
        self.setTitle(self.opts.__class__.__name__)
        for name,opt in self.opts:
            if opt.choices is not None:
                w = EnumOptionWidget(name, opt)
            elif isinstance(opt, input.StringOption):
                w = StringOptionWidget(name, opt)
            elif isinstance(opt, input.BoolOption):
                w = BoolOptionWidget(name, opt)
            else:
                w = QtGui.QLabel(name + ' = ' + str(opt.value))
            self.layout().addWidget(w)


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
        w = QtGui.QWidget()
        w.setLayout(QtGui.QVBoxLayout())
        for item in self.conf:
            w2 = OptionsWidget(item)
            w.layout().addWidget(w2)

        scroll = QtGui.QScrollArea()
        scroll.setWidget(w)
        self.setCentralWidget(scroll)


def main():
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
