# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'diagram.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from CoolProp.Plots import PropertyPlot
import CoolProp

class Ui_Dialog(QtWidgets.QDialog):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(769, 670)
        self.progressBar = QtWidgets.QProgressBar(Dialog)
        self.progressBar.setGeometry(QtCore.QRect(20, 630, 731, 23))
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName("progressBar")
        self.graphicsView = QtWidgets.QGraphicsView(Dialog)
        self.graphicsView.setGeometry(QtCore.QRect(20, 80, 731, 531))
        self.graphicsView.setObjectName("graphicsView")
        self.splitter = QtWidgets.QSplitter(Dialog)
        self.splitter.setGeometry(QtCore.QRect(620, 40, 131, 31))
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.pushButton_draw = QtWidgets.QPushButton(self.splitter)
        self.pushButton_draw.setObjectName("pushButton_draw")
        self.pushButton_clear = QtWidgets.QPushButton(self.splitter)
        self.pushButton_clear.setObjectName("pushButton_clear")
        self.splitter_2 = QtWidgets.QSplitter(Dialog)
        self.splitter_2.setGeometry(QtCore.QRect(20, 50, 231, 21))
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.comboBox = QtWidgets.QComboBox(self.splitter_2)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox_2 = QtWidgets.QComboBox(self.splitter_2)
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        self.pushButton_draw.clicked.connect(self.draw)
        self.pushButton_clear.clicked.connect(self.clear)
    def draw(self):
        type = self.comboBox_2.currentText()
        fluid = self.comboBox.currentText()
        if type == "P-h diagram":
            plot = PropertyPlot(fluid, 'ph')
            plot.calc_isolines()
            plot.show()
            self.graphicsView = plot
            self.graphicsView.show()
            
        elif type == "T-s diagram":
            plot = PropertyPlot(fluid, 'Ts', tp_limits='ORC')
            plot.calc_isolines(CoolProp.iQ, num=6)
            self.graphicsView = plot
            self.graphicsView.show()
        
    def clear(self):
        self.graphicsView.clear()


        
    
        
        

        
    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.pushButton_draw.setText(_translate("Dialog", "Draw"))
        self.pushButton_clear.setText(_translate("Dialog", "Clear"))
        self.comboBox.setItemText(0, _translate("Dialog", "air"))
        self.comboBox.setItemText(1, _translate("Dialog", "water"))
        self.comboBox.setItemText(2, _translate("Dialog", "R22"))
        self.comboBox.setItemText(3, _translate("Dialog", "R209"))
        self.comboBox.setItemText(4, _translate("Dialog", "R260"))
        self.comboBox.setItemText(5, _translate("Dialog", "직접입력"))
        self.comboBox_2.setItemText(0, _translate("Dialog", "P-h diagram"))
        self.comboBox_2.setItemText(1, _translate("Dialog", "T-s diagram"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
