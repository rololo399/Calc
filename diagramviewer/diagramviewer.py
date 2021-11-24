# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'diagramviewer.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets
import sys
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import time
import resources_qr

class Ui_Form(QWidget):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(561, 517)
        Form.setStyleSheet("background-image: url(:/background/배경.png);")
        self.gridLayout_2 = QtWidgets.QGridLayout(Form)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.comboBox = QtWidgets.QComboBox(Form)
        self.comboBox.setStyleSheet("color: qlineargradient(spread:pad, x1:0, y1:1, x2:1, y2:0, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(255, 255, 255, 255));\n"
"background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(69, 83, 100, 255), stop:1 rgba(255, 255, 255, 255));\n"
"background-color: qlineargradient(spread:pad, x1:1, y1:0, x2:1, y2:0, stop:0 rgba(69, 83, 100, 255), stop:1 rgba(255, 255, 255, 255));\n"
"color: rgb(255, 255, 255);\n"
"background-image: url(:/background/버튼.png);\n"
"font: 10pt \"맑은 고딕\";")
        self.comboBox.setObjectName("comboBox")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.gridLayout.addWidget(self.comboBox, 0, 0, 1, 1)
        self.comboBox_2 = QtWidgets.QComboBox(Form)
        self.comboBox_2.setStyleSheet("color: qlineargradient(spread:pad, x1:0, y1:1, x2:1, y2:0, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(255, 255, 255, 255));\n"
"background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(69, 83, 100, 255), stop:1 rgba(255, 255, 255, 255));\n"
"background-color: qlineargradient(spread:pad, x1:1, y1:0, x2:1, y2:0, stop:0 rgba(69, 83, 100, 255), stop:1 rgba(255, 255, 255, 255));\n"
"color: rgb(255, 255, 255);\n"
"background-image: url(:/background/버튼.png);\n"
"font: 10pt \"맑은 고딕\";")
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.gridLayout.addWidget(self.comboBox_2, 0, 1, 1, 1)
        self.pushButton = QtWidgets.QPushButton(Form)
        self.pushButton.setStyleSheet("QPushButton {\n"
"color: qlineargradient(spread:pad, x1:0, y1:1, x2:1, y2:0, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(255, 255, 255, 255));\n"
"background-image: url(:/background/버튼.png);\n"
"font: 10pt \"맑은 고딕\";\n"
"}\n"
"QPushButton:hover{\n"
"border-image: url(:/background/버튼클릭.png);\n"
"}")
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 0, 2, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(Form)
        self.pushButton_2.setStyleSheet("QPushButton {\n"
"color: qlineargradient(spread:pad, x1:0, y1:1, x2:1, y2:0, stop:0 rgba(255, 255, 255, 255), stop:1 rgba(255, 255, 255, 255));\n"
"background-image: url(:/background/버튼.png);\n"
"font: 10pt \"맑은 고딕\";\n"
"}\n"
"QPushButton:hover{\n"
"border-image: url(:/background/버튼클릭.png);\n"
"}")
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout.addWidget(self.pushButton_2, 0, 3, 1, 1)
        self.label = QtWidgets.QLabel(Form)
        self.label.setText("")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 4)
        self.gridLayout_2.addLayout(self.gridLayout, 0, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)
        self.pushButton.clicked.connect(self.label_show)
        self.pushButton_2.clicked.connect(self.label.clear)
        
        self.label.setScaledContents(True)
        
    def label_show(self):
        fluid = self.comboBox.currentText()
        type = self.comboBox_2.currentText()
        if type == "P-h diagram":
            
            if fluid == "R22" :
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/phR22.png")
                
            elif fluid == "R12":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/phR12.png")
            elif fluid == "R134a":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/phR134a.png")
            elif fluid == "R410a":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/phR410a.png")
            elif fluid == "water":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/phwater.png")
            time.sleep(1)               
            self.label.setPixmap(pixmap)
            
        elif type == "T-s diagram":
           
            if fluid == "R22" :
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/TsR22.png")
                
            elif fluid == "R12":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/TsR12.png")
            elif fluid == "R134a":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/TsR134a.png")
            elif fluid == "R410a":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/TsR410a.png")
            elif fluid == "water":
                label = QLabel(self)
                pixmap = QPixmap("C:/파이썬/diagramviewer/Tswater.png")
            time.sleep(1)   
            self.label.setPixmap(pixmap)            

        elif type == "Psychrometric Chart":
            
            label = QLabel(self)
            pixmap = QPixmap("C:/파이썬/diagramviewer/Psychro.png")
            self.label.setPixmap(pixmap)
    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.comboBox.setItemText(0, _translate("Form", "water"))
        self.comboBox.setItemText(1, _translate("Form", "R12"))
        self.comboBox.setItemText(2, _translate("Form", "R22"))
        self.comboBox.setItemText(3, _translate("Form", "R134a"))
        self.comboBox.setItemText(4, _translate("Form", "R410a"))
        self.comboBox_2.setItemText(0, _translate("Form", "P-h diagram"))
        self.comboBox_2.setItemText(1, _translate("Form", "T-s diagram"))
        self.comboBox_2.setItemText(2, _translate("Form", "Psychrometric Chart"))
        self.pushButton.setText(_translate("Form", "Draw"))
        self.pushButton_2.setText(_translate("Form", "Clear"))



if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
