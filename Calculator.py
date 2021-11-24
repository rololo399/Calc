import numpy as np
from CoolProp.CoolProp import PropsSI
import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5 import uic, QtCore, QtGui, QtWidgets
from hxlib import group
import time
import resources_qr
import diagramviewer


UI_main = "C:/파이썬/ui/Calculator.ui"
class MainWindow(QMainWindow, QWidget):
    def __init__(self):
        super().__init__()
        uic.loadUi(UI_main, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.pushButton.clicked.connect(self.Btn_dimensionless_2)
        self.pushButton_2.clicked.connect(self.Btn_Property)
        self.Btn_Convection.clicked.connect(self.diagram)
        self.Btn_Convection_2.clicked.connect(self.convection_coefficient)
        
    def Btn_Property(self):
        self.close()
        self.second = Properties()
        
    def Btn_dimensionless_2(self):
        self.close()
        self.second = Dimensionless_Calc()
        
    def diagram(self):
        self.close()
        self.second = diagram()
        
    def convection_coefficient(self):
        self.close()
        self.second = Nusselt_Calculator()
        
        
UI_diagram = "C:/파이썬/diagramviewer/diagramviewer.ui"
class diagram(QDialog):
    def __init__(self):
        super(diagram, self).__init__()
        uic.loadUi(UI_diagram, self)
        self.initUI()
        self.show()
        
    def initUI(self):
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

        
UI_Properties = "C:/파이썬/ui/Properties.ui"
class Properties(QMainWindow, QWidget):
    def __init__(self):
        super(Properties, self).__init__()
        uic.loadUi(UI_Properties, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.prop_01)
        
    def back(self):
        self.close()
        self.second = MainWindow()
        
    def prop_01(self):
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T = float(self.lineEdit.text()) + 273
        P = float(self.lineEdit_2.text()) * 1000
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        mu = round(PropsSI('viscosity', 'T', T, 'P', P, fluid), 5)
        rho = round(PropsSI('D', 'T', T, 'P', P, fluid), 3)
        cp = round(PropsSI('C', 'T', T, 'P', P, fluid), 3)
        beta = round(PropsSI('isobaric_expansion_coefficient',
                       'T', T, 'P', P, fluid), 5)
        enthalpy = round(PropsSI('H', 'T', T, 'P', P, fluid), 3)/1000
        self.textBrowser.setText(str(k))
        self.textBrowser_2.setText(str(mu))
        self.textBrowser_3.setText(str(rho))
        self.textBrowser_4.setText(str(cp))
        self.textBrowser_5.setText(str(enthalpy))
        self.textBrowser_6.setText(str(beta))
        

UI_Dimensionless = "C:/파이썬/ui/Dimensionless_Calc.ui"
class Dimensionless_Calc(QDialog, QWidget):
    def __init__(self):
        super(Dimensionless_Calc, self).__init__()
        uic.loadUi(UI_Dimensionless, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.pushButton.clicked.connect(self.Re_calculator)
        self.pushButton_2.clicked.connect(self.Gr_calculator)
        self.pushButton_4.clicked.connect(self.Ra_calculator)
        self.pushButton_3.clicked.connect(self.Pe_calculator)
        self.pushButton_5.clicked.connect(self.back)
        
    def Re_calculator(self):
        self.close()
        self.second = Reynolds_Calculator()
        
    def Gr_calculator(self):
        self.close()
        self.second = Grashof_Calculator()
        
    def Ra_calculator(self):
        self.close()
        self.second = Rayleigh_Calculator()
        
    def Pe_calculator(self):
        self.close()
        self.second = Peclect_Calculator()
        
    def back(self):
        self.close()
        self.second = MainWindow()


UI_Re = "C:/파이썬/ui/ReCalc.ui"
class Reynolds_Calculator(QMainWindow, QWidget):
    
    def __init__(self):
        super(Reynolds_Calculator, self).__init__()
        uic.loadUi(UI_Re, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.toolButton.clicked.connect(self.Reynolds)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Dimensionless_Calc()

    def Reynolds(self):
        T = float(self.lineEdit_2.text()) + 273
        P = float(self.lineEdit_3.text()) * 1000
        V = float(self.lineEdit_4.text())
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        L= float(self.lineEdit_5.text())
        
        mu = PropsSI('viscosity', 'T', T, 'P', P, fluid)
        rho = PropsSI('D', 'T', T, 'P', P, fluid)
        Re = rho*V*L/mu
        if Re < 1000:
            self.textEdit.setText(str(Re))
        elif 1000 <= Re <1000000:
            Re = round(Re/1000, 2)
            self.textEdit.setText(str(Re)+' X 10\u00b3')
        elif 1000000 <= Re <1000000000:
            Re = round(Re/1000000, 2)
            self.textEdit.setText(str(Re)+' X 10\u2076')
        else:
            Re = round(Re/1000000000, 2)
            self.textEdit.setText(str(Re)+' X 10\u2079')
        

UI_Gr = "C:/파이썬/ui/GrCalc.ui"
class Grashof_Calculator(QMainWindow, QWidget):
    
    def __init__(self):
        super(Grashof_Calculator, self).__init__()
        uic.loadUi(UI_Gr, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.toolButton.clicked.connect(self.Grashof)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Dimensionless_Calc()

    def Grashof(self):
        Gr = 0
        g = 9.81
        Ts = float(self.lineEdit_2.text()) + 273
        Ti = float(self.lineEdit_3.text()) + 273
        P = float(self.lineEdit_4.text()) * 1000
        L = float(self.lineEdit_5.text())
        Tm = (Ts + Ti)/2
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        
        beta = PropsSI('isobaric_expansion_coefficient', 'T', Tm, 'P', P, fluid)
        nu =  PropsSI('V', 'T', Tm, 'P', P, fluid)/PropsSI('D','T', Tm, 'P', P, fluid)
        Gr = g * beta * (Ts - Ti) * np.power(L,3) / np.power(nu, 2)
        if Gr < 1000:
            self.textEdit.setText(str(Gr))
        elif 1000 <= Gr <1000000:
            Gr = round(Gr/1000, 2)
            self.textEdit.setText(str(Gr)+' X 10\u00b3')
        elif 1000000 <= Gr <1000000000:
            Gr = round(Gr/1000000, 2)
            self.textEdit.setText(str(Gr)+' X 10\u2076')
        else:
            Gr = round(Gr/1000000000, 2)
            self.textEdit.setText(str(Gr)+' X 10\u2079')

        
UI_Ra = "C:/파이썬/ui/RaCalc.ui"
class Rayleigh_Calculator(QMainWindow, QWidget):
    
    def __init__(self):
        super(Rayleigh_Calculator, self).__init__()
        uic.loadUi(UI_Ra, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.toolButton.clicked.connect(self.Rayleigh)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Dimensionless_Calc()
        
    def Rayleigh(self):
        g = 9.81
        Ts = float(self.lineEdit_2.text()) + 273
        Ti = float(self.lineEdit_3.text()) + 273
        P = float(self.lineEdit_4.text()) * 1000
        L = float(self.lineEdit_5.text())
        Tm = (Ts + Ti)/2
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        
        Gr = group.Grashof(fluid, Ts, Ti, L, P)
        Pr = group.Prandtl(fluid, Tm, P)
        Ra = Gr*Pr
        if Ra < 1000:
            self.textEdit.setText(str(Ra))
        elif 1000 <= Ra <1000000:
            Ra = round(Ra/1000, 2)
            self.textEdit.setText(str(Ra)+' X 10\u00b3')
        elif 1000000 <= Ra <1000000000:
            Ra = round(Ra/1000000, 2)
            self.textEdit.setText(str(Ra)+' X 10\u2076')
        else:
            Ra = round(Ra/1000000000, 2)
            self.textEdit.setText(str(Ra)+' X 10\u2079')
        

UI_Pe = "C:/파이썬/ui/PeCalc.ui"
class Peclect_Calculator(QMainWindow, QWidget):
    
    def __init__(self):
        super(Peclect_Calculator, self).__init__()
        uic.loadUi(UI_Pe, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.toolButton.clicked.connect(self.Peclect)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Dimensionless_Calc()
        
    def Peclect(self):
        T = float(self.lineEdit_2.text()) + 273
        V = float(self.lineEdit_4.text())
        L = float(self.lineEdit_5.text())
        P = float(self.lineEdit_3.text()) * 1000
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Re = group.Reynolds(fluid, T, V, L, P)
        Pr = group.Prandtl(fluid, T, P)
        Pe = Re*Pr
        if Pe < 1000:
            self.textEdit.setText(str(Pe))
        elif 1000 <= Pe <1000000:
            Pe = round(Pe/1000, 2)
            self.textEdit.setText(str(Pe)+' X 10\u00b3')
        elif 1000000 <= Pe <1000000000:
            Pe = round(Pe/1000000, 2)
            self.textEdit.setText(str(Pe)+' X 10\u2076')
        else:
            Pe = round(Pe/1000000000, 2)
            self.textEdit.setText(str(Pe)+' X 10\u2079')
        

# =============================================================================
# Nusselt_Number
# =============================================================================

UI_Nu = "C:/파이썬/ui/NuCal.ui"
class Nusselt_Calculator(QMainWindow, QWidget):
    
    def __init__(self):
        super(Nusselt_Calculator, self).__init__()
        uic.loadUi(UI_Nu, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.forced_convection)
        self.pushButton_4.clicked.connect(self.free_convection)
        
    def back(self):
        self.close()
        self.second = MainWindow()
    
    def forced_convection(self):
        self.close()
        self.second = Nusselt_Calculator_2()
        
    def free_convection(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        

UI_Nu_2 = "C:/파이썬/ui/NuCal_2.ui"
class Nusselt_Calculator_2(QMainWindow, QWidget):
    
    def __init__(self):
        super(Nusselt_Calculator_2, self).__init__()
        uic.loadUi(UI_Nu_2, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.external_flow)
        self.pushButton_3.clicked.connect(self.internal_flow)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator()
    
    def external_flow(self):
        self.close()
        self.second = Nusselt_Calculator_ex()
        
    def internal_flow(self):
        self.close()
        self.second = Nusselt_Calculator_in()
        

        
# =============================================================================
# Nusselt_Number_external flow
# =============================================================================

UI_Nu_ex = "C:/파이썬/ui/NuCal_ex.ui"
class Nusselt_Calculator_ex(QMainWindow, QWidget):
    
    def __init__(self):
        super(Nusselt_Calculator_ex, self).__init__()
        uic.loadUi(UI_Nu_ex, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.external_flow_plate)
        self.pushButton_3.clicked.connect(self.external_flow_cylinder)
        self.pushButton_4.clicked.connect(self.external_flow_sphere)
        self.pushButton_5.clicked.connect(self.Hilpert_nd)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator()
    
    def external_flow_plate(self):
        self.close()
        self.second = external_flow_plate()
        
    def external_flow_cylinder(self):
        self.close()
        self.second = external_flow_cylinder()
        
    def external_flow_sphere(self):
        self.close()
        self.second = external_flow_sphere()
        
    def Hilpert_nd(self):
        self.close()
        self.second = Hilpert_nd()
        

UI_Nu_ex_plate = 'C:/파이썬/ui/external_flow_plate.ui'
class external_flow_plate(QMainWindow, QWidget):
  
    def __init__(self):
        super(external_flow_plate, self).__init__()
        uic.loadUi(UI_Nu_ex_plate, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.external_flow_plate)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_ex()
        
    def external_flow_plate(self):
        
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Ts = float(self.Ts.text()) + 273
        Tinf = float(self.Tinf.text()) + 273
        uinf=float(self.uinf.text())
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        Re_c=5e5
        Tf = (Ts + Tinf)/2
        Re = group.Reynolds(fluid, Tf, uinf, L, P)
        k = round(PropsSI('conductivity', 'T', Tf, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, Tf, P)
        
        if Re < Re_c:
            Nu_x = 0.3387*Re**(1/2)*Pr**(1/3)/((1+(0.0468/Pr)**(2/3))**(1/4))
            # Churchill and Ozoe, valid Pe = Re*Pr > 100
            Nu = 2*Nu_x
        else:
            A = 0.037*Re_c**(4/5) - 0.664*Re_c**(1/2)
            Nu = (0.037*Re**(4/5) - A)*Pr**(1/3)
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))


UI_Nu_ex_cylinder = 'C:/파이썬/ui/external_flow_cylinder.ui'
class external_flow_cylinder(QDialog):
  
    def __init__(self):
        super(external_flow_cylinder, self).__init__()
        uic.loadUi(UI_Nu_ex_cylinder, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.check.clicked.connect(self.external_flow_cylinder)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_ex()
        
    def external_flow_cylinder(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Ts = float(self.Ts.text()) + 273
        Tinf = float(self.Tinf.text()) + 273
        uinf=float(self.uinf.text())
        P = float(self.P.text()) * 1000
        D = float(self.D.text())
        Tf = (Ts + Tinf)/2  # film temperature
        Re = group.Reynolds(fluid, Tf, uinf, D, P)
        k = round(PropsSI('conductivity', 'T', Tf, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, Tf, P)
        d1 = 0.62*Re**(1/2)*Pr**(1/3)
        d2 = (1 + (0.4/Pr)**2/3)**(1/4)
        d3 = (1 + (Re/282000)**(5/8))**(4/5)
        Nu = 0.3 + d1/d2*d3
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/D, 2)
        self.textBrowser_2.setText(str(h))
        

UI_Nu_ex_sphere = 'C:/파이썬/ui/external_flow_sphere.ui'
class external_flow_sphere(QDialog):
  
    def __init__(self):
        super(external_flow_sphere, self).__init__()
        uic.loadUi(UI_Nu_ex_sphere, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.check.clicked.connect(self.external_flow_sphere)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_ex()

    def external_flow_sphere(self):
       # fluid, Ts, Tinf, uinf, D, P=101325
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Ts = float(self.Ts.text()) + 273
        Tinf = float(self.Tinf.text()) + 273
        uinf=float(self.uinf.text())
        P = float(self.P.text()) * 1000
        D = float(self.D.text())
        Tf = (Ts + Tinf)/2
        Re = group.Reynolds(fluid, Tf, uinf, L, P)
        k = round(PropsSI('conductivity', 'T', Tf, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, Tf, P)
        mu = PropsSI('viscosity', 'T', Tinf, 'P', P, fluid)
        mus = PropsSI('viscosity', 'T', Ts, 'P', P, fluid)
        mu_ratio=mu/mus
        Nu = 2 + (0.4*Re**(1/2) + 0.06*Re**(2/3))*Pr**(0.4)*mu_ratio**(1/4)
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))
        
        
UI_Hilpert_nd = 'C:/파이썬/UI/Hilpert_nd.ui'
class Hilpert_nd(QDialog):
  
    def __init__(self):
        super(Hilpert_nd, self).__init__()
        uic.loadUi(UI_Hilpert_nd, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.check.clicked.connect(self.Hilpert_nd)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_ex()

    def Hilpert_nd(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Tinf = float(self.Tinf.text()) + 273
        Ts = float(self.Ts.text()) + 273
        uinf=float(self.uinf.text())
        P = float(self.P.text()) * 1000
        D = float(self.D.text())
        Tf = (Ts + Tinf)/2
        Re = group.Reynolds(fluid, Tf, uinf, D, P)
        k = round(PropsSI('conductivity', 'T', Tf, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, Tf, P)
        
        if 0.4 < Re:
            C = 0.989
            m = 0.330
        elif 4.0 < Re:
            C = 0.911
            m = 0.385
        elif 40 < Re:
            C = 0.683
            m = 0.466
        elif 4000 < Re:
            C = 0.193
            m = 0.618
        elif 40000 < Re < 400000:
            C = 0.027
            m = 0.805
            
        Nu = C*Re**m*Pr**(1/3)
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/D, 2)
        self.textBrowser_2.setText(str(h))
        

# =============================================================================
# Nusselt_Number_internal flow
# =============================================================================

UI_Nu_in = "C:/파이썬/ui/NuCal_in.ui"
class Nusselt_Calculator_in(QMainWindow, QWidget):
    
    def __init__(self):
        super(Nusselt_Calculator_in, self).__init__()
        uic.loadUi(UI_Nu_in, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.thermal_entry)
        self.pushButton_3.clicked.connect(self.combined_entry)
        self.pushButton_4.clicked.connect(self.Dittus_Boelter)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_2()
    
    def thermal_entry(self):
        self.close()
        self.second = thermal_entry()
        
    def combined_entry(self):
        self.close()
        self.second = combined_entry()
        
    def Dittus_Boelter(self):
        self.close()
        self.second = Dittus_Boelter()
        
        
def Graetz(fluid, T, V, D, x, P):
    Re = Reynolds(fluid, T, V, D, P)
    k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
    Pr = group.Prandtl(fluid, T, P)
    Gz = (D/x)*Re*Pr
    return Gz, Re, Pr, k

UI_Nu_in_thermal_entry = 'C:/파이썬/UI/thermal_entry.ui'        
class thermal_entry(QDialog):
  
    def __init__(self):
        super(thermal_entry, self).__init__()
        uic.loadUi(UI_Nu_in_thermal_entry, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.check.clicked.connect(self.thermal_entry)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_in()

    def thermal_entry(self):
        # fluid, T, V, D, x, P=101325
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T = float(self.T.text()) + 273
        V = float(self.V.text())
        D = float(self.D.text())
        P = float(self.P.text()) * 1000
        x = float(self.x.text())
        Gz,Re,Pr,k = Graetz(fluid, T, V, D, x, P)
        
        Nu = 3.66 + 0.0668*Gz/(1 + 0.04*Gz**(2/3))
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/D, 2)
        self.textBrowser_2.setText(str(h))


UI_Nu_in_combined_entry = 'C:/파이썬/UI/combined_entry.ui'
class combined_entry(QDialog):
  
    def __init__(self):
        super(combined_entry, self).__init__()
        uic.loadUi(UI_Nu_in_combined_entry, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.check.clicked.connect(self.combined_entry)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_in()

    def combined_entry(self):
        # fluid, T, V, D, x, P=101325
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T = float(self.T.text()) + 273
        V = float(self.V.text())
        D = float(self.D.text())
        P = float(self.P.text()) * 1000
        x = float(self.x.text())
        Gz, Re, Pr, k = Graetz(fluid, T, V, D, x, P)
        a1 = 3.66/np.tanh(2.264*Gz^(-1/3) + 1.7*Gz^(-2/3))
        a2 = 0.0499*Gz*np.tanh(1/Gz)
        a3 = np.tanh(2.432*Pr^(1/6)*Gz^(-1/6))
        Nu = (a1 + a2)/a3
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/D, 2)
        self.textBrowser_2.setText(str(h))
        
        
UI_Nu_in_Dittus_Boelter = 'C:/파이썬/UI/Dittus_Boelter.ui'
class Dittus_Boelter(QDialog):
  
    def __init__(self):
        super(Dittus_Boelter, self).__init__()
        uic.loadUi(UI_Nu_in_Dittus_Boelter, self)
        self.setting()
        self.show()
        
    def setting(self):
        self.check.clicked.connect(self.Dittus_Boelter)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_in()

    def Dittus_Boelter(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Tm = float(self.Tm.text()) + 273
        Ts = float(self.Ts.text()) + 273
        V =float(self.V.text())
        P = float(self.P.text()) * 1000
        D = float(self.Dh.text())
        Re = group.Reynolds(fluid, Tm, V, D, P)
        k = round(PropsSI('conductivity', 'T', Tm, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, Tm, P)
        n = 0.4 if Ts > Tm else 0.3
        mu = PropsSI('viscosity', 'T', Tm, 'P', P, fluid)      
        mus = PropsSI('viscosity', 'T', Ts, 'P', P, fluid)
        mu_over_mus=mu/mus
        
        if 0.6 <= Pr <= 160:
            Nu = 0.023*Re**(4/5)*Pr**n
        elif Pr <= 16700:
            Nu = 0.027*Re**(4/5)*Pr**(1/3)*mu_over_mus**(0.14)
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/D, 2)
        self.textBrowser_2.setText(str(h))


UI_Nu_in_Siedr = 'C:/파이썬/UI/Sieder_Tate.ui'
class Sieder_Tate(QDialog):
  
    def __init__(self):
        super().__init__()
        uic.loadUi(UI, self)
        self.setting()
        
    def setting(self):
        self.check.clicked.connect(self.Sieder_Tate)

    def Sieder_Tate(self):
        # fluid, Tm, Ts, D, V, P=101325
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Tm = float(self.Tm.text())
        Ts = float(self.Ts.text())
        D = float(self.D.text())
        P = float(self.P.text())
        V = float(self.V.text())
        Re = group.Reynolds(fluid, Tm, V, D, P)
        k = round(PropsSI('conductivity', 'T', Tm, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, Tm, P)
        mu = PropsSI('viscosity', 'T', Tm, 'P', P, fluid)      
        mus = PropsSI('viscosity', 'T', Ts, 'P', P, fluid)
        mu_over_mus=mu/mus
        Nu = 0.027*Re**(4/5)*Pr**(1/3)*mu_over_mus**(0.14)
        
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/D, 2)
        self.textBrowser_2.setText(str(h)+'W/m\u00b2K')
        


# =============================================================================
# Nusselt_Number_free convection
# =============================================================================

UI_Nu_fc = "C:/파이썬/ui/NuCal_fc.ui"
class Nusselt_Calculator_fc(QMainWindow, QWidget):
    
    def __init__(self):
        super(Nusselt_Calculator_fc, self).__init__()
        uic.loadUi(UI_Nu_fc, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.pushButton_6.clicked.connect(self.back)
        self.pushButton.clicked.connect(self.fc_vertical_plate_nd)
        self.pushButton_2.clicked.connect(self.fc_plate_horizontal)
        self.pushButton_3.clicked.connect(self.fc_horizontal_cylinder_nd)
        self.pushButton_4.clicked.connect(self.fc_sphere_nd)
        self.pushButton_5.clicked.connect(self.fc_horizontal_enclosure)
        self.pushButton_7.clicked.connect(self.fc_vertical_cavity)
        self.pushButton_8.clicked.connect(self.fc_concentric_cylinder)
        self.pushButton_9.clicked.connect(self.fc_concentric_sphere)
        self.pushButton_10.clicked.connect(self.fc_tilted_cavity_nd)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator()
    
    def fc_vertical_plate_nd(self):
        self.close()
        self.second = fc_vertical_plate_nd()
        
    def fc_plate_horizontal(self):
        self.close()
        self.second = fc_plate_horizontal()
        
    def fc_horizontal_cylinder_nd(self):
        self.close()
        self.second = fc_horizontal_cylinder_nd()
        
    def fc_sphere_nd(self):
        self.close()
        self.second = fc_sphere_nd()
        
    def fc_horizontal_enclosure(self):
        self.close()
        self.second = fc_horizontal_enclosure()
        
    def fc_vertical_cavity(self):
        self.close()
        self.second = fc_vertical_cavity()
        
    def fc_concentric_cylinder(self):
        self.close()
        self.second = fc_concentric_cylinder()
        
    def fc_concentric_sphere(self):
        self.close()
        self.second = fc_concentric_sphere()
        
    def fc_tilted_cavity_nd(self):
        self.close()
        self.second = fc_tilted_cavity()


UI_fc_vertical_plate_nd = "C:/파이썬/ui/fc_vertical_plate_nd.ui"
class fc_vertical_plate_nd(QDialog):
  
    def __init__(self):
        super().__init__()
        uic.loadUi(UI_fc_vertical_plate_nd, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_vertical_plate_nd)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_vertical_plate_nd(self):
        
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1 = float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, T, P)
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = group.Rayleigh(fluid, T1, T2, L, P)
        
        if Ra < 1e9 :
           Nu = 0.68 + 0.670*Ra**(1/4)/(1+(0.492/Pr)**(9/16))**(4/9)
        else :
           Nu = (0.825 + 0.387*Ra**(1/6)/(1 + (0.492/Pr)**(9/16))**(8/27))**2
           
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))    
        

UI_fc_plate_horizontal = "C:/파이썬/ui/fc_plate_horizontal.ui"
class fc_plate_horizontal(QMainWindow, QWidget):
    
    def __init__(self):
        super(fc_plate_horizontal, self).__init__()
        uic.loadUi(UI_fc_plate_horizontal, self)
        self.initUI()
        self.show()
    
    def initUI(self):
        self.pushButton.clicked.connect(self.back)
        self.pushButton_2.clicked.connect(self.fc_plate_horizontal1_nd)
        self.pushButton_3.clicked.connect(self.fc_plate_horizontal2_nd)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_plate_horizontal1_nd(self):
        self.close()
        self.second = fc_plate_horizontal1_nd()
    
    def fc_plate_horizontal2_nd(self):
        self.close()
        self.second = fc_plate_horizontal2_nd()

UI_fc_plate_horizontal1_nd = "C:/파이썬/ui/fc_plate_horizontal1_nd.ui"
class fc_plate_horizontal1_nd(QDialog):
  
    def __init__(self):
        super(fc_plate_horizontal1_nd, self).__init__()
        uic.loadUi(UI_fc_plate_horizontal1_nd, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_plate_horizontal1_nd)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = fc_plate_horizontal()
        
    def fc_plate_horizontal1_nd(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1=float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        Pr = cv.Prandtl(fluid, T, P)
        Gr = cv.Grashof(fluid, T1, T2, L, P)
        Ra = cv.Rayleigh(fluid, T1, T2, L, P)
        
        if 1e4 < Ra < 1e7 and Pr >= 0.7:
            Nu = 0.54*Ra**(1/4)
        elif 1e7 <= Ra < 1e11:
            Nu = 0.15*Ra**(1/3)
        else:
            Nu = '계산 불가'
            
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))    
        
        
UI_fc_plate_horizontal2_nd = "C:/파이썬/ui/fc_plate_horizontal2_nd.ui"
class fc_plate_horizontal2_nd(QDialog):
  
    def __init__(self):
        super(fc_plate_horizontal2_nd, self).__init__()
        uic.loadUi(UI_fc_plate_horizontal2_nd, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_plate_horizontal2_nd)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = fc_plate_horizontal()
        
    def fc_plate_horizontal2_nd(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1=float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, T, P)
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = group.Rayleigh(fluid, T1, T2, L, P)
        Nu = 0.52*Ra**(1/5)
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))   
        

UI_fc_horizontal_cylinder_nd = "C:/파이썬/ui/fc_horizontal_cylinder_nd.ui"
class fc_horizontal_cylinder_nd(QDialog):
  
    def __init__(self):
        super(fc_horizontal_cylinder_nd, self).__init__()
        uic.loadUi(UI_fc_horizontal_cylinder_nd, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_horizontal_cylinder_nd)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_horizontal_cylinder_nd(self):
        
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1 = float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, T, P)
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = group.Rayleigh(fluid, T1, T2, L, P)
        Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))   
        

UI_fc_sphere_nd = "C:/파이썬/ui/fc_sphere_nd.ui"
class fc_sphere_nd(QDialog):
  
    def __init__(self):
        super(fc_sphere_nd, self).__init__()
        uic.loadUi(UI_fc_sphere_nd, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_sphere_nd)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_sphere_nd(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1 = float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        Pr = group.Prandtl(fluid, T, P)
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = group.Rayleigh(fluid, T1, T2, L, P)
        Nu = 2 + (0.589*Ra**(1/4))/(1+(0.469/Pr)**(9/16))**(4/9)
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
            
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))
        

UI_fc_horizontal_enclosure = "C:/파이썬/ui/fc_horizontal_enclosure.ui"
class fc_horizontal_enclosure(QDialog):
  
    def __init__(self):
        super(fc_horizontal_enclosure, self).__init__()
        uic.loadUi(UI_fc_horizontal_enclosure, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_horizontal_enclosure)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_horizontal_enclosure(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1=float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        Pr = group.Prandtl(fluid, T, P)
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = round(group.Rayleigh(fluid, T1, T2, L, P), 2)
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        Nu = 0.069*Ra**(1/3)*Pr**(0.074)
        h = round(Nu*k/L, 2)
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
        self.textBrowser_2.setText(str(h))
        
        # if Ra < 1000:
        #     self.textEdit.setText(str(Ra))
        # elif 1000 <= Ra <1000000:
        #     Ra = round(Ra/1000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u00b3')
        # elif 1000000 <= Ra <1000000000:
        #     Ra = round(Ra/1000000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u2076')
        # else:
        #     Ra = round(Ra/1000000000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u2079')
        

UI_fc_vertical_cavity = "C:/파이썬/ui/fc_vertical_cavity.ui"
class fc_vertical_cavity(QDialog):
  
    def __init__(self):
        super(fc_vertical_cavity, self).__init__()
        uic.loadUi(UI_fc_vertical_cavity, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_vertical_cavity)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_vertical_cavity(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        T1 = float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        T = (T1 + T2)/2
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        Pr = group.Prandtl(fluid, T, P)
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = round(group.Rayleigh(fluid, T1, T2, L, P), 2)
        k = round(PropsSI('conductivity', 'T', T, 'P', P, fluid), 3)
        H = float(self.H.text())
        
        if 2 <= H/L <= 10 and Pr <= 10**10 and 10**3 <= Ra <= 10**10:
            Nu = 0.22 * ((Pr * Ra / (0.2 + Pr))**0.28) * ((H / L)**(-1/4))
        elif 1 <= H/L <= 2 and 10**(-3) <= Pr <= 10**5 and 10**3 <= Ra * Pr / (0.2 + Pr):
            Nu = 0.18 * ((Pr * Ra / (0.2 + Pr))**0.29)
        elif 10 <= H/L <= 40 and 1 <= Pr <= 2 * 10**4 and 10**4 <= Ra <= 10**7:
            Nu = 0.42*(Ra**(1/4))*(Pr**(0.012))*((H/L)**(-0.3))
        elif 1 <= H/L <= 40 and 1 <= Pr <= 20 and 10**6 <= Ra <= 10**9:
            Nu = 0.046 * (Ra**(1/3))
            
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
        
        h = round(Nu*k/L, 2)
        self.textBrowser_2.setText(str(h))
        
        # if Ra < 1000:
        #     self.textEdit.setText(str(Ra))
        # elif 1000 <= Ra <1000000:
        #     Ra = round(Ra/1000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u00b3')
        # elif 1000000 <= Ra <1000000000:
        #     Ra = round(Ra/1000000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u2076')
        # else:
        #     Ra = round(Ra/1000000000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u2079')
        

UI_fc_concentric_cylinder = "C:/파이썬/ui/fc_concentric_cylinder.ui"
class fc_concentric_cylinder(QDialog):
  
    def __init__(self):
        super(fc_concentric_cylinder, self).__init__()
        uic.loadUi(UI_fc_concentric_cylinder, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.pushButton_2.clicked.connect(self.fc_concentric_cylinder)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_concentric_cylinder(self):
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Ti = float(self.lineEdit.text()) + 273
        To = float(self.lineEdit_2.text()) +273
        ri = float(self.lineEdit_3.text())
        ro = float(self.lineEdit_4.text())
        L = float(self.lineEdit_5.text())
        P = float(self.lineEdit_6.text()) * 1000
        Lc = 2*(np.log(ro/ri))**(4/3)/(ri**(-3/5) + ro**(-3/5))**(5/3)
        k = PropsSI('conductivity', 'T', Ti, 'P', P, fluid)
        Pr = group.Prandtl(fluid, Ti, P)
        Ra = round(group.Rayleigh(fluid, Ti, To, Lc, P), 2)
        k_eff = k*0.386*(Pr/(0.861 + Pr))**(1/4)*Ra**(1/4)
        if k_eff < k:
            k_eff = k
        k_eff = round(k_eff, 2)
        q = round(2*np.pi*L*k_eff*(Ti-To)/np.log(ro/ri), 2)
        self.textBrowser.setText(str(q))
        self.textBrowser_2.setText(str(k_eff))
        
        if Ra < 1000:
            Ra = round(Ra/1000, 2)
            self.textBrowser_3.setText(str(Ra))
        elif 1000 <= Ra <1000000:
            Ra = round(Ra/1000, 2)
            self.textBrowser_3.setText(str(Ra)+' X 10\u00b3')
        elif 1000000 <= Ra <1000000000:
            Ra = round(Ra/1000000, 2)
            self.textBrowser_3.setText(str(Ra)+' X 10\u2076')
        else:
            Ra = round(Ra/1000000000, 2)
            self.textBrowser_3.setText(str(Ra)+' X 10\u2079')
        
        
UI_fc_concentric_sphere = "C:/파이썬/ui/fc_concentric_sphere.ui"
class fc_concentric_sphere(QDialog):

    def __init__(self):
        super(fc_concentric_sphere, self).__init__()
        uic.loadUi(UI_fc_concentric_sphere, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.pushButton_2.clicked.connect(self.fc_concentric_sphere)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_concentric_sphere(self):
        fluid = self.comboBox.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        Ti = float(self.lineEdit.text()) + 273
        To = float(self.lineEdit_2.text()) + 273
        ri = float(self.lineEdit_3.text())
        ro = float(self.lineEdit_4.text())
        P = float(self.lineEdit_6.text()) * 1000
        Lc = (1/ri-1/ro)**(4/3)/(2**(1/3)*(ri**(-7/5)+ro**(-7/5))**(5/3))
        k = PropsSI('conductivity', 'T', Ti, 'P', P, fluid)
        Pr = group.Prandtl(fluid, Ti, P)
        Ra = round(group.Rayleigh(fluid, Ti, To, Lc, P), 2)
        k_eff = k*0.74*(Pr/(0.861 + Pr))**(1/4)*Ra**(1/4)
        if k_eff < k:
            k_eff = k
        k_eff = round(k_eff, 2)
        q = round(4*np.pi*k_eff*(Ti-To)/(1/ri-1/ro), 2)
        self.textBrowser.setText(str(q))
        self.textBrowser_2.setText(str(k_eff))
        
        if Ra < 1000:
            Ra = round(Ra/1000, 2)
            self.textBrowser_3.setText(str(Ra))
        elif 1000 <= Ra <1000000:
            Ra = round(Ra/1000, 2)
            self.textBrowser_3.setText(str(Ra)+' X 10\u00b3')
        elif 1000000 <= Ra <1000000000:
            Ra = round(Ra/1000000, 2)
            self.textBrowser_3.setText(str(Ra)+' X 10\u2076')
        else:
            Ra = round(Ra/1000000000, 2)
            self.textBrowser_3.setText(str(Ra)+' X 10\u2079')
        
        
UI_fc_tilted_cavity = "C:/파이썬/ui/fc_tilted_cavity.ui"
class fc_tilted_cavity(QDialog):
  
    def __init__(self):
        super(fc_tilted_cavity, self).__init__()
        uic.loadUi(UI_fc_tilted_cavity, self)
        self.initUI()
        self.show()
        
    def initUI(self):
        self.check.clicked.connect(self.fc_tilted_cavity)
        self.pushButton.clicked.connect(self.back)
        
    def back(self):
        self.close()
        self.second = Nusselt_Calculator_fc()
        
    def fc_tilted_cavity(self):
        fluid = self.fluid_2.currentText()
        c = fluid.split(' ')
        i = 0
        for fluid in c:
            i += 1
            if i == 1:
                break
        tau = float(self.t.text())
        T1 = float(self.T1.text()) + 273
        T2 = float(self.T2.text()) + 273
        P = float(self.P.text()) * 1000
        L = float(self.L.text())
        Gr = group.Grashof(fluid, T1, T2, L, P)
        Ra = round(group.Rayleigh(fluid, T1, T2, L, P), 2)
        k = PropsSI('conductivity', 'T', T1, 'P', P, fluid)
        a1 = 1 - 1708/(np.cos(tau)*Ra)
        if a1 < 0:
            a1 = 0
        a2 = 1 - 1708*(np.sin(1.8*tau))**(1.6)/(Ra*np.cos(tau))
        a3 = (np.cos(tau)*Ra/5830)**(1/3) - 1
        if a3 < 0:
            a3 = 0
        Nu = 1 + 1.44*a1*a2 + a3
        if Nu < 1000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu))
        elif 1000 <= Nu <1000000:
            Nu = round(Nu/1000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u00b3')
        elif 1000000 <= Nu <1000000000:
            Nu = round(Nu/1000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2076')
        else:
            Nu = round(Nu/1000000000, 2)
            self.textBrowser.setText(str(Nu)+' X 10\u2079')
        h = round(Nu*k/L, 2)

        self.textBrowser_2.setText(str(h))
        
        # if Ra < 1000:
        #     self.textEdit.setText(str(Ra))
        # elif 1000 <= Ra <1000000:
        #     Ra = round(Ra/1000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u00b3')
        # elif 1000000 <= Ra <1000000000:
        #     Ra = round(Ra/1000000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u2076')
        # else:
        #     Ra = round(Ra/1000000000, 2)
        #     self.textEdit.setText(str(Ra)+' X 10\u2079')


        
app = QApplication(sys.argv)
ex = MainWindow()
ex.show()
sys.exit(app.exec_())