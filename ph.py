# -*- coding: utf-8 -*-

import numpy as np
from CoolProp.CoolProp import PropsSI


"""
def print_prop(fluid):
    T,P,h = np.asfarray(T), np.asfarray(P), np.asfarray(h)
    n = len(T)
    print("-"*35)
    fmt = "{:>5s}{:>10s}{:>10s}{:>10s}" 
    print(fmt.format('i', 'T', 'P', 'h'))
    print(fmt.format('', 'C', 'kPa', 'kJ/kg'))
    fmt = "{:>5d}{:10.2f}{:10.2f}{:10.2f}"
    for i in range(n):
        print(fmt.format(i,T[i],P[i],h[i]))
    print("-"*35)
"""    
    

print("-"*35)
fmt = "{:>5s}{:>10s}{:>10s}{:>10s}" 
print(fmt.format('i', 'T', 'P', 'h'))
print(fmt.format('', 'C', 'kPa', 'kJ/kg'))
fmt = "{:>5d}{:10.2f}{:10.2f}{:10.2f}"
print("-"*35)




import CoolProp
from CoolProp.Plots import PropertyPlot
plot = PropertyPlot('R22', 'ph')
plot.calc_isolines()
plot.show()


def diagram(fluid, type):
    if type == 'P-h diagram':
        plot = PropertyPlot(fluid, 'ph')
        plot.calc_isolines()
        plot.show()
    elif type == 'T-s diagram':
        ts_plot = PropertyPlot(fluid, 'Ts', tp_limits='ORC')
        ts_plot.calc_isolines(CoolProp.iQ, num=6)
        ts_plot.show()
    
    

   
