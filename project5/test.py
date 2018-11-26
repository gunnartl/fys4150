import numpy as np

"""
def s_d():
    


def runge_kutta4(S, I, f):

    k1 = f(t_n, y_n)
    k2 = f(t_n + h/2, y_n + (h/2)*k1)
    k3 = f(t_n + h/2, y_n + (h/2)*k2)
    k4 = f(t_n + h,   y_n + h * k3)
    

    return y_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

"""

def runge_kutta4(y_n, t_n, h, f):
    '''
    Gir y_n+1 ved hjelp av Runge Kutta 4.
    '''
    # Konstanter
    k1 = f(t_n, y_n)
    k2 = f(t_n + h/2, y_n + (h/2)*k1)
    k3 = f(t_n + h/2, y_n + (h/2)*k2)
    k4 = f(t_n + h,   y_n + h * k3)
    
    # Regn ut og returner neste verdi
    return y_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

def f(t, y):
    '''
    Funksjon til den deriverte.
    '''
    return 0.5 * y
    
# Startverdi, skrittlengde, starttid og sluttid
y_0, h, t_0, t_end = 1, 2, 0, 10 
tider = [i for i in range(t_0, t_end+h, h)] # Liste for tiden
y_verdier = [y_0]                           # Liste for y-verdier

for i in range(0, len(tider)-1):
    y_n, t_n = y_verdier[i], tider[i]               # Hent y og t
    y_verdier.append(runge_kutta4(y_n, t_n, h, f))  # Regn ut neste verdi