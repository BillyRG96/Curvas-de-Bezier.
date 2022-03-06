# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 07:14:36 2018

@author: familia romero
"""
### APROXIMACION A LA VELOCIDAD EN CURVA POR OTRO MÉTODO DE CURVATURA
import numpy as np
import math
from random import uniform
from scipy.misc import *
from scipy.special import comb
from math import *
import csv
import matplotlib.pyplot as plt
from sympy import *
from sympy import init_session
import pandas as pd
### aQUI VA LOS DOCUMENTOS EN EXCEL
###LISTA DE EXCEL DE CADA CURVA
E = int(input('Elegir curva:\n Elije un número entre 1 y 39: \n'))
curva = ['NoNe', '1.xlsx', '2.xlsx',
         '3.xlsx', '4.xlsx', '5.xlsx',
         '6.xlsx', '7.xlsx', '8.xlsx',
         '9.xlsx', '9.1.xlsx', '10.xlsx',
         '11.xlsx', '11.1.xlsx', '12.xlsx',
         '13.xlsx', '13.1.xlsx', '14.xlsx',
         '15.xlsx', '15.1.xlsx', '16.xlsx',
         '17.xlsx', '18.xlsx', '19.xlsx',
         '20.xlsx', '20.1.xlsx', '21.xlsx',
         '22.xlsx', '23.xlsx', '23.1.xlsx',
         '24.xlsx', '25.xlsx', '26.xlsx', '27.xlsx',
         '27.1.xlsx', '28.xlsx', '29.xlsx', '30.xlsx',
         '31.xlsx', '32.xlsx']

excel = pd.ExcelFile(curva[E])
print (excel.sheet_names)
colum = excel.parse('Hoja1')
print (colum)
CX, CY = [], []
for i in colum['X']:
    i = i*100000 #Este valor depende de la escala *** averiguar.
    CX.append(i)  
print ('Cordenada x: \n ', CX)
for i in colum['Y']:
    i = i*100000 #Este valor depende de la escala *** averiguar.
    CY.append(i)
print ('coordenada y: \n ', CY)
Q = list(zip(CX, CY))
print(len(Q))
medi = len(CX) - 1
print (Q)
print('AQUI VA LA DISTANCIA ENTRE DOS PUNTOS.: ')
print(CX[0], CY[0])
print(CX[1], CY[1])
print(np.sqrt((CX[0]- CX[1])**2 + (CY[0]-CY[1])**2))
## CX y CY son las listas de puntos sobre la curva.
## AQUI EMPIEZA LA MANIPULACION DE LOS DATOS DE EXCEL CON MATEMATICAS
## Se halla la matriz de bezier con respecto a los tiempos. B(t_0)= m x n
t = [x for x in np.arange(0, (medi + 1)/ medi, 1/medi)]
T = t
m = len (Q) - 1
l = []
def Po(t): 
    for c in t:
        B = 0
        b = []
        for v in range(0, m+1):
            B = comb(m, v) * ( (c)**(v) ) * (1 - c)**(m-v)
            b.append(B)
        l.append(b)
Po(T)
print (l)
lt = np.transpose(l)
print(lt.shape)
ltt = np.transpose(lt)
print(ltt.shape)
print ("Matriz de puntos sobre la curva: \n " + str(Q))
print ("Matriz b(tiempo) \n " + str(ltt))
## Se multiplican las matrices para el metodo de interpolacion
## y así hallar los puntos de control. (Desfinity) 
mul = lt.dot(ltt)##(30, 30)
print(mul.shape)
inv = np.linalg.inv(mul)##(30, 30)
print(inv.shape)
Mmul = inv.dot(lt)
Q1 = np.asarray(Q)
print(Mmul.shape)
Desfinity = Mmul.dot(Q1)###SE CAMBIARON EL ORDEN DE LAS MATRICES.
print ("Puntos de control de la curva a determinado tiempo.\n" + str(Desfinity))
### Aqui añadi el otro codigo...................
## Este código separa las cordenadas de los puntos de control
## para mejorar su manejabilidad.
ll = list(Desfinity.tolist())
x = []
y = []
for h in Desfinity:
    for g in h:
        break
    x.append(h[0])
print (x) 
for h in Desfinity:
    for g in h:
        break
    y.append(h[1])
print (y)
n = len(ll) - 1
print  ('Los puntos de control son: ', str(ll))
## Aquí se generan los diferentes polinomios de bezier per sin ser
## reemplazados en los diferentes tiempos.
t = symbols('t')
S = np.arange(0, n + 1)
l = []
for i in S:
    def b(t): 
        B = comb(n, i) * ( (t)**(n-i) ) * (1 - t)**i
        return B
    l.append(b(t))
##Empieza a mostrar los polinomios de bezier parametricamente.
print ('El arreglo del polinomio es: ' + str(l))
h = np.multiply(x, l)
hh = np.multiply(y, l)
print ('Matriz (X por b(t)) ' + str(h))
print ('Matriz (Y por b(t)) ' + str(hh))
s = str(np.asarray(np.sum(h)))
print ('Este es el polinomio de b(t)*x_i')
print (s)
ss = str(np.asarray(np.sum(hh)))
print ('Este es el polinomio de b(t)*y_i')
print (ss)
B_t = np.array([s, ss])
print ('curva parametrica: \n' + str(B_t))
## Acá empieza el primer metodo para hallar la curvatura.
## Este es de ecuacione  diferenciales con sus respectivas derivadas.
def X(t):
    return s
dXdt = diff(s, t)
print (dXdt) 
def Y(t):
    return ss
dYdt = diff(ss, t)
print (dYdt)
def Yi(t):
    return ss
def Xi(t):
    return s
dXdt2 = diff(str(s), t, 2)
dYdt2 = diff(str(ss), t, 2)
print (dXdt2)
print (dYdt2)
def DEF(t):
    K = (dXdt * dYdt2 - dXdt2 * dYdt)/((dXdt)**2 + (dYdt)**2)**(3/2)
    return K
K = DEF(t)
k = K
print ('Apartir de acá empieza la ecuación de curvatura')
print (k)
## NUEVO METODO DE CURVATURA.
## VECTOR VELOCIDAD
print ("VECTOR DE VELOCIDAD (derivada de la trayectoria)")
vel_1 = [dXdt, dYdt]
print (vel_1)
print ("VECTOR ACELERACION (doble derivada de la trayectoria)")
ace = [dXdt2, dYdt2]
print (ace)
################################################

## Con la curvatura se puede hallar el Radio de curvatura (R). Este es 
## el método para un "t" en específico.

#tt = eval(input('ingrese el valor de t en el cual hallar curvatura \n t = '))
anot = []
inot = []
for t in np.linspace(0, 1, medi):
    j = str(k)
    #print ("Todos los tiempos- curvatura: " + str(eval(j)))
    vel_0 = str(vel_1)
    #print ("La velocidad vectorial es: " + str(eval(vel_0)))
    ace0 = str(ace)
    #print ("La aceleración vectorial es: " + str(eval(ace0)))
    j = str(k)
    #print ('la curvatura es: \n C = ' + str(eval(j)))
    R = 1/abs(eval(j))
    print ('Este valor representa el radio de curvatura. \n R = ' + str(R))
    peralte = 6 # en porcentaje
    GRA = 9.78
    miu = 0.2   
    alpha = peralte/100
    ## AHORA ES LA VELOCIDAD CON PERALTE
    def VELOCIDAD(R, alpha):
        V = math.sqrt(R*GRA*(miu*np.cos(alpha)+ np.sin(alpha))/(np.cos(alpha) - miu*np.sin(alpha)))
        return V
    velperalte = VELOCIDAD(R, miu)
    velsi = velperalte*3.6
    print ("Este es la velocidad con el peralte en km/h: ")
    print (velsi)
    anot.append(velsi)
    ## Se empieza con hallar la velocidad con los valores que se tienen.
    ## ESTA SIGUIENTE ES LA VELOCIDAD PLANA
    def VEL(miu, R):
        V = sqrt(miu*GRA*R)
        return V
    #print ('Este es el resultado el km/h sin peralte(v = sqrt(miu*g*R)): ')
    #print (VEL(miu, R)*3.6)
    ## ESTA ECUACION ES PARA HALLAR PERALTE DESPEJANDO DE LA ECUACION DE VELOCIDAD ANTERIOR
    def A_PERALTE(peralvel):
        P = math.atan(((peralvel*3.6) - miu*R*GRA)/(peralvel*3.6)**2 + R*GRA) 
        return P
    ang_De_Peralte = A_PERALTE(velperalte)
    #print ("Este es el angulo del peralte despejado de la ecuacion anterior")
    #print (ang_De_Peralte)
    def FRENAR(velsi, miu, GRA):
        F = velsi**2/(2*GRA*miu)
        return F
    Dfrenado = FRENAR(velperalte, miu, GRA)
    #print ("Este es la distancia mínima de frenado en metros.: ")
    #print ("d = " + str(Dfrenado))
    #print ("AHORA VIENE LA VELOCIDAD CON RESPECTO DE LAS CONDICIONES DE EQUILIBRIO.")
    X_c = 1.462/2
    Y_c = 2.50/2
    def VEQUILIBRIO(R, GRA, alpha):
        X_p = X_c *np.cos(alpha)
        Y_p = Y_c *np.sin(alpha)
        V = sqrt((R*GRA*(Y_p*np.sin(alpha) + X_p*np.cos(alpha)))/(X_p*np.cos(alpha)- Y_p*np.sin(alpha)))
        return V
    VelEquilibrio = VEQUILIBRIO(R, GRA, alpha)
    print ("Este es la velocidad para un carro en dicha curva: ")
    print ("V_E = " + str(VelEquilibrio))
    Velso = VelEquilibrio*3.6
    inot.append(Velso)
 
    


#####NUEVA UBICACION DE GRAFICA
for t in np.linspace(0, 1, 1000):
    fg = ss
    gf = s
    gra1, gra2 = eval(s), eval(ss)
    plt.plot(CX, CY, "ro", markersize=4)
    plt.plot(gra1 , gra2, 'co', markersize=3)   
    plt.plot(x, y, 'go')
    plt.grid(True)
    plt.title(curva[E])
    plt.xlim(min(CX)- 4.5, max(CX)+ 4.5)
    plt.ylim(min(CY)- 4.5, max(CY)+ 4.5)
    plt.plot(-7426360.2031493755, 454128.4553833721, "ko")


for i in range(len(CX)-1):
    #plt.annotate(str(inot[i]) + ' km/h', (CX[i], CY[i]))
    plt.annotate(str(anot[i]) + ' km/h', (CX[i], CY[i]))
    print("Coordenadas de los puntos" + str((CX[i], CY[i])))
    

plt.show()





