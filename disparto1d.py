#Disparo una dimension EDO segundo orden
import numpy as np

uprobu = open('uprobu.dat','w')
vprobv = open('vprobv.dat','w')
solpro = open('solpro.dat','w')


#Condiciones de Frontera
print ('Ingrese las condiciones de frontera \n')
print ('    ')
a = float(input('para y(a), indique el valor de a:    '))
print ('   ')
b = float(input('para y(b), indiqte el valor de b:    '))
print ('\n')
alpha = float(input('y(a):    '))
print ('     ')
betha = float(input('y(b):    '))
print ('\n')

#Numero de particiones
print ('Numero de particiones \n')
n = int(input('N:     '))

#Avance en x
h = (b-a)/n

#Condiciones dadas por la teoria
U10 = alpha
U20 = 0
V10 = 0
V20 = 1

x = np.arange(a,b+h,h)

#Runge Kutta para 2o orden
#Ejemplo obtenido en Burden. Pag: 676
#Problema u''
def f_1(x,U10,U20):
 return U20

def f_2(x,U10,U20):
 return (-2/x)*U20 + (2/(x**2))*U10 + (np.sin(np.log(x)))/(x**2)

U10v = [U10]
U20v = [U20]

for i in range(0,len(x)-1,1):

 k1 = h*f_1(x[i],U10,U20)

 l1 = h*f_2(x[i],U10,U20)

 k2 = h*f_1(x[i] + h/2, U10 + 0.5*k1, U20 + 0.5*l1)

 l2 = h*f_2(x[i] + h/2, U10 + 0.5*k1, U20 + 0.5*l1)

 k3 = h*f_1(x[i] + h/2, U10 + 0.5*k2, U20 + 0.5*l2) 

 l3 = h*f_2(x[i] + h/2, U10 + 0.5*k2, U20 + 0.5*l2)

 k4 = h*f_1(x[i] + h, U10 + k3, U20 + l3)
 
 l4 = h*f_2(x[i] + h, U10 + k3, U20 + l3)
 
 U10 = U10 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
 U20 = U20 + (1/6)*(l1 + 2*l2 + 2*l3 + l4)

 U10v.append(U10)
 U20v.append(U20) 

for i in range(len(x)):
    uprobu.write(f" {x[i]:.15f}    {U10v[i]:.34f}    \n")
#Metodo de P.Blanchard pag: 589
#Problema v''

def g_1(x,V10,V20):
  return V20

def g_2(x,V10,V20):
  return (-2/x)*V20 + (2/(x**2))*V10

V10v = [V10]
V20v = [V20]


for i in range(0,len(x)-1,1):

 q1 = h*g_1(x[i],V10,V20)

 w1 = h*g_2(x[i],V10,V20)

 q2 = h*g_1(x[i] + h/2, V10 + 0.5*q1, V20 + 0.5*w1)

 w2 = h*g_2(x[i] + h/2, V10 + 0.5*q1, V20 + 0.5*w1)

 q3 = h*g_1(x[i] + h/2, V10 + 0.5*q2, V20 + 0.5*w2) 

 w3 = h*g_2(x[i] + h/2, V10 + 0.5*q2, V20 + 0.5*w2)

 q4 = h*g_1(x[i] + h, V10 + q3, V20 + w3)
 
 w4 = h*g_2(x[i] + h, V10 + q3, V20 + w3)
 
 V10 = V10 + (1/6)*(q1 + 2*q2 + 2*q3 + q4)
 V20 = V20 + (1/6)*(w1 + 2*w2 + 2*w3 + w4)

 V10v.append(V10)
 V20v.append(V20)

for i in range(len(x)):
    vprobv.write(f" {x[i]:.15f}    {V10v[i]:.34f}    \n")
    

#Coeficiente constante
print ('\n')
mu = (betha-U10v[-1])/V10v[-1]

Yv = []

for j in range(0,len(x),1):
  Y = U10v[j] + mu*V10v[j]
  Yv.append(Y)

for i in range(len(x)):
    solpro.write(f" {x[i]:.15f}     {Yv[i]:.34f}   \n")

