import numpy as np
import matplotlib.pyplot as plt

#function needed to create the time domain
def createT(tmin,tmax,N):
    t = np.linspace(tmin,tmax,N)
    return t

#functions needed for the euler method
def odeEuler(f,y0,t):
    y = np.zeros(len(t))
    y[0] = y0
    for n in range(0,len(t)-1):
        y[n+1] = y[n] + f(y[n],t[n])*(t[n+1] - t[n])
    return y

def EulerMethod(tmin,tmax,N,y0,func,title,disp=True):
    t = createT(tmin,tmax,N)
    f = func
    y = odeEuler(f,y0,t)
    if disp == True:
        plt.plot(t,y,'b.-')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.title(title + str(" (N = ") + str(N) + str(" )"))
        plt.grid(True)
        plt.show()
    else: 
        pass
    return t, y

#functions needed for RK2 method
def RKsecondorder(f,x0,t):
    x = np.zeros(len(t))
    x[0] = x0
    for i in range(0,len(t)-1):
        h = t[i+1] - t[i]
        k1 = h * f( x[i], t[i] ) / 2.0
        x[i+1] = x[i] + h * f( x[i] + k1, t[i] + h / 2.0 )
    return x

def RungeKuttaMethodSecondOrder(tmin,tmax,N,y0,func,title):
    t = createT(tmin,tmax,N)
    f = func
    y = RKsecondorder(f,y0,t)
    plt.plot(t,y,'b.-')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title(title + str(" (N = ") + str(N) + str(" )"))
    plt.grid(True)
    plt.show()
    return t,y

#functions needed for RK4 method
def RKfourthorder(f,x0,t):
    x = np.zeros(len(t))
    x[0] = x0
    for i in range(0,len(t)-1):
        dt = t[i+1] - t[i]
        dx1 = dt*f(x[i],t[i])
        dx2 = dt*f(x[i]+0.5*dx1,t[i]+0.5*dt)
        dx3 = dt*f(x[i]+0.5*dx2,t[i]+0.5*dt)
        dx4 = dt*f(x[i]+dx3,t[i]+dt) 
        x[i+1] = x[i] + (1.0/6.0)*(dx1 + 2*dx2 + 2*dx3 + dx4)
    return x
 
def RungeKuttaMethodFourthOrder(tmin,tmax,N,y0,func,title):
    t = createT(tmin,tmax,N)
    f = func
    y = RKfourthorder(f,y0,t)
    plt.plot(t,y,'b.-')
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title(title + str(" (N = ") + str(N) + str(" )"))
    plt.grid(True)
    plt.show()   
    return t,y

#answer to 9b
func1 = lambda y,t: 2*y/t**2
y1 = EulerMethod(1,10,10,1,func1,r'$ x\' = \frac{2x}{t^2} $')

func2 = lambda y,t: 2*y**2/t**2
y2 = EulerMethod(1,1.5,10,1,func2,r'$ x\' = \frac{2x^2}{t^2} $')

#answer to 9c
y1 = EulerMethod(1,10,100,1,func1,r'$ x\' = \frac{2x}{t^2} $')
y2 = EulerMethod(1,1.5,100,1,func2,r'$ x\' = \frac{2x^2}{t^2} $')

#answers to 10
tmin = 0
tmax = 10
x0 = 0
t0 = 0
Narray = [20,50,100,1000]
func3 = lambda y, t: - y**4 + np.sin(2*t)

#answer to 10a
e1 =  EulerMethod(tmin,tmax,Narray[0],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')
e2 =  EulerMethod(tmin,tmax,Narray[1],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')
e3 =  EulerMethod(tmin,tmax,Narray[2],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')

#answer to 10b
rk2i = RungeKuttaMethodSecondOrder(tmin,tmax,Narray[0],x0,func3,r'$ x\' = -x^4 + sin(2*t) $') 
rk2ii = RungeKuttaMethodSecondOrder(tmin,tmax,Narray[1],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')
rk2iii = RungeKuttaMethodSecondOrder(tmin,tmax,Narray[2],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')

#answer to 10c
rk4i = RungeKuttaMethodFourthOrder(tmin,tmax,Narray[0],x0,func3,r'$ x\' = -x^4 + sin(2*t) $') 
rk4ii = RungeKuttaMethodFourthOrder(tmin,tmax,Narray[1],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')
rk4iii = RungeKuttaMethodFourthOrder(tmin,tmax,Narray[2],x0,func3,r'$ x\' = -x^4 + sin(2*t) $')

#answer to 10d
x1,y1 = e1[0],e1[1]
x2,y2 = rk2i[0], rk2i[1]
x3,y3 = rk4i[0], rk4i[1]
x4,y4 = EulerMethod(tmin,tmax,Narray[2],x0,func3,r'$ x\' = -x^4 + sin(2*t) $',False)

plt.plot(x1,y1,"b--",label="Euler N=20")
plt.plot(x2,y2,"r--",label="RK2")
plt.plot(x3,y3,"kx",label="RK4")
plt.plot(x4,y4,"g-",label="Euler N=1000")
plt.title("Comparative Plot")
plt.xlabel("t")
plt.ylabel("x")
plt.legend()
plt.grid(True)
plt.show()

