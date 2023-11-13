##
# Vincent Van Schependom
# r0976455

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
x, y = sp.symbols('x y')
alpha, beta, gamma, delta = 1.1, 0.4, 0.4, 0.1


def evenwichtspunten(prooi,roofdier):

    stelsel = sp.Matrix([
        prooi(x,y),
        roofdier(x,y)
    ])

    oplossing = sp.solve(stelsel,x,y)

    print(f'Er zijn twee punten waar er een populatie-evenwicht is: {oplossing}')

    return oplossing

def jacobiaan(prooi,roofdier):

    def j(x0,y0): 
        return sp.Matrix([
            [sp.diff(prooi(x,y),x),sp.diff(prooi(x,y),y)],
            [sp.diff(roofdier(x,y),x),sp.diff(roofdier(x,y),y)]
        ]).subs(x,x0).subs(y,y0)

    return j

def real(imag):

    return (imag + sp.conjugate(imag))

def soortevenwicht(J,ptn):

    symbolen = list()

    zadelpunten = list()
    stabielepunten = list()
    spiraalpunten = list()

    for punt in ptn:

        eigenwaarde1 = J(punt[0],punt[1]).eigenvects()[0][0]
        eigenwaarde2 = J(punt[0],punt[1]).eigenvects()[1][0]

        if real(eigenwaarde1) < 0 and real(eigenwaarde2) < 0:

            symbolen.append('+')
            print(f'Het punt {punt} is stabiel')

            stabielepunten.append(punt)

        elif (real(eigenwaarde1) and real(eigenwaarde2) > 0) or (real(eigenwaarde1) > 0 and real(eigenwaarde2) < 0):

            symbolen.append('+')
            print(f'Het punt {punt} is een zadelpunt')

            zadelpunten.append(punt)
        
        elif sp.Eq(sp.conjugate(eigenwaarde1),eigenwaarde2):

            symbolen.append('x')
            print(f'Het punt {punt} is een spiraalpunt')

            spiraalpunten.append(punt)
        
    return symbolen

def punten(J,ptn):

    zadelpunten = list()
    stabielepunten = list()
    spiraalpunten = list()

    for punt in ptn:

        eigenwaarde1 = J(punt[0],punt[1]).eigenvects()[0][0]
        eigenwaarde2 = J(punt[0],punt[1]).eigenvects()[1][0]

        if real(eigenwaarde1) < 0 and real(eigenwaarde2) < 0:

            stabielepunten.append(punt)

        elif (real(eigenwaarde1) and real(eigenwaarde2) > 0) or (real(eigenwaarde1) > 0 and real(eigenwaarde2) < 0):

            zadelpunten.append(punt)
        
        elif sp.Eq(sp.conjugate(eigenwaarde1),eigenwaarde2):

            spiraalpunten.append(punt)
        
    return zadelpunten, stabielepunten, spiraalpunten

def euler(prooi, roofdier, x0, y0, dt, t_max):

    t = 0

    x = [x0]
    y = [y0]

    while t < t_max:

        x.append(x[-1] + dt*prooi(x[-1],y[-1]))
        y.append(y[-1] + dt*roofdier(x[-1],y[-1]))

        t+=dt

    return np.array(x), np.array(y)

def figuur1(prooi, roofdier, x0, y0, dt, t_max):

    xx = np.arange(0,t_max+dt,dt)
    y1 = euler(prooi,roofdier,x0,y0,dt,t_max)[0]
    y2 = euler(prooi,roofdier,x0,y0,dt,t_max)[1]

    plt.figure()
    plt.suptitle('Lotka-Volterra Modle - Euler Methode')
    plt.plot(xx,y1,label='Prooidieren (x)')
    plt.plot(xx,y2,label='Roofdieren (y)')
    plt.xlim(0,t_max)
    plt.legend()
    plt.xlabel('Tijd')
    plt.ylabel('Populatie')
    plt.savefig('figuur1.pdf')
    plt.show()
    plt.close()

    return None

def figuur2(prooi, roofdier, x0, y0, dt, t_max):

    x_euler, y_euler = euler(prooi, roofdier, 10, 10, 0.001, 100)

    J = jacobiaan(prooi,roofdier)
    ptn = evenwichtspunten(prooi, roofdier)
    zadelpunten, stabielepunten, spiraalpunten = punten(J, ptn)

    # zadelpunten = [[(0.00),(0.00)]]
    # spiraalpunten = [[(4.00),(2.75)]]
    # stabielepunten = []

    plt.figure()
    plt.plot(x_euler,y_euler,linewidth=0.5)

    for x,y in zadelpunten:

        plt.scatter(x,y, marker='+')
    
    for x,y in spiraalpunten:

        plt.scatter(x,y, marker='x')

    for x,y in stabielepunten:

        plt.scatter(x,y, marker='o')

    plt.suptitle('Lotka-Volterra Model - Euler Methode')
    plt.xlabel('Prooidierpopulatie')
    plt.ylabel('Roofdierpopulatie')
    plt.savefig('figuur2.pdf')
    plt.show()
    plt.close()

    return None