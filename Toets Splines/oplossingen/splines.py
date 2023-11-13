##
# Schrijf hier je naam
# Het script splines.py berekent de beste spline door een rij punten.

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


def rechte(begin, eind):
    def koorde(x): return sp.Rational(
        (eind[1]-begin[1])/(eind[0]-begin[0]))*(x-begin[0])+begin[1]
    return koorde


def parabool(begin, eind, rico):
    a, b, c = sp.symbols('a b c')
    def parabool(x): return a*x**2+b*x+c
    def afgeleide(x): return 2*a*x+b
    parameters = sp.solve([parabool(
        begin[0])-begin[1], parabool(eind[0])-eind[1], afgeleide(begin[0])-rico], a, b, c)

    def parabool(x): return parameters[a]*x**2+parameters[b]*x+parameters[c]
    return parabool


def koorden(X, Y):
    functies = []
    for k in range(len(X)-1):
        koorde = rechte([X[k], Y[k]], [X[k+1], Y[k+1]])
        functies += [koorde]
    return functies


def splines(X, Y, rico):
    x = sp.symbols('x')
    m = rico
    functies = []
    for k in range(len(X)-1):
        functies += [parabool([X[k], Y[k]], [X[k+1], Y[k+1]], m)]
        m = sp.diff(functies[-1](x), x).subs(x, X[k+1])
    return functies


def optimaal(X, Y):
    x, m = sp.symbols('x m')
    rechten = koorden(X, Y)
    parabolen = splines(X, Y, m)
    S = 0
    for k in range(len(X)-1):
        S += sp.integrate((parabolen[k](x)-rechten[k]
                          (x))**2, (x, X[k], X[k+1]))
    print(S)
    m0 = sp.solve(sp.diff(S, m), m)[0]
    mi = np.linspace(float(m0)-10, float(m0)+10, 500)
    Si = sp.lambdify(m, S, 'numpy')(mi)
    plt.figure()
    plt.plot(mi, Si)
    plt.savefig('oppervlakte.pdf')
    plt.close()
    return m0


def figuur(X, Y):
    x = sp.symbols('x')
    m = optimaal(X, Y)
    parabolen = splines(X, Y, m)
    plt.figure()
    plt.plot(X, Y)
    for k in range(len(X)-1):
        a = X[k]
        b = X[k+1]
        f = parabolen[k]
        def Df(x0): return sp.diff(f(x), x).subs(x, x0)
        xx = np.linspace(a, b, 50)
        yy = sp.lambdify(x, f(x), 'numpy')(xx)
        plt.plot(xx, yy, 'r', lw=2)
        plt.plot([a, a+(b-a)/3], [f(a), f(a)+Df(a)*(b-a)/3], 'g--')
        plt.plot([b, b-(b-a)/3], [f(b), f(b)-Df(b)*(b-a)/3], 'g--')
    plt.plot(X, Y, 'ko', ms=6)
    plt.savefig('splines.pdf')
    plt.close()
