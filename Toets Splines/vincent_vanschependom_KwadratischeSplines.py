import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

x = sp.symbols('x')


def rechte(begin, einde):
    """
    De functie genereert het voorschrift van de rechte tussen twee punten begin(x0,y0) en einde(x1,y1).

    Parameters
    ----------
    begin : list[int1,int2]
        Het beginpunt is een lijst met een x-coördinaat en een y-coördinaat.
    einde : list[int1,int2]
        Het eindpunt is een lijst met een x-coördinaat en een y-coördinaat.

    Returns
    -------
    function
        De functie die de rechte door het begin- en eindpunt beschrijft.
    """
    x1 = begin[0]
    y1 = begin[1]

    x2 = einde[0]
    y2 = einde[1]

    rico = sp.Rational((y2-y1)/(x2-x1))

    def g(x): return rico*(x-x1) + y1

    return g


def parabool(begin, einde, rico):
    """
    De functie genereert het voorschrift van de parabool tussen twee punten begin(x0,y0) en einde(x1,y1), met als richtingscoëfficiënt in het beginpunt de parameter rico.

    Parameters
    ----------
    begin : list[int1,int2]
        Het beginpunt is een lijst met een x-coördinaat en een y-coördinaat.
    einde : list[int1,int2]
        Het eindpunt is een lijst met een x-coördinaat en een y-coördinaat.
    rico : int/sympy.core.symbol.Symbol
        De rico voor de parabool in het beginpunt.

    Returns
    -------
    function
        De functie die de parabool door het begin- en eindpunt beschrijft met passende rico.
    """
    a, b, c = sp.symbols('a, b, c')

    def f(x): return a*(x**2)+b*x+c
    def df(x0): return sp.diff(f(x), x).subs(x, x0)

    oplossing = sp.solve([sp.Eq(df(begin[0]), rico),
                          sp.Eq(f(begin[0]), begin[1]),
                          sp.Eq(f(einde[0]), einde[1])], a, b, c)

    a = oplossing[a]
    b = oplossing[b]
    c = oplossing[c]

    return f


def koorden(X, Y):
    """
    Deze functie berekent een lijst van functies van rechten door twee opeenvolgende punten, gegeven een
    lijst X van x-waarden en een lijst Y van y-waarden.

    Parameters
    ----------
    X : list[int]
        Een lijst met x-waarden.
    Y : list[int]
        Een lijst met y-waarden.

    Returns
    -------
    list[function]
        Een lijst van functies van rechten door twee opeenvolgende punten, gegeven een
        lijst X van x-waarden en een lijst Y van y-waarden.
    """
    gi = list()

    for i in range(len(X)-1):
        gi.append(rechte([X[i], Y[i]], [X[i+1], Y[i+1]]))

    return gi


def splines(X, Y, rico):
    """
    Deze functie berekent een lijst van functies van parabolen door twee opeenvolgende punten, gegeven
    een lijst X van x-waarden, een lijst Y van y-waarden en de afgeleide rico in het eerste punt.

    Parameters
    ----------
    X : list[int]
        Een lijst met x-waarden.
    Y : list[int]
        Een lijst met y-waarden.
    rico : int/sympy.core.symbol.Symbol
        De rico in het eerste punt. Kan ook een symbool zijn.

    Returns
    -------
    list
        Een lijst van functies van parabolen door twee opeenvolgende punten, gegeven
        een lijst X van x-waarden, een lijst Y van y-waarden en de afgeleide rico in het eerste punt.
    """
    fi = list()
    rico_i = rico

    for i in range(len(X)-1):
        fi.append(parabool([X[i], Y[i]], [X[i+1], Y[i+1]], rico_i))
        rico_i = sp.diff(fi[i](x), x).subs(x, X[i+1])

    return fi


def optimaal(X, Y):
    """
    Deze functie berekent de kwadratische oppervlakte S(m) tussen de koorden en de splines in functie van de parameter m. 
    De procedure print de vergelijking S(m) af, tekent de grafiek (m, S(m)) en heeft als uitvoer de waarde waarvoor deze oppervlakte extremaal is.

    Parameters
    ----------
    X : list[int]
        Een lijst met x-waarden.
    Y : list[int]
        Een lijst met y-waarden.

    Returns
    -------
    sympy.core.numbers.Rational
        De optimale waarde voor m.
    """
    m = sp.symbols('m')
    gi = koorden(X, Y)
    fi = splines(X, Y, m)

    som = 0

    for i in range(len(X)-1):
        som += sp.integrate((fi[i](x)-gi[i](x))**2, (x, X[i], X[i+1]))

    def s(x0): return som.subs(m, x0)
    def ds(x0): return sp.diff(s(x), x).subs(x, x0)

    extremum = sp.solve(ds(m), m)[0]

    print(s(m))

    S = sp.lambdify(x, s(x), 'numpy')

    xx = np.linspace(-100, 100, 1000)
    yy = S(xx)
    plt.figure()
    plt.plot(xx, yy)
    plt.show()
    plt.savefig('oppervlakteVanSchependomVincent.pdf')

    return extremum


def figuur(X, Y):
    """
    Deze functie maakt een grafiek van de punten, koorden, splines en raaklijnen in het geval de kwadratische
    oppervlakte minimaal is.

    Parameters
    ----------
    X : _type_
        _description_
    Y : _type_
        _description_
    """
    m = optimaal(X, Y)

    gi = koorden(X, Y)
    fi = splines(X, Y, m)

    plt.figure()

    for i in range(len(gi)):

        xx1 = np.linspace(X[i], X[i+1], 2)
        yy1 = gi[i](xx1)
        plt.plot(xx1, yy1, 'b')

        xx2 = np.linspace(X[i], X[i+1], 100)
        yy2 = fi[i](xx2)
        plt.plot(xx2, yy2, 'r')

        afg = sp.diff(fi[i](x), x).subs(x, X[i])

        xx3 = [X[i]-.5, X[i]+.5]
        yy3 = [fi[i](X[i])-.5*afg, fi[i](X[i])+.5*afg]

        plt.plot(xx3, yy3, 'g--')

    afg = sp.diff(fi[-1](x), x).subs(x, X[-1])

    xx4 = [X[-1]-.5, X[-1]+.5]
    yy4 = [fi[-1](X[-1])-.5*afg, fi[-1](X[-1])+.5*afg]

    plt.plot(xx4, yy4, 'g--')

    plt.scatter(X, Y, c='k')
    plt.savefig('splinesVanSchependomVincent.pdf')
