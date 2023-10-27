import sympy as sp
x = sp.symbols('x', real=True)

def f(x): return (x**3)/(x**2-1)
def g(x): return (x**3-x**2+x-1)/(x**3+x**2-x-1)

def asymptootOneindig(f):

    hor = sp.limit(f(x),x,sp.oo)

    if hor.is_finite:

        return f"Er is een horizontale asymptoot met vergelijking y = {hor}"

    else:

        a = sp.limit(f(x)/x,x,sp.oo)
        b = sp.limit(f(x)-a*x,x,sp.oo)

        if a.is_finite and b.is_finite:

            def SA(x0): return ((a*x)+b).subs(x,x0)
            return f"Er is een schuine asymptoot met vergelijking y = {SA(x)}"
        
        else:

            return "Geen asymptoten"
        

def klasseerNulpunten(f):

    nulpTeller = set(sp.solve(sp.numer(f(x))))
    nulpNoemer = set(sp.solve(sp.denom(f(x))))

    nulpFunctie = nulpTeller.difference(nulpNoemer)
    probleemGevallen = nulpTeller.union(nulpNoemer) - nulpFunctie

    perforaties = []
    asymptoten = []
    nulpunten = list(nulpFunctie)

    for n in probleemGevallen:

        # Het is een perforatie
        if sp.limit(f(x),x,n).is_finite:

            perforaties.append(n)

        else:

            asymptoten.append(n)

    for nulpunt in nulpunten: print(f"{nulpunt} -- nulpunt")
    for perforatie in perforaties: print(f"{perforatie} -- perforatie")
    for asymptoot in asymptoten: print(f"{asymptoot} -- asymptoot")

    # print([nulpunten, perforaties, asymptoten])


def klasseerKritiekePunten(f):

    def df(x0): return sp.diff(f(x),x).subs(x,x0)
    def d2f(x0): return sp.diff(df(x),x).subs(x,x0)
    def d3f(x0): return sp.diff(d2f(x),x).subs(x,x0)

    nulpEerste = sp.solve(df(x))
    nulpTweede = sp.solve(d2f(x))

    minima, maxima, buigpunten = list(), list(), list()

    # Nulpunten van de eerste afgeleide
    for n in nulpEerste:

        # Tweede afgeleide positief
        if d2f(n) > 0:

            print(f"{n} -- minimum")
            minima.append(n)
        
        # Tweede afgeleide negatief
        elif d2f(n) < 0:

            print(f"{n} -- maximum")
            maxima.append(n)

    # Nulpunten van de tweede afgeleide
    for n in nulpTweede:

        if d3f(n) != 0:

            print(f"{n} -- buigpunt")
            buigpunten.append(n)

            def raaklijn(x0): return df(x0)*(x - x0) + f(x0)
            print(f"De buigraaklijn heeft vergelijking y = {raaklijn(n)}")

    # print([minima, maxima, buigpunten])
    

def functieVerloop(f):

    print("\n----------------------------")
    print(f"Analyse van de grafiek met vergelijking y = {f(x)}")
    print(asymptootOneindig(f))
    print("Merkwaardige punten:")
    klasseerNulpunten(f)
    klasseerKritiekePunten(f)


functieVerloop(f)
functieVerloop(g)

# print(asymptootOneindig(f))
# klasseerNulpunten(f)

# print("")

# print(asymptootOneindig(g))
# klasseerNulpunten(g)

# print("")

# klasseerKritiekePunten(f)
# klasseerKritiekePunten(g)