import sympy as sp

def asymptootOneindig(f):
    """
    Bepaalt de horizontale of schuine asymptoot van een functie bij x -> oneindig. Een passende boodschap wordt afgedrukt en de functie die de asymptoot beschrijft wordt teruggegeven.

    Parameters
    ----------
    f : function
        De wiskundige functie waarvan de asymptoot wordt bepaald.

    Returns
    -------
    function or None
        De vergelijking van de horizontale of schuine asymptoot, of None als er geen asymptoot is.

    Example
    -------
    >>> def f(x): return x**2/(x + 1)
    >>> asymptootOneindig(f)
    Er is een schuine asymptoot met vergelijking y= x - 1

    Notes
    -----
    De functie maakt gebruik van de SymPy-bibliotheek voor het bepalen van limieten en het uitvoeren van
    wiskundige bewerkingen. Zorg ervoor dat je SymPy hebt geïnstalleerd en geïmporteerd.

    """
    x = sp.symbols('x', real=True)
    limiet = sp.limit(f(x),x,sp.oo)
    if limiet.is_finite:
        def asymptoot(x): return limiet
        print('Er is een horizontale asymptoot met vergelijking y=', asymptoot(x))
    else:
        a = sp.limit(f(x)/x,x,sp.oo)
        b = sp.limit(f(x)-a*x,x,sp.oo)
        if a.is_finite and b.is_finite:
            def asymptoot(x): return a*x+b
            print('Er is een schuine asymptoot met vergelijking y=', asymptoot(x))
        else:
            asymptoot = None
            print('Er is geen horizontale of schuine asymptoot')
    return asymptoot

def klasseerNullen(f):
    """
    Bepaalt en klasseert de nullen van een gebroken functie in nulpunten, perforaties en asymptoten.

    Parameters
    ----------
    f : function
        De wiskundige functie waarvan de nullen worden geanalyseerd en geklasseerd.

    Returns
    -------
    list of lists
        Een lijst die drie sublijsten bevat. De eerste lijst bevat de nulpunten, de tweede lijst de perforaties,
        en de derde lijst de verticale asymptoten.

    Example
    -------
    >>> def f(x): return x**2/(x + 1)
    >>> klasseerNullen(f)
    0 -- nulpunt
    -1 -- asymptoot
    [[0], [], [-1]]

    Notes
    -----
    De functie maakt gebruik van de SymPy-bibliotheek voor het oplossen van vergelijkingen en het bepalen
    van limieten. Zorg ervoor dat je SymPy hebt geïnstalleerd en geïmporteerd.

    """
    x = sp.symbols('x', real=True)
    nulpuntenTeller = sp.solve(sp.numer(f(x)),x)
    nulpuntenNoemer = sp.solve(sp.denom(f(x)),x)
    nullen = set(nulpuntenTeller + nulpuntenNoemer)
    nulpunten = []
    perforaties = []
    asymptoten = []
    for k in nullen:
        if k not in nulpuntenNoemer:
            print(k,'-- nulpunt')
            nulpunten += [k]
        elif k in nulpuntenTeller and sp.limit(f(x),x,k).is_finite:
            print(k,'-- perforatie')
            perforaties += [k]
        else:
            print(k,'-- asymptoot')
            asymptoten += [k]
    return [nulpunten,perforaties,asymptoten]

def raaklijn(f, a):
    """
    Bepaalt de vergelijking van de raaklijn aan een gegeven functie `f` op het punt `x=a`.

    Parameters
    ----------
    f : function
        De wiskundige functie waaraan de raaklijn wordt bepaald.
    a : numeric
        Het punt waarop de raaklijn aan de functie `f` wordt bepaald.

    Returns
    -------
    function
        De vergelijking van de raaklijn aan de functie `f` op het punt `x=a`.

    Example
    -------
    >>> def f(x): return x**2
    >>> x = sp.symbols('x', real=True)
    >>> raaklijn(f,1)(x)
    2x−1

    Notes
    -----
    De functie maakt gebruik van de SymPy-bibliotheek voor het differentiëren en evalueren van
    wiskundige uitdrukkingen. Zorg ervoor dat je SymPy hebt geïnstalleerd en geïmporteerd.

    """
    x = sp.symbols('x', real=True)
    def Df(x0): 
        return sp.diff(f(x), x).subs(x, x0)
    def rkl(x): 
        return Df(a) * (x - a) + f(a)
    return rkl

def klasseerKritiekePunten(f):
    """
    Bepaalt en klasseert de kritieke punten van een gegeven functie `f` in minima, maxima en buigpunten.

    Parameters
    ----------
    f : function
        De wiskundige functie waarvan de kritieke punten worden geanalyseerd en geklasseerd.

    Returns
    -------
    list of lists
        Een lijst die drie sublijsten bevat. De eerste lijst bevat de minima, de tweede lijst de maxima,
        en de derde lijst de buigpunten.

    Example
    -------
    >>> def f(x): return x**3-3*x**2+2*x
    >>> klasseerKritiekePunten(f)
    1 - sqrt(3)/3 -- maximum
    sqrt(3)/3 + 1 -- minimum
    1 -- buigpunt
    De buigraaklijn heeft vergelijking y= 1 - x
    [[sqrt(3)/3 + 1], [1 - sqrt(3)/3], [1]]

    Notes
    -----
    De functie maakt gebruik van de SymPy-bibliotheek voor het differentiëren en oplossen van vergelijkingen.
    Zorg ervoor dat je SymPy hebt geïnstalleerd en geïmporteerd. Daarnaast wordt de eerder gedefinieerde
    `raaklijn` functie gebruikt om de vergelijking van de buigraaklijn te bepalen.

    """
    x = sp.symbols('x', real=True)
    def Df(x0): 
        return sp.diff(f(x), x).subs(x, x0)
    def D2f(x0): 
        return sp.diff(Df(x), x).subs(x, x0)
    def D3f(x0): 
        return sp.diff(D2f(x), x).subs(x, x0)

    kritiekePunten = sp.solve(Df(x), x)
    nulpuntenTweedeAfgeleide = sp.solve(D2f(x), x)

    maxima = []
    minima = []
    buigpunten = []

    for k in kritiekePunten:
        if D2f(k) > 0:
            print(k, '-- minimum')
            minima += [k]
        elif D2f(k) < 0:
            print(k, '-- maximum')
            maxima += [k]

    for k in nulpuntenTweedeAfgeleide:
        if D3f(k) != 0:
            print(k, '-- buigpunt')
            print('De buigraaklijn heeft vergelijking y=', raaklijn(f, k)(x))
            buigpunten += [k]

    return [minima, maxima, buigpunten]

def functieVerloop(f):
    """
    Analyseert het verloop van een gegeven functie `f` en drukt informatie af over zijn grafiek,
    zoals de vergelijking, horizontale/schuine asymptoten, merkwaardige punten, en kritieke punten.

    Parameters
    ----------
    f : function
        De wiskundige functie die geanalyseerd wordt.

    Returns
    -------
    None
        De functie retourneert niets. In plaats daarvan drukt het relevante informatie af over de functie `f`.

    Notes
    -----
    Deze functie is afhankelijk van andere functies zoals `asymptootOneindig`, `klasseerNullen`, en 
    `klasseerKritiekePunten` voor zijn analyse. Zorg ervoor dat deze functies ook in je code zijn gedefinieerd.
    De functie maakt ook gebruik van de SymPy-bibliotheek. Zorg ervoor dat je SymPy hebt geïnstalleerd 
    en geïmporteerd.

    """
    x, a = sp.symbols('x a', real=True)
    
    print('Analyse van de grafiek met vergelijking')
    print('y =', f(x))
    asymptootOneindig(f)
    print('Merkwaardige punten:')
    klasseerNullen(f)
    klasseerKritiekePunten(f)
    
    return None
