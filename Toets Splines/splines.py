# -*- coding: utf-8 -*-
import sympy as sp

x = sp.symbols('x')

def rechte(begin, einde):
    
    x1 = begin[0]
    x2 = begin[1]
    y1 = einde[0]
    y2 = einde[1]
    
    rico = int((y2-y1)/(x2-x1))
    
    def g(x): return rico*(x-x1) + y1
    
    return g

g = rechte([1,3],[3,7])

print(g(x))

def parabool(begin, einde, rico):
    
    a, b, c = sp.symbols('a, b, c')
    
    def f(x): return a*x**2+b*x+c
    def df(x0): return sp.diff(f(x),x).subs(x,x0)
    
    oplossing = sp.solve([sp.Eq(df(begin[0]),rico), 
              sp.Eq(f(begin[0]),begin[1]), 
              sp.Eq(f(einde[0]),einde[1])], a,b,c)
    
    a = oplossing[a]
    b = oplossing[b]
    c = oplossing[c]
    
    return f

f = parabool([1,3],[3,7],1)

print(f(x))
    
    
    
    
    
    
    
    
    
    