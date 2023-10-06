import sympy as sp

print("")
print("")


# functievoorschrift
def f(x): return (x**3)/(x**2-1)
def sign(x): return sp.sign(x)

# symbolen definiëren
x, y = sp.symbols('x, y')

kritischePunten = []

# nulwaarden van de noemer
nulwNoemer = sp.solve(sp.denom(f(x)), x)
# nulwaarden van de teller
nulwTeller = sp.solve(sp.numer(f(x)), x)
# nulwaarden van f
nulwFunctie = sp.solve(f(x))

print("Nulwaarden noemer: ", nulwNoemer)
print("Nulwaarden teller: ", nulwTeller)

print("Nulpunten van f:")
for nulwaarde in nulwFunctie:
    kritischePunten.append(nulwaarde)
    print('( 0 ,', f(nulwaarde),')')

print('')

for nulwaarde in nulwNoemer:
    kritischePunten.append(nulwaarde)

kritischePunten.sort()

tekenschema = []

i = 0;
while i < len(kritischePunten) :

    if i == 0:

        if sign(f(kritischePunten[i]-0.001)) == 1:
            tekenschema.append("+")
        if sign(f(kritischePunten[i]-0.001)) == -1:
            tekenschema.append("-")

    if kritischePunten[i] in nulwNoemer:
        tekenschema.append("|")
    else:
        tekenschema.append("0")

    if sign(f(kritischePunten[i]+0.001)) == 1:
        tekenschema.append("+")
    if sign(f(kritischePunten[i]+0.001)) == -1:
        tekenschema.append("-")

    i+=1


print('Tekenschema:')
print(tekenschema)


print('')

verticaleAsymptoten = []
horizontaleAsymptoten = []

# overloop alle nulwaarden van de noemer, dus potentiële VA
for nulwaarde in nulwNoemer :

    linkerLimiet = sp.limit(f(x), x, nulwaarde, dir="-")
    rechterLimiet = sp.limit(f(x), x, nulwaarde, dir="+")

    if linkerLimiet == -sp.oo or linkerLimiet == sp.oo or rechterLimiet ==-sp.oo or rechterLimiet == sp.oo:
        
        verticaleAsymptoten.append(nulwaarde)
        print('De rechte met vergelijking x =', nulwaarde, "is een verticale asymptoot")



# horizontale asymptoot
positieveLimietHorAsymp = sp.limit(f(x), x, sp.oo)
negatieveLimietHorAsymp = sp.limit(f(x), x, -sp.oo)

# als de limiet niet naar oneindig gaat is het een HA
if positieveLimietHorAsymp != sp.oo and positieveLimietHorAsymp != -sp.oo:

    horizontaleAsymptoten.append(positieveLimietHorAsymp)
    print('f heeft een horizontale asymptoot voor x gaande naar', sp.oo)

# als de limiet niet naar oneindig gaat is het een HA
if negatieveLimietHorAsymp != sp.oo and negatieveLimietHorAsymp != -sp.oo:

    horizontaleAsymptoten.append(negatieveLimietHorAsymp)
    print('f heeft een horizontale asymptoot voor x gaande naar', -sp.oo)


# loop over de VA
for a in verticaleAsymptoten:

    # voorwaarde 1
    if sp.limit((f(x)/x), x, sp.oo) == a:

        # als er HA zijn
        if len(horizontaleAsymptoten) != 0:
            
            for b in horizontaleAsymptoten:

                # voorwaarde 2
                if sp.limit(f(x)-a*x) == b:

                    rechte = a*x+b
                    print(rechte, "is een SA")

        # als er geen HA zijn
        else:

            rechte = a*x
            print("De rechte met vergelijking y =", rechte, "is een SA")


# afgeleide functie van f
def afgeleide(x0): return sp.diff(f(x),x).subs(x,x0)

