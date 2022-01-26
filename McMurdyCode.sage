def IsoAdd(phi1,phi2,Wrs1,Wrs2):
    f1=phi1[0]
    g1=phi1[1]
    f2=phi2[0]
    g2=phi2[1]
    mm=(g2-g1)/(f2-f1)
    XX=mm^2*Wrs1-Wrs2[2]-f1-f2
    YY=-mm*(XX-f1)-g1


    return [XX,YY]

def EndoAdd(phi1,phi2,Wrs):

    return IsoAdd(phi1,phi2,Wrs,Wrs)

def IsoDouble(phi,Wrs1,Wrs2):
    f=phi[0]
    g=phi[1]
    mm=((Wrs2.derivative())(x=f))/(2*g*Wrs1)
    XX=mm^2*Wrs1-Wrs2[2]-2*f
    YY=-mm*(XX-f)-g
    return [XX,YY]

def IsoAdd(phi1,phi2,Wrs1,Wrs2):
    f1=phi1[0]
    g1=phi1[1]
    f2=phi2[0]
    g2=phi2[1]
    mm=(g2-g1)/(f2-f1)
    XX=mm^2*Wrs1-Wrs2[2]-f1-f2
    YY=-mm*(XX-f1)-g1
    return [XX,YY]

def EndoAdd(phi1,phi2,Wrs):
    return IsoAdd(phi1,phi2,Wrs,Wrs)

def IsoDouble(phi,Wrs1,Wrs2):
    f=phi[0]
    g=phi[1]
    mm=((Wrs2.derivative())(x=f))/(2*g*Wrs1)
    XX=mm^2*Wrs1-Wrs2[2]-2*f
    YY=-mm*(XX-f)-g
    return [XX,YY]

def IsoTriple(phi,Wrs1,Wrs2):
    return IsoAdd(IsoDouble(phi,Wrs1,Wrs2),phi,Wrs1,Wrs2)

def EndoTriple(phi,Wrs):
    return IsoTriple(phi,Wrs,Wrs)

def EndoDouble(phi,Wrs):
    return IsoDouble(phi,Wrs,Wrs)

def IsoOpp(Phi):
    return [Phi[0],-Phi[1]]


#Note: IsoCompose will compose phi2(phi1).
def IsoCompose(phi1,phi2):
    f1=phi1[0]
    g1=phi1[1]
    f2=phi2[0]
    g2=phi2[1]
    return [f2(x=f1),g1*g2(x=f1)]

def IsEndo(Phi,Wrs):
    return (Phi[1])^2*Wrs==Wrs(x=Phi[0])

def DupConjugates(Wrs):
    rr=Wrs.roots(multiplicities=False) #x-coordinates of 2-torsion points
    aa=Wrs[2] #coefficient of x^2 in the Weierstrass eq Wrs
    conjList=[x]
    for i in range(3):
        conjList.append(-aa-rr[i]-x+Wrs/(x-rr[i])^2)
    return conjList #Returns a list of x + (2-torsion pt)

#The function Trip Conjugates assumes Wrs is in *short weierstrass form* (i.e., x^3 + ax + b)
def TripConjugates(Wrs):
    rr = (3*x^4 + 6*Wrs[1]*x^2 + 12*Wrs[0]*x - Wrs[1]^2).roots(multiplicities = False)
    aa = Wrs[2]
    conjList=[x]
    for i in range(4):
        conjList.append(-aa-rr[i]-x+Wrs/(x-rr[i])^2)
    return conjList #Returns a list of x-coordinates of points: x + (3-torsion pt), for half othe nontrivial 3-torsion points.
#Since this only considers x-coordiantes, we don't actually add all of the 3-torsionpoints- Just four of them- not their negatives.

def DupPolyRoots(poly,Wrs):
    Ef2=EndoDouble([x,1],Wrs)[0] #gives the x-coordinate of the point [2](x,1)
    DupConjList=DupConjugates(Wrs)
    normOfPoly=1
    for i in range(4):
        normOfPoly=normOfPoly*poly(x=DupConjList[i])
    #print(Ef2.numerator().factor())
    xi=(Ef2.numerator()).roots()[0][0]
    temp=normOfPoly
    outPoly=0
    for i in range(1+poly.degree()):
        nextCoeff=temp(x=xi)
        outPoly=outPoly+nextCoeff*x^i
        temp=(temp-nextCoeff)/Ef2
    outPoly=outPoly/nextCoeff
    return outPoly

def fourthRoot(PP):
    dd=PP.degree()
    n=dd//4
    rootPoly=x^n
    for i in range(n):
        rootPoly=rootPoly+(x^(n-i-1)/4)*((PP-rootPoly^4)[dd-i-1])
    return rootPoly

def DivideBy2(phi,Wrs):
    FF=phi[0]
    GG=phi[1]
    Fnum=FF.numerator()
    Fden=FF.denominator()
    FnumLead=Fnum[Fnum.degree()]
    FdenLead=Fden[Fden.degree()]
    cF=FnumLead/FdenLead
    PP=(1/FnumLead)*Fnum
    QQ=(1/FdenLead)*(Fden/Wrs).numerator()
    pp=DupPolyRoots(PP,Wrs)
    qq=DupPolyRoots(QQ,Wrs)
    pp0=fourthRoot(pp)
    qq0=fourthRoot(qq)
    ff=4*cF*pp0/qq0
    Ef2b=EndoDouble([x,1],Wrs)[1]
    gg=GG/Ef2b(x=ff)
    return [ff,gg]

def Iso2fromKernel(Wrs,xVal):
    tt=Wrs/(x-xVal)^2
    aa=Wrs[2]
    bb=Wrs[1]
    cc=Wrs[0]
    Qroots=(x^4-2*bb*x^2-8*cc*x+bb^2-4*aa*cc-xVal*4*Wrs).roots()
    x1=Qroots[0][0]
    x2=Qroots[1][0]
    bet1=tt(x=x1)
    bet2=tt(x=x2)
    Wrs2=x*(x-bet1)*(x-bet2)

    EE=EllipticCurve(K,[0,Wrs2[2],0,Wrs2[1],Wrs2[0]])
    jj=EE.j_invariant()
    phi1=[tt,(x-x1)*(x-x2)/(x-xVal)^2]
    XX=(x-bet1)*(x-bet2)/x
    YY=(x^2-bet1*bet2)/x^2
    phi2=[K(1)/K(4)*XX+xVal,K(1)/K(8)*YY]
    return [jj,Wrs2,phi1,phi2]

def smaller_field(phi):
    phi = [[phi[0].numerator(), phi[0].denominator()],[phi[1].numerator(), phi[1].denominator()]]
    phi_small = []
    for i in range(len(phi)):
        lst = []
        for j in range(len(phi[0])):
            coeff = [K(a[0]) for a in phi[i][j]]
            #print(coeff[2])
            deg = len(coeff)
            #RR.<x> = PolynomialRing(Ksub)
            RR.<x> = PolynomialRing(K)
            poly = RR(sum((coeff[i])*x^(deg - i) for i in range(deg)))
            lst.append(poly)
        phi_small.append(lst)
    return phi_small









