def IsoTriple(phi,Wrs1,Wrs2):
    return IsoAdd(IsoDouble(phi,Wrs1,Wrs2),phi,Wrs1,Wrs2)

#The function Trip Conjugates assumes Wrs is in *short weierstrass form* (i.e., x^3 + ax + b)
def TripConjugates(Wrs):
    rrfunction = 3*x^4 + 6*Wrs[1]*x^2 + 12*Wrs[0]*x - Wrs[1]^2
    rr = rrfunction.roots(multiplicities = False)
    a2 = Wrs[2]
    conjList=[t]
    yy = [y0,y1,y2,y3]
    for i in range(4):
        conjList.append((yy[i]-Y)^2/(t - rr[i])^2 - a2 - rr[i]-t)
        conjList.append((yy[i]+Y)^2/(t - rr[i])^2 - a2 - rr[i]-t)
    return conjList,rr #Returns a list of x-coordinates of points: x + (3-torsion pt), for half othe nontrivial 3-torsion points.
#Since this only considers x-coordinates, we don't actually add all of the 3-torsionpoints- Just four of them- not their negatives.

def SubstituteYsquare(Yvar,polyn,Wrs):
    Wrs = sum(Wrs[i]*t^i for i in range(4))
    maxYdeg = polyn.degree(Yvar)
    newpoly = 0
    for j in range(maxYdeg+1):
        coeff_now = polyn.coefficient({Yvar:j})
        if ZZ(j).mod(2) == 0:
            newpoly = newpoly+coeff_now*(Wrs)^(ZZ(ZZ(j)/2))
        if ZZ(j).mod(2) != 0 and coeff_now !=0:
            print("ERROR")
    return newpoly

def Substituteyisquare(Yvar,polyn,Wrs):
    maxYdeg = polyn.degree(Yvar)
    newpoly = 0
    for j in range(maxYdeg+1):
        coeff_now = polyn.coefficient({Yvar:j})
        if ZZ(j).mod(2) == 0:
            newpoly = newpoly+coeff_now*(Wrs)^(ZZ(ZZ(j)/2))
        if ZZ(j).mod(2) != 0 and coeff_now !=0:
            print("ERROR")
    return newpoly

def TripPolyRoots(poly,Wrs):
    degpoly = poly.degree()
    Ef3 = EndoTriple([x,1],Wrs)[0]
    ConjListAndrr = TripConjugates(Wrs)
    TripConjList = ConjListAndrr[0]
    rrconj = ConjListAndrr[1]
    ylist = [y0,y1,y2,y3]
    normOfPoly = poly
    temppolynlist = []
    temppolydlist = []
    for i in range(4):
        temppoly = poly(x = TripConjList[2*i+1])*poly(x = TripConjList[2*i+2])
        temppolyn1 = temppoly.numerator()
        temppolyd1 = temppoly.denominator()
        temppolynlist.append(SubstituteYsquare(Y,temppolyn1,Wrs))
        temppolydlist.append(SubstituteYsquare(Y,temppolyd1,Wrs))
    for i in range(4):
        for j in range(len(rrconj)):
            temppolynlist[i] = Substituteyisquare(ylist[j],temppolynlist[i],Wrs(x=rrconj[j]))
            temppolydlist[i] = Substituteyisquare(ylist[j],temppolydlist[i],Wrs(x=rrconj[j]))
    num = prod(temppolynlist)
    den = prod(temppolydlist)
    polyt = poly(x=t)
    norm = polyt*num/den
    normOfPoly = norm(t=x)
    xi = (Ef3.numerator()).roots()[0][0]
    temp = normOfPoly
    outPoly = 0
    for i in range(1+poly.degree()):
        nextCoeff=temp(x=xi)
        outPoly=outPoly+nextCoeff*x^i
        temp=(temp-nextCoeff)/Ef3
    outPoly=outPoly/nextCoeff
    return outPoly

def TripPolyRoots_comments_11_1(poly,Wrs):
    degpoly = poly.degree()
    Ef3 = EndoTriple([x,1],Wrs)[0]
    ConjListAndrr = TripConjugates(Wrs)
    TripConjList = ConjListAndrr[0]
    rrconj = ConjListAndrr[1]
    ylist = [y0,y1,y2,y3]
    print("STEP 4:")
    normOfPoly = poly
    print("N(x) <- ",normOfPoly)
    print("D(x) <- ","1")
    temppolynlist = []
    temppolydlist = []
    print("STEP 5:")
    print("Begin the For loop for i = 1,...,4:")
    for i in range(4):
        temppoly = poly(x = TripConjList[2*i+1])*poly(x = TripConjList[2*i+2])
        temppolyn1 = temppoly.numerator()
        temppolyd1 = temppoly.denominator()
        temppolynlist.append(SubstituteYsquare(Y,temppolyn1,Wrs))
        temppolydlist.append(SubstituteYsquare(Y,temppolyd1,Wrs))
        print("STEP 10 of For loop:")
        print("N_i(x) for i = ",i,": ",temppolynlist[i-1])
        print("D_i(x) for i = ",i,": ",temppolydlist[i-1])
    for i in range(4):
        for j in range(len(rrconj)):
            temppolynlist[i] = Substituteyisquare(ylist[j],temppolynlist[i],Wrs(x=rrconj[j]))
            temppolydlist[i] = Substituteyisquare(ylist[j],temppolydlist[i],Wrs(x=rrconj[j]))
    num = prod(temppolynlist)
    den = prod(temppolydlist)
    polyt = poly(x=t)
    norm = polyt*num/den
    print("STEP 11:")
    print("N_P(x) = ", norm)
    normOfPoly = norm(t=x)
    xi = (Ef3.numerator()).roots()[0][0]
    temp = normOfPoly
    outPoly = 0
    print("STEP 12:")
    print("Begin the second For loop:")
    for i in range(1+poly.degree()):
        nextCoeff=temp(x=xi)
        print("NextCoeff = ",nextCoeff)
        outPoly=outPoly+nextCoeff*x^i
        print("STEP 14 of For loop:")
        print("p(x) = ",outPoly)
        temp=(temp-nextCoeff)/Ef3
    outPoly=outPoly/nextCoeff
    print("STEP 17:")
    print("p(x) = ",outPoly)
    return outPoly

def ninthRoot(PP):
    dd = PP.degree()
    if dd == 0:
        return 1
    n = dd//9
    rootPoly = x^n
    for i in range(n):
        rootPoly = rootPoly + (x^(n - i - 1)/9)*((PP - rootPoly^9)[dd-i-1])
    if rootPoly^9==PP:
        return rootPoly
    else:
        print("error, not a 9th power")
        return False

def DivideBy3(phi,Wrs):
    FF = phi[0]
    GG = phi[1]
    b2 = Wrs[2]*4
    b4 = 2*Wrs[1]
    b6 = 4*Wrs[0]
    b8 = 4*Wrs[2]*Wrs[0]-Wrs[1]^2
    ThreeDivPoly = (x^4 + b2/3*x^3 + b4*x^2 + b6*x + b8/3)
    Fnum=FF.numerator()
    Fden=FF.denominator()
    FnumLead=Fnum[Fnum.degree()]
    FdenLead=Fden[Fden.degree()]
    cF=FnumLead/FdenLead
    PP=(1/FnumLead)*Fnum
    QQ=(1/FdenLead)*(Fden/(ThreeDivPoly^2)).numerator()
    MultBy3 = EndoTriple([x,1],Wrs)
    PProots = PP.roots()
    pp=TripPolyRoots(PP,Wrs)
    if QQ.degree() == 0:
        qq = 1
    else:
        qq=TripPolyRoots(QQ,Wrs)
    pp0 = ninthRoot(pp)
    if qq == 1:
        qq0 = 1
    else: 
        qq0 = ninthRoot(qq)
    ff = 9*cF*pp0/qq0
    Ef3b = MultBy3[1]
    gg = GG/Ef3b(x = ff)
    return [ff,gg]

def DivideBy3_comments11_1(phi,Wrs):
    FF = phi[0]
    GG = phi[1]
    b2 = Wrs[2]*4
    b4 = 2*Wrs[1]
    b6 = 4*Wrs[0]
    b8 = 4*Wrs[2]*Wrs[0]-Wrs[1]^2
    ThreeDivPoly = (x^4 + b2/3*x^3 + b4*x^2 + b6*x + b8/3)
    Fnum=FF.numerator()
    Fden=FF.denominator()
    FnumLead=Fnum[Fnum.degree()]
    FdenLead=Fden[Fden.degree()]
    cF=FnumLead/FdenLead
    PP=(1/FnumLead)*Fnum
    QQ=(1/FdenLead)*(Fden/(ThreeDivPoly^2)).numerator()
    MultBy3 = EndoTriple([x,1],Wrs)
    PProots = PP.roots()
    pp=TripPolyRoots_comments_11_1(PP,Wrs)
    if QQ.degree() == 0:
        qq = 1
    else:
        qq=TripPolyRoots(QQ,Wrs)
    pp0 = ninthRoot(pp)
    if qq == 1:
        qq0 = 1
    else: 
        qq0 = ninthRoot(qq)
    ff = 9*cF*pp0/qq0
    Ef3b = MultBy3[1]
    gg = GG/Ef3b(x = ff)
    return [ff,gg]

def DivideBy3_comments(phi,Wrs):
    FF = phi[0]
    GG = phi[1]
    b2 = Wrs[2]*4
    b4 = 2*Wrs[1]
    b6 = 4*Wrs[0]
    b8 = 4*Wrs[2]*Wrs[0]-Wrs[1]^2
    ThreeDivPoly = (x^4 + b2/3*x^3 + b4*x^2 + b6*x + b8/3)
    Fnum=FF.numerator()
    Fden=FF.denominator()
    FnumLead=Fnum[Fnum.degree()]
    FdenLead=Fden[Fden.degree()]
    cF=FnumLead/FdenLead
    PP=(1/FnumLead)*Fnum
    QQ=(1/FdenLead)*(Fden/(ThreeDivPoly^2)).numerator()
    print("STEP 1:")
    print("c_F = ",cF)
    print("P(x) = ",PP)
    print("Q(x) = ",QQ)
    MultBy3 = EndoTriple([x,1],Wrs)
    print("STEP 2:")
    print("Psi_(E_1728,3)(x) = ",MultBy3)
    PProots = PP.roots()
    pp=TripPolyRoots(PP,Wrs)
    print("STEP 3:")
    print("p(x) = ",pp)
    if QQ.degree() == 0:
        qq = 1
    else:
        qq=TripPolyRoots(QQ,Wrs)
    print("STEP 4:")
    print("q(x) = ",qq)
    pp0 = ninthRoot(pp)
    if qq == 1:
        qq0 = 1
    else: 
        qq0 = ninthRoot(qq)
    print("STEP 5:")
    print("p_0(x) = ",pp0)
    print("q_0(x) = ",qq0)
    ff = 9*cF*pp0/qq0
    print("STEP 6:")
    print("f(x) = ",ff)
    Ef3b = MultBy3[1]
    gg = GG/Ef3b(x = ff)
    print("g(x) = ",gg)
    print("STEP 7:")
    return [ff,gg]