'''
This is the code to calculate the analytical form of P_n(x), Q^(+/-)_n(x)
with 1 <= x <= 2 only, and calculate the density bound with x=2 up to given n.
'''


from sympy import symbols, Integral, Eq, integrate, expand
import math
x, c, c2 = symbols('x c TMP')

def constexpr(f):
    expr = c2 * 1
    return expr.subs(c2, f)

maxn = 10


ps = [2 / c - 1]
# calculate P_n(x), x<=2
for i in range(1, maxn + 1):
    lf = ps[i-1]
    newf = integrate(lf, (c, 1, c2)).subs(c2, c) * (2 / c)
    print("P", i, expand(newf))
    ps.append(newf)

# calculate Q_n^+(x), x<=2
qp = [constexpr(0)]
for i in range(1, maxn + 1):
    lf = qp[i-1]
    newf = integrate(lf, (c, 1, c2)).subs(c2, c)
    newf = (newf + ps[i-1]) / c
    print("Qplus", i, expand(newf))
    qp.append(newf)

# calculate Q_n^-(x), 1<x<2
qm = [constexpr(0)]
for i in range(1, maxn + 1):
    lf = qm[i-1]
    newf = integrate(lf, (c, 1, c2)).subs(c2, c)
    if (i == 1):
        newf = (newf + 1) / c
    else:
        newf = newf / c
    #  newf = (newf + ps[i-1].subs(c, c-1)) / c
    print("Qminus-1", i, expand(newf))
    qm.append(newf)

# calculate Q_n^-(x), x = 2
qm2 = [constexpr(0), constexpr(0)]
for i in range(2, maxn + 1):
    lf = qm[i-1]
    newf = integrate(lf, (c, 1, c2)).subs(c2, c)
    if (i == 2):
        newf = (newf + 1) / c
    else:
        newf = newf / c
    print("Qminus-2", i, expand(newf))
    qm2.append(newf)

# calculate density bound
for i in range(1, maxn + 1):
    ans = 0
    resp = 1
    resq = 1
    for j in range(2, i + 1):
        #  pv = qp[j].evalf(subs={c: 2})
        pv = qp[j].subs(c, 2)
        #  qv = qm2[j].evalf(subs={c: 2})
        qv = qm2[j].subs(c, 2)
        print(pv, qv)
        ans += (pv + qv) / j / 2
        resp -= pv
        resq -= qv
    ans += (resp + resq) / (i + 1) / 2
    if i > 1:
        print(i, ans, ans.evalf())

