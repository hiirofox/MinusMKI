import sympy as sp


def iir_to_svf_symbolic(n):
    q, r = sp.symbols("q r")

    # 普通 IIR 系数:
    # A(q) = 1 + a1 q + ... + an q^n
    # B(q) = b0 + b1 q + ... + bn q^n
    a = sp.symbols(f"a1:{n+1}")
    b = sp.symbols(f"b0:{n+1}")

    A = 1 + sum(a[i - 1] * q**i for i in range(1, n + 1))
    B = sum(b[i] * q**i for i in range(n + 1))

    # r = Hi = 1 / (1 - q), therefore q = 1 - 1/r
    Ar = sp.expand(r**n * A.subs(q, 1 - 1/r))
    Br = sp.expand(r**n * B.subs(q, 1 - 1/r))

    # 分母常数项，用来归一化到 1 - c1 r - ... - cn r^n
    K = sp.Poly(Ar, r).coeff_monomial(r**0)

    den = sp.expand(Ar / K)
    num = sp.expand(Br / K)

    # den = 1 - c1 r - c2 r^2 - ...
    c = [
        sp.simplify(-sp.Poly(den, r).coeff_monomial(r**k))
        for k in range(1, n + 1)
    ]

    # num = d0 + d1 r + ... + dn r^n
    d = [
        sp.simplify(sp.Poly(num, r).coeff_monomial(r**k))
        for k in range(n + 1)
    ]

    return {
        "A(q)": A,
        "B(q)": B,
        "Ar(r)": Ar,
        "Br(r)": Br,
        "K": K,
        "c": c,
        "d": d,
    }


result = iir_to_svf_symbolic(2)

for i, expr in enumerate(result["d"]):
    print(f"d{i} = {sp.simplify(expr)}")

for i, expr in enumerate(result["c"], start=1):
    print(f"c{i} = {sp.simplify(expr)}")