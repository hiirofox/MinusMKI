import sympy as sp


def iir_to_svf_symbolic_ap(n):
    u, q, r, k = sp.symbols("u q r k")

    a = sp.symbols(f"a1:{n+1}")
    b = sp.symbols(f"b0:{n+1}")

    A = 1 + sum(a[i - 1] * u**i for i in range(1, n + 1))
    B = sum(b[i] * u**i for i in range(n + 1))

    # 一阶全通保角变换:
    # u <- (k + q) / (1 + k q)
    u_ap = (k + q) / (1 + k*q)

    # SVF 积分器变量:
    # r = 1 / (1 - q)
    # q = 1 - 1/r
    q_svf = 1 - 1/r

    # 合成变换
    u_r = sp.factor(sp.together(u_ap.subs(q, q_svf)))

    # 取 u_r 的分母
    u_num, u_den = sp.fraction(u_r)

    # 关键修正：
    # 先 together，再取 numerator，确保 Ar / Br 是关于 r 的多项式
    Ar = sp.factor(sp.together(u_den**n * A.subs(u, u_r)).as_numer_denom()[0])
    Br = sp.factor(sp.together(u_den**n * B.subs(u, u_r)).as_numer_denom()[0])

    Ar_poly = sp.Poly(Ar, r)
    Br_poly = sp.Poly(Br, r)

    # 分母常数项，用来归一化
    K = Ar_poly.coeff_monomial(r**0)

    den = sp.expand(Ar / K)
    num = sp.expand(Br / K)

    den_poly = sp.Poly(den, r)
    num_poly = sp.Poly(num, r)

    c = [
        sp.factor(sp.simplify(-den_poly.coeff_monomial(r**i)))
        for i in range(1, n + 1)
    ]

    d = [
        sp.factor(sp.simplify(num_poly.coeff_monomial(r**i)))
        for i in range(n + 1)
    ]

    return {
        "c": c,
        "d": d,
        "u_r": u_r,
        "K": K,
    }

def iir_to_zdf_state_symbolic_ap(n):
    u, s, k = sp.symbols("u s k")

    a = sp.symbols(f"a1:{n+1}")
    b = sp.symbols(f"b0:{n+1}")

    A = 1 + sum(a[i - 1] * u**i for i in range(1, n + 1))
    B = sum(b[i] * u**i for i in range(n + 1))

    q_s = s / (1 + s)

    # 一阶全通保角变换:
    # u <- (k + q) / (1 + k*q)
    u_s = sp.factor(sp.together((k + q_s) / (1 + k*q_s)))

    u_num, u_den = sp.fraction(u_s)

    As = sp.factor(sp.together(u_den**n * A.subs(u, u_s)).as_numer_denom()[0])
    Bs = sp.factor(sp.together(u_den**n * B.subs(u, u_s)).as_numer_denom()[0])

    As_poly = sp.Poly(As, s)
    Bs_poly = sp.Poly(Bs, s)

    # 全通版本下，常数项不是 1，而是 A(k)，所以要归一化
    K = As_poly.coeff_monomial(s**0)

    den = sp.expand(As / K)
    num = sp.expand(Bs / K)

    den_poly = sp.Poly(den, s)
    num_poly = sp.Poly(num, s)

    p = [sp.Integer(1)]
    for i in range(1, n + 1):
        p.append(sp.factor(sp.simplify(den_poly.coeff_monomial(s**i))))

    c = [
        sp.factor(sp.simplify(p[i] / p[i - 1]))
        for i in range(1, n + 1)
    ]

    d = [
        sp.factor(sp.simplify(num_poly.coeff_monomial(s**i) / p[i]))
        for i in range(n + 1)
    ]

    return {
        "c": c,
        "d": d,
        "p": p[1:],
        "K": K,
    }



def emit_c_coeff_code(result, n, scalar_type="double"):
    # 顺序：d0, d1, ..., dn, c1, ..., cn
    names = (
        [f"d{i}" for i in range(n + 1)] +
        [f"c{i}" for i in range(1, n + 1)]
    )

    exprs = result["d"] + result["c"]

    # 公共子表达式消除
    temps, reduced = sp.cse(
        exprs,
        symbols=sp.numbered_symbols("t"),
        optimizations="basic",
        order="canonical",
    )

    lines = []

    for sym, expr in temps:
        lines.append(f"{scalar_type} {sym} = {sp.ccode(expr)};")

    for name, expr in zip(names, reduced):
        lines.append(f"{name} = {sp.ccode(expr)};")

    return "\n".join(lines)

stages=2

result = iir_to_zdf_state_symbolic_ap(stages)

print(emit_c_coeff_code(result, stages))