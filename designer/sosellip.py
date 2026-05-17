import math
import numpy as np
import sympy as sp
from scipy import signal


def ba2svf2_coeffs(b, a):
    """
    Convert one ordinary biquad section:

        H(q) = (b0 + b1*q + b2*q^2)
             / (1  + a1*q + a2*q^2)

    into your SVFKernel2order coefficients:

        x = in - z1 - z2
        out = d0*x + d1*z1 + d2*z2
        z2 += c2*z1
        z1 += c1*x

    Returns:
        [d0, d1, d2, c1, c2]
    """
    b = np.asarray(b, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)

    b = b / a[0]
    a = a / a[0]

    b0, b1, b2 = b
    _, a1, a2 = a

    c1 = a1 + 2.0
    c2 = (a1 + a2 + 1.0) / c1

    d0 = b0
    d1 = (2.0 * b0 + b1) / c1
    d2 = (b0 + b1 + b2) / (a1 + a2 + 1.0)

    return np.array([d0, d1, d2, c1, c2], dtype=np.float64)


def design_elliptic_svf_sos(
    num_sos=3,
    sample_rate=48000.0,
    cutoff=10000.0,
    rp=0.1,
    rs=60.0,
):
    """
    Design a 2*num_sos order elliptic lowpass, then convert each SOS
    to [d0, d1, d2, c1, c2].
    """
    order = 2 * num_sos

    sos = signal.ellip(
        N=order,
        rp=rp,
        rs=rs,
        Wn=cutoff,
        btype="lowpass",
        analog=False,
        output="sos",
        fs=sample_rate,
    )

    coeffs = []

    for section in sos:
        b = section[0:3]
        a = section[3:6]

        coeffs.append(ba2svf2_coeffs(b, a))

    return np.asarray(coeffs, dtype=np.float64), sos


def format_cpp_2d_array(name, values, c_type="float"):
    rows = []

    for row in values:
        body = ", ".join(f"{v:.9g}f" for v in row)
        rows.append(f"        {{ {body} }}")

    return (
        f"constexpr static {c_type} {name}"
        f"[{len(values)}][{len(values[0])}] = {{\n"
        + ",\n".join(rows)
        + "\n    };"
    )


def pow_to_mul(expr, max_power=4):
    return expr.replace(
        lambda x: (
            isinstance(x, sp.Pow)
            and x.exp.is_Integer
            and 1 < x.exp <= max_power
        ),
        lambda x: sp.Mul(*([x.base] * int(x.exp)), evaluate=False),
    )


def generate_svf_apply_apf_cpp(scalar_type="float"):
    """
    Generate C++ code for:

        void SVFApplyAPF(float* dst, const float* src, float k)

    src layout:
        src[0] = ed0
        src[1] = ed1
        src[2] = ed2
        src[3] = ec1
        src[4] = ec2

    dst layout:
        dst[0] = d0
        dst[1] = d1
        dst[2] = d2
        dst[3] = c1
        dst[4] = c2
    """
    k = sp.Symbol("k")

    ed0 = sp.Symbol("src[0]")
    ed1 = sp.Symbol("src[1]")
    ed2 = sp.Symbol("src[2]")
    ec1 = sp.Symbol("src[3]")
    ec2 = sp.Symbol("src[4]")

    # Original section:
    #
    # H(t) = (ed0 + ed1*p1*t + ed2*p2*t^2)
    #      / (1   +     p1*t +     p2*t^2)
    #
    # p1 = ec1
    # p2 = ec1 * ec2
    p1 = ec1
    p2 = ec1 * ec2

    kp = 1 + k
    km = 1 - k

    # APF warped state variable:
    #
    # t = (k + (1+k)*s) / (1-k)
    #
    # After multiplying numerator and denominator by (1-k)^2:
    #
    # Den(s) = den0 + den1*s + den2*s^2
    # Num(s) = num0 + num1*s + num2*s^2
    den0 = km**2 + p1 * k * km + p2 * k**2
    den1 = kp * (p1 * km + 2 * p2 * k)
    den2 = p2 * kp**2

    num0 = ed0 * km**2 + ed1 * p1 * k * km + ed2 * p2 * k**2
    num1 = kp * (ed1 * p1 * km + 2 * ed2 * p2 * k)
    num2 = ed2 * p2 * kp**2

    # Target section:
    #
    # H(s) = (d0 + d1*P1*s + d2*P2*s^2)
    #      / (1  +    P1*s +    P2*s^2)
    #
    # P1 = c1
    # P2 = c1*c2
    d0 = sp.simplify(num0 / den0)
    d1 = sp.simplify(num1 / den1)
    d2 = sp.simplify(num2 / den2)

    c1 = sp.simplify(den1 / den0)
    c2 = sp.simplify(den2 / den1)

    exprs = [d0, d1, d2, c1, c2]

    temps, reduced = sp.cse(
        exprs,
        symbols=sp.numbered_symbols("t"),
        optimizations=None,
        order="canonical",
    )

    lines = []
    lines.append(f"void SVFApplyAPF({scalar_type}* dst, const {scalar_type}* src, {scalar_type} k)")
    lines.append("{")

    for sym, expr in temps:
        expr = pow_to_mul(expr)
        lines.append(f"    const {scalar_type} {sym} = {sp.ccode(expr)};")

    for i, expr in enumerate(reduced):
        expr = pow_to_mul(expr)
        lines.append(f"    dst[{i}] = {sp.ccode(expr)};")

    lines.append("}")

    return "\n".join(lines)


def main():
    sample_rate = 48000.0
    ellip_cutoff = 10000.0

    # 2 * num_sos order elliptic filter.
    # num_sos = 3 -> 6th order.
    num_sos = 3

    rp = 0.1
    rs = 60.0

    coeffs, sos = design_elliptic_svf_sos(
        num_sos=num_sos,
        sample_rate=sample_rate,
        cutoff=ellip_cutoff,
        rp=rp,
        rs=rs,
    )

    np.set_printoptions(precision=17, suppress=False)

    print("// ===== SOS from scipy =====")
    print(repr(sos))
    print()

    print("// ===== SVF section coefficients: d0,d1,d2,c1,c2 =====")
    print(format_cpp_2d_array("ellpCoeffs", coeffs, c_type="float"))
    print()

    print("// ===== Generated APF application function =====")
    print(generate_svf_apply_apf_cpp(scalar_type="float"))


if __name__ == "__main__":
    main()