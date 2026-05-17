import math
import numpy as np
import sympy as sp
from scipy import signal


def svf2_to_ba(d, c1, c2):
    """
    SVFKernel2order coeffs -> ordinary biquad b,a.

    Kernel:
        x = in - z1 - z2
        out = d0*x + d1*z1 + d2*z2
        z2 += c2*z1
        z1 += c1*x

    Returns:
        b = [b0,b1,b2]
        a = [1,a1,a2]
    """
    d0, d1, d2 = d

    b0 = d0
    b1 = -2.0 * d0 + c1 * d1
    b2 = d0 - c1 * d1 + c1 * c2 * d2

    a1 = c1 - 2.0
    a2 = 1.0 - c1 + c1 * c2

    b = np.array([b0, b1, b2], dtype=np.float64)
    a = np.array([1.0, a1, a2], dtype=np.float64)

    return b, a


def eval_biquad_mag_at(b, a, w):
    """
    Evaluate |B(q)/A(q)| at q = exp(-j*w).
    """
    q = np.exp(-1j * w)

    b0, b1, b2 = b
    a0, a1, a2 = a

    num = b0 + b1 * q + b2 * q * q
    den = a0 + a1 * q + a2 * q * q

    return abs(num / den)


def eval_mode_cascade_mag_at(rows, mode_offset, w):
    """
    rows layout:
        [0]  c1
        [1]  c2
        [2]  lp_d0
        [3]  lp_d1
        [4]  lp_d2
        [5]  bp_d0
        [6]  bp_d1
        [7]  bp_d2
        [8]  hp_d0
        [9]  hp_d1
        [10] hp_d2

    mode_offset:
        LP = 2
        BP = 5
        HP = 8
    """
    mag = 1.0

    for row in rows:
        c1 = row[0]
        c2 = row[1]

        d = row[mode_offset:mode_offset + 3]

        b, a = svf2_to_ba(d, c1, c2)

        mag *= eval_biquad_mag_at(b, a, w)

    return mag


def apply_cutoff_gain_matching(
    rows,
    sample_rate,
    ellip_cutoff,
    target="lp",
    eps=1.0e-12,
):
    """
    Match LP/BP/HP gain at ellip_cutoff.

    target:
        "min"  : safest, only attenuates louder modes to the quietest one.
        "bp"   : match LP/HP to BP at cutoff.
        1.0    : force all modes to unity at cutoff.
        number : force all modes to that magnitude.

    Returns:
        compensated rows, info dict.
    """
    rows = np.asarray(rows, dtype=np.float64).copy()

    num_sos = rows.shape[0]

    w = 2.0 * math.pi * ellip_cutoff / sample_rate

    lp_mag = eval_mode_cascade_mag_at(rows, 2, w)
    bp_mag = eval_mode_cascade_mag_at(rows, 5, w)
    hp_mag = eval_mode_cascade_mag_at(rows, 8, w)

    if target == "min":
        target_mag = min(lp_mag, bp_mag, hp_mag)
    elif target == "lp":
        target_mag = lp_mag
    else:
        target_mag = float(target)

    mode_infos = {
        "lp_mag_before": lp_mag,
        "bp_mag_before": bp_mag,
        "hp_mag_before": hp_mag,
        "target_mag": target_mag,
    }

    modes = [
        ("lp", 2, lp_mag),
        ("bp", 5, bp_mag),
        ("hp", 8, hp_mag),
    ]

    for name, offset, mag in modes:
        if mag < eps:
            section_gain = 1.0
        else:
            total_gain = target_mag / mag
            section_gain = total_gain ** (1.0 / num_sos)

        rows[:, offset:offset + 3] *= section_gain

        mode_infos[f"{name}_section_gain"] = section_gain
        mode_infos[f"{name}_total_gain"] = section_gain ** num_sos

    mode_infos["lp_mag_after"] = eval_mode_cascade_mag_at(rows, 2, w)
    mode_infos["bp_mag_after"] = eval_mode_cascade_mag_at(rows, 5, w)
    mode_infos["hp_mag_after"] = eval_mode_cascade_mag_at(rows, 8, w)

    return rows, mode_infos

def poly_eval_q(coeffs, q):
    """
    coeffs: b0, b1, ..., bn for B(q)
    """
    acc = 0.0 + 0.0j
    p = 1.0 + 0.0j

    for c in coeffs:
        acc += c * p
        p *= q

    return acc


def positive_root(x, n):
    x = float(np.real_if_close(x))

    if x < 0.0:
        raise ValueError(f"negative gain cannot be evenly rooted safely: {x}")

    return x ** (1.0 / n)


def b_to_svf_d(b, a):
    """
    Convert one biquad:

        H(q) = (b0 + b1*q + b2*q^2)
             / (1  + a1*q + a2*q^2)

    to:

        x = in - z1 - z2
        out = d0*x + d1*z1 + d2*z2
        z2 += c2*z1
        z1 += c1*x

    Returns d0,d1,d2,c1,c2.
    """
    b = np.asarray(b, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)

    b = b / a[0]
    a = a / a[0]

    b0, b1, b2 = b
    _, a1, a2 = a

    c1 = a1 + 2.0
    c2 = (1.0 + a1 + a2) / c1

    p2 = c1 * c2

    d0 = b0
    d1 = (b1 + 2.0 * b0) / c1
    d2 = (b0 + b1 + b2) / p2

    return np.array([d0, d1, d2, c1, c2], dtype=np.float64)


def design_same_pole_multimode_svf_sos(
    num_sos=3,
    sample_rate=48000.0,
    ellip_cutoff=10000.0,
    rp=0.1,
    rs=60.0,
):
    """
    Generate a 2*num_sos order fixed-pole multimode filter.

    Output layout per SOS section:

        [0]  c1
        [1]  c2

        [2]  lp_d0
        [3]  lp_d1
        [4]  lp_d2

        [5]  bp_d0
        [6]  bp_d1
        [7]  bp_d2

        [8]  hp_d0
        [9]  hp_d1
        [10] hp_d2
    """
    order = 2 * num_sos

    # We use scipy elliptic only to get the denominator / pole placement.
    sos = signal.ellip(
        N=order,
        rp=rp,
        rs=rs,
        Wn=ellip_cutoff,
        btype="lowpass",
        analog=False,
        output="sos",
        fs=sample_rate,
    )

    den_sections = []
    a_global = np.array([1.0], dtype=np.float64)

    for sec in sos:
        a = np.asarray(sec[3:6], dtype=np.float64)
        a = a / a[0]

        den_sections.append(a)
        a_global = np.convolve(a_global, a)

    # Global canonical mode numerators.
    #
    # LP shape:
    #   B_lp(q) = g_lp * (1 + q)^order
    #
    # HP shape:
    #   B_hp(q) = g_hp * (1 - q)^order
    #
    # BP shape:
    #   B_bp(q) = g_bp * (1 - q^2)^num_sos

    A_dc = poly_eval_q(a_global, 1.0)
    A_nyq = poly_eval_q(a_global, -1.0)

    lp_total_gain = np.real_if_close(A_dc) / (2.0 ** order)
    hp_total_gain = np.real_if_close(A_nyq) / (2.0 ** order)

    lp_g = positive_root(lp_total_gain, num_sos)
    hp_g = positive_root(hp_total_gain, num_sos)

    # BP normalize around ellip_cutoff.
    w0 = 2.0 * math.pi * ellip_cutoff / sample_rate
    q0 = np.exp(-1j * w0)

    A_w0 = poly_eval_q(a_global, q0)
    BP_shape_w0 = (1.0 - q0 * q0) ** num_sos

    bp_total_gain = abs(A_w0) / max(abs(BP_shape_w0), 1.0e-30)
    bp_g = bp_total_gain ** (1.0 / num_sos)

    # Same numerator section shapes for all sections.
    b_lp_sec = np.array([lp_g, 2.0 * lp_g, lp_g], dtype=np.float64)
    b_bp_sec = np.array([bp_g, 0.0, -bp_g], dtype=np.float64)
    b_hp_sec = np.array([hp_g, -2.0 * hp_g, hp_g], dtype=np.float64)

    rows = []

    for a in den_sections:
        lp = b_to_svf_d(b_lp_sec, a)
        bp = b_to_svf_d(b_bp_sec, a)
        hp = b_to_svf_d(b_hp_sec, a)

        c1 = lp[3]
        c2 = lp[4]

        # Sanity check: all modes must share c1/c2.
        if abs(bp[3] - c1) > 1e-10 or abs(hp[3] - c1) > 1e-10:
            raise RuntimeError("c1 mismatch")
        if abs(bp[4] - c2) > 1e-10 or abs(hp[4] - c2) > 1e-10:
            raise RuntimeError("c2 mismatch")

        rows.append([
            c1, c2,
            lp[0], lp[1], lp[2],
            bp[0], bp[1], bp[2],
            hp[0], hp[1], hp[2],
            ])
        
    rows = np.asarray(rows, dtype=np.float64)

    rows, gain_info = apply_cutoff_gain_matching(
        rows,
        sample_rate=sample_rate,
        ellip_cutoff=ellip_cutoff,
        target="min",
    )

    print("// cutoff gain matching")
    for key, value in gain_info.items():
        print(f"// {key}: {value}")

    return rows, sos, a_global


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
    Generate:

        void SVFApplyAPF(float* dst, const float* src, float k)

    src layout:
        src[0] = d0
        src[1] = d1
        src[2] = d2
        src[3] = c1
        src[4] = c2

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

    p1 = ec1
    p2 = ec1 * ec2

    kp = 1 + k
    km = 1 - k

    den0 = km**2 + p1 * k * km + p2 * k**2
    den1 = kp * (p1 * km + 2 * p2 * k)
    den2 = p2 * kp**2

    num0 = ed0 * km**2 + ed1 * p1 * k * km + ed2 * p2 * k**2
    num1 = kp * (ed1 * p1 * km + 2 * ed2 * p2 * k)
    num2 = ed2 * p2 * kp**2

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

    num_sos = 3
    rp = 0.1
    rs = 60.0

    coeffs, sos, a_global = design_same_pole_multimode_svf_sos(
        num_sos=num_sos,
        sample_rate=sample_rate,
        ellip_cutoff=ellip_cutoff,
        rp=rp,
        rs=rs,
    )

    np.set_printoptions(precision=17, suppress=False)

    print("// ===== generated fixed-pole multimode coefficients =====")
    print("// layout per row:")
    print("// c1,c2, lp_d0,lp_d1,lp_d2, bp_d0,bp_d1,bp_d2, hp_d0,hp_d1,hp_d2")
    print(format_cpp_2d_array("ellpCoeffs", coeffs, c_type="float"))
    print()

    print("// ===== optional debug: scipy SOS =====")
    print(repr(sos))
    print()

    print("// ===== generated APF warp function =====")
    print(generate_svf_apply_apf_cpp(scalar_type="float"))


if __name__ == "__main__":
    main()