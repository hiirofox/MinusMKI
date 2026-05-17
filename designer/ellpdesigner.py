import numpy as np
from scipy import signal


def format_c_array(name, values, c_type="double"):
    body = ", ".join(f"{v:.17g}" for v in values)
    return f"const {c_type} {name}[{len(values)}] = {{ {body} }};"


def main():
    # ===== Filter spec =====
    fs = 48000.0
    fc = 10000.0

    # 改这里即可
    order = 6

    rp = 0.1     # passband ripple, dB
    rs = 60.0    # stopband attenuation, positive dB value

    # ===== Design elliptic lowpass IIR =====
    b, a = signal.ellip(
        N=order,
        rp=rp,
        rs=rs,
        Wn=fc,
        btype="lowpass",
        analog=False,
        output="ba",
        fs=fs,
    )

    # 保险起见，归一化到 a[0] == 1
    b = np.asarray(b, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)

    b = b / a[0]
    a = a / a[0]

    # ===== Python array output =====
    np.set_printoptions(precision=17, suppress=False)

    print("b = np.array(")
    print(repr(b))
    print(", dtype=np.float64)")
    print()

    print("a = np.array(")
    print(repr(a))
    print(", dtype=np.float64)")
    print()

    # ===== C array output =====
    print(format_c_array("b", b))
    print(format_c_array("a", a))


if __name__ == "__main__":
    main()