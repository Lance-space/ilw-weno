"""Profile figures for the nine examples (robust to E-less exponent output)."""
import os, sys, re, math
import numpy as np
sys.path.insert(0, r"C:\Users\w10\Documents\deepseek-harness\default-workspace\ilw_weno\tools")
from pdfplot import Canvas

B = r"C:\Users\w10\Documents\deepseek-harness\default-workspace\ilw_weno"
FIG = os.path.join(B, "paper", "figs")
NUM = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+(?:[EeDd][-+]?\d+)?|[-+]?\d*\.\d+[-+]\d+")

CASES = [("ex6_1", "1__heat_1D", 1, True), ("ex6_2", "2__Barenblatt_m=2", 1, True),
         ("ex6_3", "3__two_boxes", 1, False), ("ex6_4", "4__BL_IBVP", 1, False),
         ("ex6_5", "5__degen_1D", 1, False), ("ex6_6", "6__heat_2D", 2, True),
         ("ex6_7", "7__PME_2D_gauss", 2, False), ("ex6_8", "8__BL_2D", 2, False),
         ("ex6_9", "9__degen_2D", 2, False)]


def nums(line):
    out = []
    for tok in line.replace(",", " ").split():
        m = re.fullmatch(r"([-+]?\d*\.?\d+)([-+])(\d+)", tok)   # 1.234-102
        if m and "e" not in tok.lower():
            out.append(float(m.group(1) + "e" + m.group(2) + m.group(3)))
            continue
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out


def load(tag):
    rows = []
    started = False
    for line in open(os.path.join(B, "runs", tag, "profile.plt")):
        if line.startswith(" ZONE") or line.startswith(" VARIABLES") or line.startswith(" TITLE"):
            started = True
            continue
        if not started:
            continue
        v = nums(line)
        if v:
            rows.append(v)
    return rows


for name, tag, dim, exact in CASES:
    rk = load(f"nine_{tag}_RK3-WENO6")
    il = load(f"nine_{tag}_ILW3-WENO6")
    if dim == 1:
        x = np.array([r[0] for r in rk])
        u_rk = np.array([r[1] for r in rk])
        u_ex = np.array([r[2] for r in rk]) if exact else None
        u_il = np.array([r[1] for r in il])
    else:
        n = int(round(math.sqrt(len(rk))))
        x = np.array([r[0] for r in rk]).reshape(n, n)[n // 2, :]
        u_rk = np.array([r[2] for r in rk]).reshape(n, n)[n // 2, :]
        u_ex = np.array([r[3] for r in rk]).reshape(n, n)[n // 2, :] if exact else None
        u_il = np.array([r[2] for r in il]).reshape(n, n)[n // 2, :]
    ys = list(u_rk) + list(u_il) + (list(u_ex) if exact else [])
    ymin, ymax = float(min(ys)), float(max(ys))
    if ymax - ymin < 1e-12:
        ymax = ymin + 1.0
    c = Canvas(300, 190)
    c.axes(46, 34, 275, 158)
    c.data_range(float(min(x)), float(max(x)), ymin - 0.06 * (ymax - ymin), ymax + 0.12 * (ymax - ymin))
    c.frame()
    if exact:
        c.polyline(list(zip(x, u_ex)), lw=1.1)
    c.polyline(list(zip(x, u_rk)), lw=0.8, dash=(4, 2))
    c.polyline(list(zip(x, u_il)), lw=0.8, dash=(1, 2))
    c.text((min(x) + max(x)) / 2, 18, "x", 8, "center")
    c.text(15, 96, "u", 8, "center", angle=90)
    yy = 150
    if exact:
        c.fig_line(205, yy, 228, yy, lw=1.1); c.text(231, yy - 2.4, "exact", 6.5); yy -= 11
    c.fig_line(205, yy, 228, yy, lw=0.8, dash=(4, 2)); c.text(231, yy - 2.4, "RK3", 6.5); yy -= 11
    c.fig_line(205, yy, 228, yy, lw=0.8, dash=(1, 2)); c.text(231, yy - 2.4, "ILW3", 6.5)
    c.save(os.path.join(FIG, name + ".pdf"))
    print("%-7s n=%d  u in [%.3g, %.3g]" % (name, len(x), ymin, ymax))
