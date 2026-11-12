import os, sys, json
import numpy as np
sys.path.insert(0, r"C:\Users\w10\Documents\deepseek-harness\default-workspace\ilw_weno\tools")
from figex import load

B = r"C:\Users\w10\Documents\deepseek-harness\default-workspace\ilw_weno"
P = os.path.join(B, "paper")
CASES = [("3__two_boxes", 12.0 / 160, 1), ("4__BL_IBVP", 1.0 / 400, 1),
         ("5__degen_1D", 4.0 / 400, 1), ("7__PME_2D_gauss", 20.0 / 80, 2),
         ("8__BL_2D", 3.0 / 120, 2), ("9__degen_2D", 3.0 / 120, 2)]
diff = {}
for tag, dx, dim in CASES:
    col = 1 if dim == 1 else 2
    rk3 = np.array([r[col] for r in load("nine_%s_RK3-WENO6" % tag)])
    il3 = np.array([r[col] for r in load("nine_%s_ILW3-WENO6" % tag)])
    rk4 = np.array([r[col] for r in load("nine_%s_RK4-WENO8" % tag)])
    il4 = np.array([r[col] for r in load("nine_%s_ILW4-WENO8" % tag)])
    cell = dx * dx if dim == 2 else dx
    d3, d4 = abs(il3 - rk3), abs(il4 - rk4)
    diff[tag] = (float(d3.sum() * cell), float(d3.max()), float(d4.sum() * cell), float(d4.max()))
    print("%-16s ILW3-RK3 %.2e %.2e | ILW4-RK4 %.2e %.2e" % ((tag,) + diff[tag]))
json.dump(diff, open(os.path.join(B, "diff_results.json"), "w"), indent=1)
