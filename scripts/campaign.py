import os, subprocess, re, sys, shutil, json, time

BASE = r"C:\Users\w10\Documents\deepseek-harness\default-workspace\ilw_weno"
EXE  = os.path.join(BASE, "build", "ilw", "ilw.exe")
RUNS = os.path.join(BASE, "runs")

def run_case(tag, nml, keep=True):
    d = os.path.join(RUNS, tag)
    os.makedirs(d, exist_ok=True)
    shutil.copy(EXE, os.path.join(d, "ilw.exe"))
    with open(os.path.join(d, "input.nml"), "w") as f:
        f.write(nml)
    t0 = time.time()
    p = subprocess.run([os.path.join(d, "ilw.exe")], cwd=d,
                       capture_output=True, text=True)
    wall = time.time() - t0
    out = p.stdout
    rows = []
    if os.path.exists(os.path.join(d, "results.csv")):
        for line in open(os.path.join(d, "results.csv")).read().splitlines()[1:]:
            parts = line.split(",")
            if len(parts) < 11:
                continue
            rows.append(dict(level=int(parts[0]), nx=int(parts[1]), ny=int(parts[2]),
                             e1=float(parts[3]), einf=float(parts[4]),
                             r1=float(parts[5]), rinf=float(parts[6]),
                             nstep=int(parts[7]), cpu=float(parts[8]),
                             umin=float(parts[9]), mass=float(parts[10])))
    return dict(tag=tag, nml=nml, wall=wall, rows=rows, stdout=out,
                diverged=("divergence" in out))

def nml(**kw):
    keys = ["idim","iorder","nq","icase","ischeme","nlev","nx0","ny0","cfl","tend","tini","pmem","ipos","ibc","kconv","sconv","iops","nsrc"]
    s = "&para_nml\n"
    for k, v in kw.items():
        assert k in keys, k
        s += f" {k} = {v},\n"
    s += "/\n"
    return s

if __name__ == "__main__":
    print(run_case("smoke", nml(idim=1, iorder=6, nq=4, icase=1, ischeme=3,
                                nlev=2, nx0=40, cfl=0.4, tend=2.0)))
