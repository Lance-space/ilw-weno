# ILW-WENO

One-step improved Lax-Wendroff (ILW) time discretization for high-order WENO
discretizations of nonlinear degenerate parabolic equations.

The code accompanies the manuscript

> L. Li, J. Zhu, *One-step improved Lax-Wendroff time discretization for
> high-order WENO schemes applied to nonlinear degenerate parabolic equations*,
> submitted.

and reproduces all the tables and figures of its Section 6.

## What is implemented

The spatial operator is a mapped WENO reconstruction of the flux, sixth order
(`WENO6`) or eighth order (`WENO8`), written as a difference of half-point fluxes;
for the convection-diffusion problems a fifth-order WENO-JS discretization of the
first-order term is added.

Four time integrators are available and share the same spatial operator:

| `ischeme` | method | evaluations of `L` per step | evaluations of `C` per step |
|---|---|---|---|
| 3  | SSP-RK3        | 3 | 0 |
| 13 | ILW3           | 1 | 3 |
| 4  | classical RK4  | 4 | 0 |
| 14 | ILW4           | 1 | 4 |

`L` is the WENO operator and `C` the cheap centred difference of the same order.
The ILW step is

    u^{n+1} = u^n + dt [ L(u^n) + sum_{j>=2} b_j C(V_j) - (sum_{j>=2} b_j) C(u^n) ]
    V_j     = u^n + dt [ sum_{l<j} a_jl C(V_l) + sigma_j d ],   sigma_j = sum_{l<j} a_jl

with the operator difference `d = L(u^n) - C(u^n)` frozen at the start of the step.

## Two ways to run it

### 1. One file per example (recommended)

`examples/ex6_k/<scheme>/` contains a self-contained Fortran program for example
`k` with one of the four schemes; the six modules are inlined and the parameters
are assigned in the source, so nothing else is needed.

    f90 /optimize:4 /real_size:64 ilw_ex6_1_ILW3_WENO6.f90      (Compaq/Intel)
    gfortran -O3 -fdefault-real-8 ilw_ex6_1_ILW3_WENO6.f90      (gfortran)

The program prints, for each grid level, the `L1` and `Linf` errors, the
convergence rate, the number of steps and the CPU time, and writes `results.csv`
and `profile.plt`.

### 2. One executable, parameters in a namelist

    cd src
    f90 /optimize:4 /real_size:64 comm.f90 para.f90 funcs.f90 ops.f90 timeint.f90 main.f90
    ./ilw.exe          (reads input.nml from the working directory)

The namelist fields are `idim, iorder, nq, icase, ischeme, nlev, nx0, ny0, cfl,
tend, tini, pmem, ipos, ibc, kconv, sconv, iops, nsrc, ieq, epsbl, iblflux, ubcl,
ubcr`.  The `input.nml` files used for the paper are those of
`examples/ex6_k/<scheme>/` in the alternative distribution of the paper; the
standalone programs contain the same values.

## The nine examples

| # | equation | domain | data | `T` | grid |
|---|---|---|---|---|---|
| 1 | `u_t = u_xx` | `[-pi,pi]` periodic | `u0 = exp(-t0) sin x`, `t0 = 0.1` | 2 | `N = 20,40,80,160` |
| 2 | `u_t = (u^m)_xx` | `[-6,6]`, `u = 0` at both ends | Barenblatt at `t = 1`, `m = 2,3,5,8` | 2 | `N = 200, 400` |
| 3 | `u_t = (u^6)_xx` | `[-6,6]`, `u = 0` | two boxes, heights 1 and 2 | 0.8 | `N = 160` |
| 4 | Buckley-Leverett, `eps = 0.01` | `[0,1]`, inflow `u(0,t) = 1` | `u0 = 1-x` on `[0,1/3]` | 0.2 | `N = 400` |
| 5 | strongly degenerate | `[-2,2]`, `u = 0` | two pulses of height `+-1` | 0.7 | `N = 400` |
| 6 | `u_t = u_xx + u_yy` | `[-pi,pi]^2` periodic | `u0 = exp(-2 t0) sin(x+y)` | 1 | `N^2 = 10^2 ... 80^2` |
| 7 | `u_t = (u^2)_xx + (u^2)_yy` | `[-10,10]^2` periodic | two Gaussian caps | 0.5 | `80^2` |
| 8 | 2D Buckley-Leverett, `eps = 0.01` | `[-1.5,1.5]^2` periodic | unit disc | 0.5 | `120^2` |
| 9 | 2D strongly degenerate | `[-1.5,1.5]^2` periodic | two discs of `+-1` | 0.5 | `120^2` |

`CFL = 0.4` throughout, except example 2 with `m = 5`, which uses `CFL = 0.2`.
`nq = 4` is paired with `WENO6` and `nq = 6` with `WENO8`.

## Expected output

`examples/ex6_1/ILW3-WENO6` prints

    level  1  nx  20   L1 4.7686E-07
    level  2  nx  40   L1 7.8337E-09   rate L1 5.93
    level  3  nx  80   L1 1.2401E-10   rate L1 5.98
    level  4  nx 160   L1 1.9502E-12   rate L1 5.99   3243 steps

which is the accuracy table of example 1 of the paper.

## Results

`results/nine_results.json` holds, for each of the nine examples and each of the
four schemes, the number of steps, the errors and the computing time reported in
Table 1 of the paper; `results/diff_results.json` holds the pointwise differences
between the ILW and RK solutions quoted in the text.  The scripts in `scripts/`
read the `profile.plt` files written by the runs and regenerate those numbers.

## Acknowledgement of support

L. Li is supported by the Natural Science Foundation of Henan Province
(252300423517) and NSFC (12501548).  J. Zhu is supported by NSFC (12472292) and
the Science Challenge Project (TZ2025007).

## License

MIT, see `LICENSE`.

## Contact

Liang Li, liliangnuaa@163.com; Jun Zhu, zhujun@nuaa.edu.cn.
