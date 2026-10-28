"""lnD_ren_parts_common_vgeneral_phi_v2.py -- the six figures of the final
result ln D' and the numbers they show: every term of its numerical form,

    ln D'_{dim-reg, MSbar} = Sigma_dang + sum_{n=2}^{n_max} I_n(s2_max)
                             + sum_{n=2}^{n_max} T_n(s2_max) + Zcal
                             + ln D|_{k=1} + ln D|_{k=2},

split into its three parts, each against either cutoff (Lambda^2 of every
s^2-integral, or n_c of every sum over the partial waves), drawn the SAME way,
with the same legends and axis texts, one picture per part and cutoff:

    Numerical part               the computed part (blue) up to the cutoff,
                                 the cutoff (dotted) and the blue curve
                                 evaluated there (big dot); against n_c every
                                 curve of the three parts is points at the
                                 integers (a sum over whole waves has no
                                 value in between), against Lambda^2 the
                                 computed part is its data points joined by
                                 a line and the rest lines;
    Analytical part              what lies beyond the cutoffs -- the high-s^2
                                 tail, the high-n tail and the dangerous waves
                                 beyond s2_max -- against Lambda^2 each running
                                 towards its complete value (dotted); against
                                 n_c, as points at the integers, the sums of
                                 the high-s^2 tail and of the dangerous waves
                                 END at n_max and at n = 1, where they reach
                                 their complete values (big dots), only the
                                 high-n tail runs on (dotted);
    Numerical + analytical part  the computed part and, on top of it, the
                                 computed part plus everything the cutoffs
                                 miss (black), running onto ln D' (red); from
                                 a tenth of s2_max on (Lambda^2) and from
                                 n_max/2 on (n_c).

The plot codes plot_lnD_ren_<part>_vs_<s2|n> call contributions (the numbers,
computed ONCE, with all the gates below) and part_figure (the drawing).

TERM <-> PIPELINE (the seven-term assembly of assemble_lnD_vgeneral_phi_v2)
    Sigma_dang     (T1) Sigma_dang, the two dangerous waves n = 0, 1 (poles
                   subtracted, five fictitious modes added, their traces
                   dropped above s2_max)
    I_n(s2_max)    (T3) I_n_s2_max = -Int_0^{s2_max} ds^2 delta_geq3 of wave n
    T_n(s2_max)    (T5) T_n_s2_max, the high-s^2 tail of wave n beyond s2_max
    Zcal           (T6) Zcal, the high-n tail: the fit of I_n(s^2=inf),
                   a/nu^3 + b/nu^5 + c/nu^7 (nu = n + 1), summed over n > n_max
    ln D|_{k=1}    (T2)+(T4) A1_fin, the analytic mu-dependent tadpole term
    ln D|_{k=2}    (T7) -A2_fin/2, the analytic mu-dependent fish term

THE CURVES (contributions)
AGAINST Lambda^2, the upper limit of the s^2-integral (log axis):
    blue    Sigma_dang(Lambda^2) + sum_{n=2}^{n_max} I_n(Lambda^2)
            + ln D|_{k=1} + ln D|_{k=2}, up to s2_max where the data end.
            Sigma_dang(Lambda^2) is Sigma_dang integrated up to Lambda^2 (the
            two dangerous waves with their poles subtracted, the five
            fictitious modes at m_f^2 = m2fict): Delta_0(Lambda^2)
            + Delta_1(Lambda^2) - 5 ln(m_f^2 + Lambda^2) + ln|omega^2_{0,0}|,
            the Delta's by the production sector rule (running_sector_integral;
            below the third sector node the quadratic through the three
            lowest nodes, integrated exactly: low_s2_integral)
            carried onto the band nodes; I_n(Lambda^2) is I_n(s2_max) with
            s2_max -> Lambda^2 (the production quadrature); the two
            mu-dependent terms are constants;
    purple  what Sigma_dang still adds beyond s2_max, from s2_max on: its pole
            and fictitious terms there (the dangerous traces dropped as in
            (T1)), in closed form ln((omega^2_{0,0} + Lambda^2)/(omega^2_{0,0}
            + s2_max)) + 4 ln((omega^2_{1,0} + Lambda^2)/(omega^2_{1,0} + s2_max))
            - 5 ln((m_f^2 + Lambda^2)/(m_f^2 + s2_max));
    green   sum_n T_n, from s2_max on (the analytic prediction exists only
            above s2_max): sum_n T_n(s2_max) (1 - g_n(Lambda^2)), g_n the part
            of each wave's high-s^2 tail still above Lambda^2;
    orange  Zcal, its complete value alone (dotted): the high-n tail is a
            sum over the waves n > n_max of their WHOLE s^2-integrals, which
            no s^2-cutoff touches -- it does not run in Lambda^2.
AGAINST n_c, the cutoff of the sum over the partial waves:
    blue    the two mu-dependent terms plus the dangerous part up to s2_max --
            at n = 0 the negative-mode wave Delta_0 - ln(m_f^2 + s2_max)
            + ln|omega^2_{0,0}|, at n = 1 plus the zero-mode wave
            Delta_1 - 4 ln(m_f^2 + s2_max) (Sigma_dang wave by wave, each wave
            with its share 1 and 4 of the fictitious modes, integrated up to
            s2_max) -- then + sum_{m=2}^{n} I_m(s2_max) up to n_max: exactly the
            computed part, the same number at n_max as against Lambda^2 at
            s2_max;
    purple  what the two dangerous waves add beyond s2_max (their pole and
            fictitious terms, ln(m_f^2 + s2_max) - ln(omega^2_{0,0} + s2_max) at
            n = 0 and 4 [ln(m_f^2 + s2_max) - ln(omega^2_{1,0} + s2_max)] at
            n = 1), summed over n;
    green   sum_{m=2}^{n} T_m(s2_max), n = 2..n_max;
    orange  sum_{m=n_max+1}^{n} I_m(s^2=inf) from the fit, n > n_max: the
            analytic continuation, -> Zcal.
THE BLACK CURVE (numerical_plus_missed_curves): the blue curve plus
everything the cutoffs miss, each as far as it reaches at the running
cutoff --
    against Lambda^2: blue + Zcal for Lambda^2 <= s2_max, above it
                      blue(s2_max) + Zcal + sum_n T_n(Lambda^2)
                      + Sigma_dang beyond s2_max(Lambda^2);
    against n_c:      blue + Sigma_dang beyond s2_max (n = 0, 1)
                      + sum_{m=2}^{n} T_m(s2_max) for n <= n_max, above it
                      blue(n_max) + Sigma_dang beyond s2_max + sum_n T_n(s2_max)
                      + sum_{m=n_max+1}^{n} I_m(s^2=inf);
both run onto ln D'.

GATES (fail-loud, CHECK_TOL): the two dangerous waves add up to the stored
Sigma_dang; the running sector integrals end on the stored Delta's; the
running band sum ends on sum_n I_n(s2_max); I_n(s2_max) + T_n(s2_max) is the
stored I_n(s^2=inf); the stored T_n(s2_max) are today's completion of the
band tables; the end values against both cutoffs reproduce the stored ln D'
(and plot_common.certify_lnD has rebuilt the stored ln D' from its parts
before anything is computed).
"""
import os
import sys

import numpy as np
from scipy.integrate import simpson

import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = os.path.normpath(os.path.join(_HERE, '..', 'computation'))
for _p in (_HERE, _COMP):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pipeline_quadrature_vgeneral_phi_v2 import (I_n_from_delta, low_s2_strip,
                                                 high_s2_tail,
                                                 high_s2_remainder_ratios)
from plot_common_vgeneral_phi_v2 import (DOTTED_LINE_WIDTH, POINT_SIZE,
                                         band_tables_multikey, certify_lnD,
                                         load_run)

# one colour per contribution, the same in every figure
COLOR_BLUE, COLOR_T, COLOR_Z, COLOR_RESULT = ('tab:blue', 'tab:green',
                                              'tab:orange', 'tab:red')
COLOR_DANG = 'tab:purple'      # what Sigma_dang adds beyond s2_max
COLOR_TOTAL = 'black'          # the numerical part plus what the cutoffs miss

#: the curves against Lambda^2 run to LAMBDA2_TOP_FACTOR * s2_max, those
#: against n_c to N_PLOT_FACTOR * n_max
LAMBDA2_TOP_FACTOR = 1.0e4
N_PLOT_FACTOR = 5
#: tolerance of the bookkeeping gates (absolute, on numbers of order 1-100)
CHECK_TOL = 1e-9

# the figure symbols, written once
D_SUB = r'D_{\underset{\overline{MS}}{\text{dim-reg}}}'
LN_D_PRIME = r'\ln\left(D^{\prime}_{\underset{\overline{MS}}{\text{dim-reg}}}\right)'
ADDBACKS = (rf'\left.\ln\left({D_SUB}\right)\right|_{{k=1}} + '
            rf'\left.\ln\left({D_SUB}\right)\right|_{{k=2}}')
#: the dangerous waves split at s2_max, in the notation
#: Sigma_dang = -Int_0^inf ds^2 sigma_dang(s^2) - 5 ln(m_f^2) + ln|omega^2_{0,0}|:
#: the part up to s2_max (inside the computed part) and the part beyond it
DANG_UP_TO = (r'-\int_0^{s^2_{max}} ds^2\, \sigma_{dang}(s^2) - 5\ln(m_f^2)'
              r' + \ln(|\omega^2_{0,0}|)')
DANG_BEYOND = r'-\int_{s^2_{max}}^{\infty} ds^2\, \sigma_{dang}(s^2)'
#: the cutoff of every sum over the partial waves, upright (the index n of a
#: wave stays italic)
N_CUTOFF = r'\mathrm{n}_\mathrm{c}'
#: the label of the computed part, the same in every final-result figure (the
#: panel title says which cutoff runs); on two lines, so that it fits the
#: width of the picture
NUMERICAL_PART_LABEL = (
    rf'${DANG_UP_TO}$' '\n'
    rf'$+ \sum_{{n=2}}^{{n_{{max}}}} I_n(s^2_{{max}}) + {ADDBACKS}$')
#: what the cutoffs miss, in the same notation
MISSED = (rf'+ \sum_{{n=2}}^{{n_{{max}}}} T_n(s^2_{{max}}) {DANG_BEYOND}'
          r' + \mathcal{Z}')


# --------------------------------------------------------------------------- #
#  The numbers                                                                #
# --------------------------------------------------------------------------- #
def running_sector_integral(s2, delta_sub, pole):
    """-Int_0^{s2_j} delta_sub ds2 at every node s2_j (j >= 2) of a sector
    grid, by the production rule of that grid (pipeline_quadrature): a plain
    grid by I_n_from_delta on the truncated grid; a pole-clustered grid (the
    pole inside it) by the pole-split rule -- Simpson in ln|s2 - pole| per
    side, the trapezoid sliver across the pole, low_s2_strip below the first
    node.  At the last node it is the stored Delta_n (the caller
    gates).  Entries 0, 1 are NaN (the rules need three nodes)."""
    s2, delta_sub = np.asarray(s2, float), np.asarray(delta_sub, float)
    out = np.full(len(s2), np.nan)
    if not (s2[0] < pole < s2[-1]):
        for j in range(2, len(s2)):
            out[j] = float(I_n_from_delta(s2[:j + 1], delta_sub[:j + 1]))
        return out
    i_above = int(np.searchsorted(s2, pole))       # first node above the pole
    s2_below, delta_below = s2[:i_above], delta_sub[:i_above]
    strip = low_s2_strip(s2, delta_sub)
    for j in range(2, len(s2)):
        if j < i_above:
            integral = simpson(delta_sub[:j + 1] * (s2[:j + 1] - pole),
                               x=np.log(pole - s2[:j + 1]))
        else:
            integral = (simpson(delta_below * (s2_below - pole),
                                x=np.log(pole - s2_below))
                        + 0.5 * (delta_below[-1] + delta_sub[i_above])
                        * (s2[i_above] - s2_below[-1]))
            if j > i_above:
                s2_above, delta_above = s2[i_above:j + 1], delta_sub[i_above:j + 1]
                integral += simpson(delta_above * (s2_above - pole),
                                    x=np.log(s2_above - pole))
        out[j] = -(integral + strip)
    return out


def low_s2_integral(s2, delta, x):
    """Int_0^x delta ds2 for every x below the third node s2[2] of a sector
    grid, where the sector rules (they need three nodes) do not reach: the
    exact integral of the quadratic through the three lowest nodes -- the
    rule pipeline_quadrature.low_s2_strip uses below the first node (there
    x = s2[0]), run up to x."""
    x0, x1, x2 = float(s2[0]), float(s2[1]), float(s2[2])
    x = np.asarray(x, float)

    def moment(a, b):                  # Int_0^x (s2-a)(s2-b) ds2
        return x**3 / 3.0 - 0.5 * (a + b) * x**2 + a * b * x

    w0 = moment(x1, x2) / ((x0 - x1) * (x0 - x2))    # Lagrange weights of the
    w1 = moment(x0, x2) / ((x1 - x0) * (x1 - x2))    # three lowest nodes
    w2 = moment(x0, x1) / ((x2 - x0) * (x2 - x1))
    return w0 * delta[0] + w1 * delta[1] + w2 * delta[2]


def gate(name, value, target):
    """Fail loud when a bookkeeping identity does not hold."""
    if not (np.isfinite(value) and abs(value - target) <= CHECK_TOL * max(1.0, abs(target))):
        raise RuntimeError(f'[ABORT] {name}: {value!r} != {target!r} '
                           f'(dev {value - target:+.3e}).')


def contributions(args):
    """Every term of ln D', running in Lambda^2 and in n, with all the
    gates of the module docstring -- the ONE owner of these numbers.  Returns
    a dict of the curves and the end values."""
    # the run and its stored ln D', rebuilt from its parts and checked first
    run = load_run(args.data_dir, args.tag, args.bounce_npz,
                   need=('lnD_ren', 'band_integrals'))
    lnD_ren, parts = certify_lnD(run)
    lnD_ren_npz = run['lnD_ren']
    band_n = np.asarray(lnD_ren_npz['band_n'], int)              # band waves 2..n_max
    I_inf = np.asarray(lnD_ren_npz['I_n_inf'], float)     # completed per-wave I_n(inf)
    if np.any(~np.isfinite(I_inf)):
        raise RuntimeError('[ABORT] non-finite I_n(inf) in the lnD_ren npz '
                           '(fail-loud, no nansum).')
    n_max, s2_max = int(lnD_ren_npz['n_max']), float(lnD_ren_npz['s2_max'])
    # the production fit a/nu^3 + b/nu^5 + c/nu^7 of I_n(inf)
    a, b, c = (float(lnD_ren_npz['fit_a']), float(lnD_ren_npz['fit_b']),
               float(lnD_ren_npz['fit_c']))
    Zcal = float(lnD_ren_npz['Zcal'])
    background = run['bg']
    nu = band_n + 1.0

    I_n_s2_max = np.asarray(lnD_ren_npz['I_n_s2_max'], float)
    T_n_s2_max = np.asarray(lnD_ren_npz['T_n_s2_max'], float)
    gate('I_n(s2_max) + T_n(s2_max) = I_n(inf)',
         float(np.max(np.abs(I_n_s2_max + T_n_s2_max - I_inf))), 0.0)
    T_sum = float(np.sum(T_n_s2_max))
    ln_D_k1, ln_D_k2 = float(parts['A1_fin']), float(parts['minus_half_A2_fin'])
    omega2_00 = float(lnD_ren_npz['omega2_00'])
    omega2_10 = float(lnD_ren_npz['omega2_10'])
    deg1, m2_f = float(lnD_ren_npz['deg1']), float(lnD_ren_npz['m2fict'])
    Sigma_dang = float(parts['Sigma_dang'])

    # ---- the dangerous waves: Sigma_dang wave by wave -----------------------
    # each wave with its share of the fictitious modes (1 and deg1), split at
    # s2_max: the part up to s2_max belongs to the computed part, the part
    # beyond it (the pole and fictitious terms there, the dangerous traces
    # dropped as in (T1)) is what the s2_max cutoff misses
    Delta = (float(lnD_ren_npz['Delta_0']),
             float(lnD_ren_npz['Delta_1']))
    Sigma_dang_n0 = Delta[0] - np.log(m2_f + s2_max) + np.log(abs(omega2_00))
    Sigma_dang_n1 = Delta[1] - deg1 * np.log(m2_f + s2_max)
    dang_beyond_n0 = np.log(m2_f + s2_max) - np.log(omega2_00 + s2_max)
    dang_beyond_n1 = deg1 * (np.log(m2_f + s2_max) - np.log(omega2_10 + s2_max))
    gate('the two dangerous waves, up to and beyond s2_max, = Sigma_dang',
         Sigma_dang_n0 + Sigma_dang_n1 + dang_beyond_n0 + dang_beyond_n1, Sigma_dang)

    # ---- against Lambda^2 ----------------------------------------------------
    # the band waves n = 2..n_max, running in Lambda^2 (production quadrature)
    tab = band_tables_multikey(args.data_dir, args.tag, args.bounce_npz,
                               n_max=n_max, s2_cap=s2_max)
    s2, delta_geq3 = np.asarray(tab['s2'], float), np.asarray(tab['delta_geq3'], float)
    Lambda2 = s2[2:]
    I_band_running = np.array([float(np.sum(I_n_from_delta(s2[:j + 1], delta_geq3[:, :j + 1])))
                               for j in range(2, len(s2))])
    gate('running band sum at s2_max = sum_n I_n(s2_max)',
         I_band_running[-1], float(np.sum(I_n_s2_max)))
    # the two dangerous waves, running in Lambda^2 on their own sector grids,
    # carried onto the band nodes (both grids end on s2_max); below the third
    # sector node, where the running integral starts, the quadratic through
    # the three lowest nodes is integrated from 0 to Lambda^2 (np.interp
    # would hold the value at the third node there)
    Delta_running = []
    for n, omega2 in ((0, omega2_00), (1, omega2_10)):
        sector_npz = run.get(f'sector_n{n}')
        if sector_npz is None:
            raise SystemExit(f'[ABORT] sector_n{n} npz not found -- the dangerous '
                             f'wave n = {n} is drawn from it.')
        s2_sector = np.asarray(sector_npz['s2_sector'], float)
        delta_sub = np.asarray(sector_npz['delta_sub'], float)
        Delta_sector = running_sector_integral(s2_sector, delta_sub, -omega2)
        gate(f'running sector integral n = {n} at s2_max = Delta_{n}',
             Delta_sector[-1], Delta[n])
        Delta_n = np.interp(np.log(Lambda2), np.log(s2_sector[2:]),
                            Delta_sector[2:])
        low = Lambda2 < s2_sector[2]
        Delta_n[low] = -low_s2_integral(s2_sector, delta_sub, Lambda2[low])
        Delta_running.append(Delta_n)
    Sigma_dang_running = (Delta_running[0] + Delta_running[1]
                          - (1.0 + deg1) * np.log(m2_f + Lambda2) + np.log(abs(omega2_00)))
    numerical_part_vs_s2 = (Sigma_dang_running + I_band_running
                            + ln_D_k1 + ln_D_k2)
    Sigma_dang_beyond = Sigma_dang - Sigma_dang_running[-1]
    # the high-s^2 tail, from s2_max on; the stored T_n(s2_max) must be the
    # completion the band tables' endpoint column gives today
    T_n_today, _ = high_s2_tail(nu, s2_max, background['m_min2'], background['R'],
                                      background['tr_U_int3'], delta_geq3[:, -1], quiet=True)
    gate('stored T_n(s2_max) = completion of the band tables',
         float(np.max(np.abs(T_n_s2_max - T_n_today))), 0.0)
    Lambda2_T = np.geomspace(s2_max, LAMBDA2_TOP_FACTOR * s2_max, 70)
    g_T, _ = high_s2_remainder_ratios(nu, s2_max, Lambda2_T, background['m_min2'],
                                      background['R'], background['tr_U_int3'])
    T_running = (1.0 - g_T) @ T_n_s2_max
    # what Sigma_dang adds beyond s2_max, running from s2_max on: the pole and
    # fictitious terms of Sigma_dang integrated from s2_max to Lambda^2 (the
    # dangerous traces dropped there, as in (T1)); its Lambda^2 -> inf
    # limit is Sigma_dang_beyond
    def dang_beyond(L2):
        return (np.log((omega2_00 + L2) / (omega2_00 + s2_max))
                + deg1 * np.log((omega2_10 + L2) / (omega2_10 + s2_max))
                - (1.0 + deg1) * np.log((m2_f + L2) / (m2_f + s2_max)))
    dang_running = dang_beyond(Lambda2_T)
    gate('Sigma_dang beyond s2_max = its closed-form limit',
         (1.0 + deg1) * np.log(m2_f + s2_max) - np.log(omega2_00 + s2_max)
         - deg1 * np.log(omega2_10 + s2_max), Sigma_dang_beyond)
    # the high-n tail does not run in Lambda^2: every wave n > n_max enters
    # with its whole s^2-integral (the fit), so against Lambda^2 it is the
    # constant Zcal; the black curve's grid is the band nodes followed by the
    # tail grid
    Lambda2_Z = np.concatenate([Lambda2, Lambda2_T[1:]])
    sum_s2 = numerical_part_vs_s2[-1] + Sigma_dang_beyond + T_sum + Zcal
    gate('against Lambda^2: blue + Sigma_dang beyond + sum T_n + Zcal = ln D\'', sum_s2, lnD_ren)
    print(f'[lnD_ren_parts] against Lambda^2: ln D\' = '
          f'{numerical_part_vs_s2[-1]:+.6f} (blue at '
          f's2_max) {Sigma_dang_beyond:+.6f} (Sigma_dang beyond s2_max) {T_sum:+.6f} '
          f'(sum T_n) {Zcal:+.6f} (Zcal) = {sum_s2:+.6f}')

    # ---- against n -----------------------------------------------------------
    n_numerical_part = np.concatenate([[0, 1], band_n])
    numerical_part_vs_n = ln_D_k1 + ln_D_k2 + np.concatenate(
        [[Sigma_dang_n0, Sigma_dang_n0 + Sigma_dang_n1],
         Sigma_dang_n0 + Sigma_dang_n1 + np.cumsum(I_n_s2_max)])
    gate('the computed part at n_max = the computed part at s2_max',
         numerical_part_vs_n[-1], numerical_part_vs_s2[-1])
    dang_beyond_cumulative = np.array([dang_beyond_n0, dang_beyond_n0 + dang_beyond_n1])
    T_cumulative = np.cumsum(T_n_s2_max)
    n_tail = np.arange(n_max + 1, N_PLOT_FACTOR * n_max + 1)
    Z_partial = np.cumsum(a / (n_tail + 1.0)**3 + b / (n_tail + 1.0)**5
                          + c / (n_tail + 1.0)**7)
    sum_n = (numerical_part_vs_n[-1] + dang_beyond_cumulative[-1] + T_sum
             + Zcal)
    gate('against n: blue + Sigma_dang beyond + sum T_n + Zcal = ln D\'', sum_n, lnD_ren)
    print(f'[lnD_ren_parts] against n: ln D\' = '
          f'{numerical_part_vs_n[-1]:+.6f} (blue at '
          f'n_max) {dang_beyond_cumulative[-1]:+.6f} (Sigma_dang beyond s2_max) '
          f'{T_sum:+.6f} (sum T_n) {Zcal:+.6f} (Zcal) = {sum_n:+.6f}')
    return dict(lnD_ren=lnD_ren, s2_max=s2_max, n_max=n_max, Lambda2=Lambda2,
                numerical_part_vs_s2=numerical_part_vs_s2,
                Lambda2_T=Lambda2_T, T_running=T_running,
                dang_running=dang_running, Lambda2_Z=Lambda2_Z,
                Sigma_dang_beyond=Sigma_dang_beyond, T_sum=T_sum, Zcal=Zcal,
                sum_s2=sum_s2, n_numerical_part=n_numerical_part,
                numerical_part_vs_n=numerical_part_vs_n, band_n=band_n,
                T_cumulative=T_cumulative, n_tail=n_tail, Z_partial=Z_partial,
                sum_n=sum_n, dang_beyond_cumulative=dang_beyond_cumulative)


def numerical_plus_missed_curves(contributions):
    """The blue curve plus everything the cutoffs miss: against Lambda^2 on
    the band nodes followed by the tail grid (Lambda2_Z; the high-n tail
    Zcal is a constant there), against n on n = 0..N_PLOT_FACTOR n_max.
    Returns numerical_plus_missed_vs_s2, n_all, numerical_plus_missed_vs_n."""
    numerical_plus_missed_vs_s2 = np.concatenate([
        contributions['numerical_part_vs_s2'] + contributions['Zcal'],
        contributions['numerical_part_vs_s2'][-1] + contributions['Zcal']
        + contributions['T_running'][1:] + contributions['dang_running'][1:]])
    # at n = 0 and from n = 1 on
    dang_beyond = contributions['dang_beyond_cumulative']
    n_all = np.concatenate([contributions['n_numerical_part'],
                            contributions['n_tail']]).astype(float)
    numerical_plus_missed_vs_n = np.concatenate([
        contributions['numerical_part_vs_n'][:2] + dang_beyond,
        contributions['numerical_part_vs_n'][2:] + dang_beyond[-1]
        + contributions['T_cumulative'],
        contributions['numerical_part_vs_n'][-1] + dang_beyond[-1]
        + contributions['T_sum'] + contributions['Z_partial']])
    return numerical_plus_missed_vs_s2, n_all, numerical_plus_missed_vs_n


def lambda2_window(contributions):
    """The Lambda^2 axis of every final-result figure: the band nodes to the
    end of the tail grid."""
    return (contributions['Lambda2'][0] / 1.5,
            contributions['Lambda2_T'][-1] * 1.5)


def n_window(contributions):
    """The n axis of every final-result figure (the points at n = 0, 1 off
    the spine)."""
    return -0.01 * contributions['n_tail'][-1], contributions['n_tail'][-1]


def analytic_window(contributions):
    """The y axis of the analytical-part figures: the three parts and zero,
    whatever their signs."""
    values = (contributions['T_sum'], contributions['Zcal'],
              contributions['Sigma_dang_beyond'], 0.0)
    lowest, highest = min(values), max(values)
    return (min(1.18 * lowest, -0.12 * highest),
            max(-0.12 * lowest, 1.5 * highest))


# --------------------------------------------------------------------------- #
#  The drawing                                                                #
# --------------------------------------------------------------------------- #
#: how a curve is drawn (sizes: plot_common, the one set): against Lambda^2,
#: a continuous cutoff, the computed part is its data points joined by a line
#: and every other curve a LINE; against n_c, which counts partial waves, a
#: value exists only at the integers n_c = 0, 1, 2, ... -- every curve is
#: POINTS alone
LINE = dict(fmt='-')
POINTS = dict(fmt='o', ms=POINT_SIZE, mew=0)
POINTS_JOINED = dict(fmt='o-', ms=POINT_SIZE, mew=0)
#: a value reached at a cutoff: a big dot there
BIG_POINT_SIZE = 8.0
#: the complete value of an analytical contribution: approached for the
#: cutoff -> infinity (a dotted line) or reached AT a cutoff (a big dot there)
APPROACHED = dict(at=None)


def reached_at(x):
    return dict(at=x)


def plot_styled(ax, x, y, style, **kw):
    """ax.plot with one of the styles above."""
    style = dict(style)
    ax.plot(x, y, style.pop('fmt'), **style, **kw)


def cutoff_view(contributions, cutoff):
    """Everything the three panels need against one cutoff, 's2' or 'n'."""
    (numerical_plus_missed_vs_s2, n_all,
     numerical_plus_missed_vs_n) = numerical_plus_missed_curves(contributions)
    if cutoff == 's2':
        s2_max = contributions['s2_max']
        return dict(cutoff=r'$\Lambda^2$', xlabel=r'$\Lambda^2$',
                    log_x=True, window=lambda2_window(contributions),
                    total_window=(s2_max / 10.0,
                                  contributions['Lambda2_Z'][-1]),
                    x_numerical_part=contributions['Lambda2'],
                    numerical_part=contributions['numerical_part_vs_s2'],
                    x_numerical_plus_missed=contributions['Lambda2_Z'],
                    numerical_plus_missed=numerical_plus_missed_vs_s2,
                    numerical_style=POINTS_JOINED, missed_style=LINE,
                    x_cut=s2_max,
                    cut=r'$s^2_{max}$', cut_value=rf'$s^2_{{max}} = {s2_max:g}$',
                    T=(contributions['Lambda2_T'], contributions['T_running'],
                       LINE),
                    Z=None,            # a constant against Lambda^2: dotted alone
                    dang=(contributions['Lambda2_T'],
                          contributions['dang_running'], LINE),
                    T_complete=APPROACHED, Z_complete=APPROACHED,
                    dang_complete=APPROACHED)
    n_max = contributions['n_max']
    return dict(cutoff=rf'${N_CUTOFF}$', xlabel=rf'${N_CUTOFF}$',
                log_x=False, window=n_window(contributions),
                total_window=(0.5 * n_max, n_all[-1]),
                x_numerical_part=contributions['n_numerical_part'].astype(float),
                numerical_part=contributions['numerical_part_vs_n'],
                x_numerical_plus_missed=n_all,
                numerical_plus_missed=numerical_plus_missed_vs_n,
                numerical_style=POINTS, missed_style=POINTS, x_cut=n_max,
                cut=r'$n_{max}$', cut_value=rf'$n_{{max}} = {n_max}$',
                T=(contributions['band_n'], contributions['T_cumulative'],
                   POINTS),
                Z=(contributions['n_tail'], contributions['Z_partial'], POINTS),
                dang=(np.array([0.0, 1.0]),
                      contributions['dang_beyond_cumulative'], POINTS),
                T_complete=reached_at(n_max), Z_complete=APPROACHED,
                dang_complete=reached_at(1.0))


def x_axis(ax, v, window):
    """The cutoff axis of a panel."""
    if v['log_x']:
        ax.set_xscale('log')
    ax.set_xlim(*window)
    ax.grid(True, which='both')     # the one grid style of plot_common


def draw_numerical(ax, contributions, v):
    """Numerical part: the computed part (blue) up to the cutoff, the cutoff
    (dotted) and the blue curve evaluated there (big dot)."""
    plot_styled(ax, v['x_numerical_part'], v['numerical_part'],
                v['numerical_style'], color=COLOR_BLUE, zorder=4,
                label=NUMERICAL_PART_LABEL)
    ax.axvline(v['x_cut'], color='0.35', ls=':', lw=DOTTED_LINE_WIDTH, zorder=3,
               label=v['cut_value'])
    ax.plot(v['x_numerical_part'][-1], v['numerical_part'][-1], 'o',
            color=COLOR_BLUE, ms=BIG_POINT_SIZE, zorder=5,
            label=(rf'blue curve evaluated at {v["cut"]} '
                   rf'$= {v["numerical_part"][-1]:.6f}$'))
    x_axis(ax, v, v['window'])
    # 8 units above and below the curve, a tenth of its range if that is less
    numerical_part = v['numerical_part']
    pad = min(8.0, 0.1 * (numerical_part.max() - numerical_part.min()))
    ax.set_ylim(numerical_part.min() - pad, numerical_part.max() + pad)
    ax.set_ylabel(rf'Numerical contribution to ${LN_D_PRIME}$')


def draw_analytical(ax, contributions, v):
    """Analytical part: what lies beyond the cutoffs, each part followed by
    its complete value -- approached for the cutoff -> infinity (dotted) or
    reached at a cutoff (big dot) -- in the order dangerous waves, high-s^2
    tail, high-n tail, then the cutoff (dotted) with its value, in one legend
    column; a part that does not run in this cutoff (Zcal against Lambda^2)
    is its dotted constant alone, labelled by its name."""
    parts = (
        (v['dang'], contributions['Sigma_dang_beyond'], COLOR_DANG,
         rf'${DANG_BEYOND}$', DANG_BEYOND, v['dang_complete']),
        (v['T'], contributions['T_sum'], COLOR_T,
         r'$\sum_{n=2}^{n_{max}} T_n$: the high-$s^2$ tail',
         r'\sum_{n=2}^{n_{max}} T_n(s^2_{max})', v['T_complete']),
        (v['Z'], contributions['Zcal'], COLOR_Z,
         r'$\mathcal{Z}$: the high-$n$ tail', r'\mathcal{Z}', v['Z_complete']))
    for curve, value, color, label, name, complete in parts:
        if curve is not None:
            x, y, style = curve
            plot_styled(ax, x, y, style, color=color, zorder=4, label=label)
            value_label = rf'Complete value of: ${name} = {value:+.4f}$'
        else:                       # no running curve: the constant is the part
            value_label = label
        if complete['at'] is None:
            ax.axhline(value, color=color, ls=':', lw=DOTTED_LINE_WIDTH, zorder=2,
                       label=value_label)
        else:
            ax.plot(complete['at'], value, 'o', color=color, ms=BIG_POINT_SIZE,
                    zorder=5, label=value_label)
    # the cutoff (dotted) with its value, last
    ax.axvline(v['x_cut'], color='0.35', ls=':', lw=DOTTED_LINE_WIDTH, zorder=3,
               label=v['cut_value'])
    ax.axhline(0.0, color='0.3')
    x_axis(ax, v, v['window'])
    ax.set_ylim(*analytic_window(contributions))
    ax.set_ylabel(rf'Analytical contribution to ${LN_D_PRIME}$')


def draw_total(ax, contributions, v):
    """Numerical + analytical part: the computed part (blue) and, on top of
    it, the computed part plus everything the cutoffs miss (black), running
    onto ln D' (red)."""
    lnD_ren = contributions['lnD_ren']
    plot_styled(ax, v['x_numerical_part'], v['numerical_part'],
                v['numerical_style'], color=COLOR_BLUE, zorder=4,
                label=NUMERICAL_PART_LABEL)
    plot_styled(ax, v['x_numerical_plus_missed'], v['numerical_plus_missed'],
                v['missed_style'], color=COLOR_TOTAL, zorder=5,
                label=rf'Numerical part ${MISSED}$')
    ax.axvline(v['x_cut'], color='0.35', ls=':', lw=DOTTED_LINE_WIDTH, zorder=3,
               label=v['cut_value'])
    ax.axhline(lnD_ren, color=COLOR_RESULT, ls='--', zorder=3,
               label=rf'${LN_D_PRIME} = {lnD_ren:.6f}$')
    lo = v['total_window'][0]
    x_axis(ax, v, v['total_window'])
    # both curves inside the window and ln D', whichever side the curves
    # approach it from
    numerical_part = v['numerical_part'][v['x_numerical_part'] >= lo]
    numerical_plus_missed = v['numerical_plus_missed'][
        v['x_numerical_plus_missed'] >= lo]
    bottom = min(numerical_part.min(), numerical_plus_missed.min(), lnD_ren)
    top = max(numerical_part.max(), numerical_plus_missed.max(), lnD_ren)
    span = top - bottom
    ax.set_ylim(bottom - 0.15 * span, top + 0.08 * span)
    ax.set_ylabel(rf'Total contribution to ${LN_D_PRIME}$')


#: the three parts: (title, drawer, legend columns) -- the analytical part
#: lists each curve followed by its limit in one column, the summed part the
#: curves with the lines in two
PARTS = {'numerical': ('Numerical part', draw_numerical, 1),
         'analytical': ('Analytical part', draw_analytical, 1),
         'numerical_plus_analytical': ('Numerical + analytical part', draw_total, 2)}


def part_figure(contributions, cutoff, part):
    """One part of the final result against one cutoff ('s2' or 'n'): a
    single panel titled '<part> vs. Cutoff <Lambda^2 | n_c>', its legend below.
    Returns the figure."""
    v = cutoff_view(contributions, cutoff)
    title, draw, ncol = PARTS[part]
    # a picture <= 8 in wide (sizes: plot_common, the one set)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    draw(ax, contributions, v)
    ax.set_title(f'{title} vs. Cutoff {v["cutoff"]}')
    ax.set_xlabel(v['xlabel'])
    # below the panel, clear of its x label (the pad in font sizes)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), borderaxespad=3.3,
              ncol=ncol, handlelength=1.8, columnspacing=1.2)
    return fig
