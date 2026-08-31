"""
Drop-in replacements for getDoublingTime and getLoss.

Fixes:
  * doubling time from the mu_max exponential-phase window (log-linear fit),
    not from global max/min OD over the whole run
  * survival curve read as a probability, not a percentage
  * estimating equation in the substituted (exposure-age) form, which is
    better conditioned: the exp((beta-L)*T) prefactor cancels analytically
    instead of being computed and divided away
  * loss rate converted to per-division as q = L / beta, not L * d
  * explicit guards and a self-widening bracket, so a failure names its cause
    instead of surfacing as "f(a) and f(b) must have different signs"

Paste these two methods over the existing ones in findPlasmidLossRate.
Requires: import math, numpy as np, pandas as pd, scipy as sp
"""

import math
import numpy as np
import pandas as pd
import scipy as sp


def getDoublingTime(self, col="Top10 / pGFPck Δ6", min_points=2):
    """Doubling time from the exponential-phase window.

    Exponential phase is the contiguous window of time points maximising the
    specific growth rate mu, estimated by regressing ln(OD600) on time.
    Returns d = ln(2) / mu_max in hours.
    """
    df = self.avg_odpcb_df
    if "Hours" not in df.columns:
        raise KeyError("expected a 'Hours' column in the OD data; found %r" % list(df.columns))
    if col not in df.columns:
        raise KeyError("strain column %r not in OD data; found %r" % (col, list(df.columns)))

    t = pd.to_numeric(df["Hours"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(t) & np.isfinite(y) & (y > 0)
    t, y = t[keep], y[keep]
    if t.size < min_points:
        raise ValueError("need at least %d usable OD points, got %d" % (min_points, t.size))

    ly = np.log(y)
    best = None
    for i in range(t.size):
        for j in range(i + min_points, t.size + 1):
            slope = np.polyfit(t[i:j], ly[i:j], 1)[0]
            if best is None or slope > best[0]:
                best = (slope, t[i], t[j - 1])

    mu, t0, t1 = best
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError("no window with positive growth rate; check the OD column")
    self.growth_window = (t0, t1)
    self.mu_max = mu
    return math.log(2) / mu


def getLoss(self, row: pd.Series, bracket=(1e-12, 50.0), max_widen=6) -> dict:
    """Per-division plasmid loss probability for one replicate, by treatment.

    Solves   R = 1 / ( L * INT_0^T exp(-(beta - L) s) p(s) ds )   for L,
    then returns q = L / beta.
    """
    beta = math.log(2) / self.d_time
    rates = {}

    for treatment in row.index:
        treat = treatment.split(".")[0]
        prop_dark = row[treatment]
        if pd.isna(prop_dark):
            continue

        prop = float(prop_dark)
        if not (0.0 < prop < 1.0):
            raise ValueError(
                "{}: observed plasmid-free proportion is {:g}. The model needs "
                "0 < p < 1. p = 0 makes R infinite and p = 1 makes R zero; in "
                "either case the root-finder has no sign change to bracket and "
                "scipy reports \"f(a) and f(b) must have different signs\". "
                "A replicate with zero plasmid-free colonies has no point "
                "estimate: report it as 0/n with a one-sided upper bound "
                "instead of passing it through the model.".format(treatment, prop)
            )
        R = (1.0 - prop) / prop

        hours = pd.to_numeric(self.avg_via_df["Hours"], errors="coerce")
        surv = pd.to_numeric(self.avg_via_df[treat], errors="coerce") / 100.0  # percent -> probability
        ok = surv.notna() & hours.notna()
        s_obs = hours[ok].to_numpy(dtype=float)
        p_obs = surv[ok].to_numpy(dtype=float)
        if s_obs.size < 2:
            raise ValueError("%s: fewer than two usable survival time points" % treat)
        if p_obs.max() > 1.0 + 1e-9:
            raise ValueError(
                "%s: survival values exceed 1 after dividing by 100 (max %.4g). "
                "The kill-curve file is probably already on a 0-1 scale, so the "
                "/100 is being applied twice." % (treat, p_obs.max())
            )

        T = float(s_obs.max())
        p = lambda s: float(np.interp(s, s_obs, p_obs))

        def f(L):
            integral = sp.integrate.quad(
                lambda s: math.exp(-(beta - L) * s) * p(s), 0.0, T, limit=200
            )[0]
            if integral <= 0.0:
                return float("inf")
            return 1.0 / (L * integral) - R

        lo, hi = bracket
        f_lo, f_hi = f(lo), f(hi)
        widened = 0
        while f_lo * f_hi > 0 and widened < max_widen:
            # f is monotone decreasing in L, so a same-sign pair means the root
            # is outside the bracket: push whichever end is on the wrong side.
            if f_hi > 0:
                hi *= 10.0
                f_hi = f(hi)
            else:
                lo /= 100.0
                f_lo = f(lo)
            widened += 1

        if not (np.isfinite(f_lo) and np.isfinite(f_hi)):
            raise ValueError(
                "%s: objective is not finite at the bracket ends "
                "(f(%g) = %r, f(%g) = %r). beta = %.4f /h, T = %.4f h, R = %.4g. "
                "A very large beta overflows exp((beta - L) * T); check that "
                "getDoublingTime returned a sensible doubling time (got %.5g h)."
                % (treat, lo, f_lo, hi, f_hi, beta, T, R, self.d_time)
            )
        if f_lo * f_hi > 0:
            raise ValueError(
                "%s: no sign change on [%g, %g] after widening "
                "(f(lo) = %.6g, f(hi) = %.6g). beta = %.4f /h, T = %.4f h, "
                "observed proportion = %.6g, R = %.6g."
                % (treat, lo, hi, f_lo, f_hi, beta, T, prop, R)
            )

        L = sp.optimize.brentq(f, lo, hi)
        rates[treatment] = L / beta  # per-division loss probability

    return rates
