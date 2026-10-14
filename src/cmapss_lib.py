"""Paper 4 library: C-MAPSS loading, health indicator, rules and scoring (see PROTOCOL_CMAPSS.md)."""
import sys, json, numpy as np, pandas as pd
from scipy.stats import norm, wilcoxon, binomtest, f as Fdist
from p3lib import rule, persist, score, cost, cp_ci, BASE_END
RUNIN = 10
SENS = [2, 3, 4, 7, 8, 9, 11, 12, 13, 14, 15, 17, 20, 21]
SUBSETS = ['FD001', 'FD002', 'FD003', 'FD004']
Cp, Cf, Cfa0, EL = 14500.0, 88000.0, 2500.0, 1053.2
AGE = {'1.5': 71227.37673442051, '2.0': 55859.9366910734, '3.0': 40241.562925333004}

def load(DATA):
    E = []
    for sub in SUBSETS:
        cols = ['unit', 'cycle', 's1', 's2', 's3'] + [f'x{i}' for i in range(1, 22)]
        d = pd.read_csv(f'{DATA}/train_{sub}.txt', sep=r'\s+', header=None, names=cols)
        d['regime'] = d.s1.round().astype(int) if sub in ('FD002', 'FD004') else 0
        X = [f'x{i}' for i in SENS]
        early = d[(d.cycle > RUNIN) & (d.cycle <= BASE_END)]
        mu = early.groupby('regime')[X].mean(); sd = early.groupby('regime')[X].std().replace(0, np.nan)
        Z = (d[X] - mu.loc[d.regime].values) / sd.loc[d.regime].values
        for u, g in d.groupby('unit'):
            z = Z.loc[g.index].values; base = z[RUNIN:BASE_END]
            m, s = np.nanmean(base, 0), np.nanstd(base, 0, ddof=1)
            ok = np.isfinite(s) & (s > 0)
            zz = (z[:, ok] - m[ok]) / s[ok]
            D = np.sqrt(np.nanmean(zz ** 2, 1))
            # Hotelling T2 (sensitivity)
            S = np.cov(zz[RUNIN:BASE_END].T); Si = np.linalg.pinv(S); c0 = zz - zz[RUNIN:BASE_END].mean(0)
            T2 = np.einsum('ij,jk,ik->i', c0, Si, c0)
            life = float(g.cycle.max()); t = g.cycle.values.astype(float)
            E.append(dict(dataset=sub, unit=int(u), t=t, life=life, ttg=(life - t) / life, axes=[D], T2=T2, p=int(ok.sum())))
    return E

def mk_p(x, w):
    n = len(x); p = np.full(n, np.nan)
    if n < w: return p
    W = np.lib.stride_tricks.sliding_window_view(x, w); s = np.zeros(W.shape[0])
    for i in range(w - 1): s += np.sign(W[:, i + 1:] - W[:, [i]]).sum(axis=1)
    var = w * (w - 1) * (2 * w + 5) / 18.0
    z = np.where(s > 0, (s - 1) / np.sqrt(var), np.where(s < 0, (s + 1) / np.sqrt(var), 0.0))
    p[w - 1:] = 1 - norm.cdf(z); return p

def alarm(e, method, sig=None, **kw):
    r = e['axes'][0] if sig is None else sig
    return rule(r, method, **kw)

def evaluate(E, method, minlead=0.01, abs_lead=None, sig_key=None, k=3, **kw):
    rows = []
    for e in E:
        sig = e[sig_key] if sig_key else None
        a = alarm(e, method, sig, k=k, **kw)
        s = score(e, a, 0.0 if abs_lead else minlead)
        if abs_lead:
            frac = s['sus_frac']; s['det'] = bool(s['det'] and frac * e['life'] >= abs_lead); s['sus_frac'] = frac if s['det'] else 0.0
        rows.append(dict(dataset=e['dataset'], unit=e['unit'], life=e['life'], **s))
    return rows

def summarize(rows, cfa=Cfa0):
    n = len(rows); d = sum(r['det'] for r in rows); nu = [r['nuis'] for r in rows]
    leads = [r['sus_frac'] for r in rows if r['det']]; leadc = [r['sus_frac'] * r['life'] for r in rows if r['det']]
    return dict(n=n, detected=d, ci=cp_ci(d, n), nuis_total=int(sum(nu)), nuis_units=int(sum(x > 0 for x in nu)), nuis_mean=float(np.mean(nu)),
                early_alarm_units=int(sum(r['early_alarm'] for r in rows)), median_lead_frac=float(np.median(leads)) if leads else None,
                median_lead_cycles=float(np.median(leadc)) if leadc else None, cost_pess=cost(rows, cfa), cost_opt=cost(rows, cfa, 'opt'))

def paired(r0, r1):
    a0 = np.array([r['nuis'] for r in r0]); a1 = np.array([r['nuis'] for r in r1])
    d0 = np.array([r['det'] for r in r0]); d1 = np.array([r['det'] for r in r1]); lost, gained = int((d0 & ~d1).sum()), int((~d0 & d1).sum())
    return dict(wilcoxon_p=float(wilcoxon(a0, a1).pvalue) if np.any(a0 != a1) else 1.0, lost=lost, gained=gained,
                mcnemar_p=float(binomtest(lost, lost + gained).pvalue) if lost + gained else 1.0)

