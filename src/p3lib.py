"""Paper 3 library: data loading, alarm rules M0-M7, scoring and empirical cost (see PROTOCOL.md)."""
import sys, json, itertools, numpy as np, pandas as pd
from scipy.stats import beta as Beta, wilcoxon, binomtest

RUNIN, BASE_END, MIN_SNAP, FINAL, MINLEAD = 10, 60, 120, 0.15, 0.01
Cp, Cf, Cfa0, EL = 14500.0, 88000.0, 2500.0, 1053.2
AGE = {'1.5': 71227.37673442051, '2.0': 55859.9366910734, '3.0': 40241.562925333004}   # paper 2, results/pooled.json
RTF = 1000 * Cf / EL

def cp_ci(k, n, a=0.05):
    lo = 0.0 if k == 0 else Beta.ppf(a / 2, k, n - k + 1); hi = 1.0 if k == n else Beta.ppf(1 - a / 2, k + 1, n - k)
    return [float(lo), float(hi)]

def persist(b, k):
    if k <= 1: return b.copy()
    c = np.convolve(b.astype(int), np.ones(k, int))[:len(b)]
    return c >= k

# ---------------- data: bearing -> dict(t, frac_left, life, axes: list of ratio arrays)
def load_failed(POOLED):
    u = pd.read_csv(POOLED, low_memory=False)
    B = []
    for (ds, unit), g in u.groupby(['dataset', 'unit'], sort=False):
        if g.n_snap.iloc[0] < MIN_SNAP: continue
        axes = []
        for ax, s in g.groupby('axis', sort=False):
            s = s.sort_values('idx'); r = s.rms.values
            axes.append(r / r[RUNIN:BASE_END].mean())
        s = g[g.axis == g.axis.iloc[0]].sort_values('idx')
        n = min(len(a) for a in axes)
        B.append(dict(dataset=ds, unit=unit, t=s.t_h.values[:n], life=float(s.life_h.iloc[0]), ttg=s.ttg_frac.values[:n],
                      axes=[a[:n] for a in axes]))
    return B

def load_ims_nonfailed(IMSD):
    out = []
    spec = {'T1': ('test1.csv', [1, 2]), 'T2': ('test2.csv', [2, 3, 4]), 'T3': ('test3.csv', [1, 2, 4])}
    for run, (f, bears) in spec.items():
        d = pd.read_csv(f'{IMSD}/{f}', parse_dates=['timestamp']).sort_values(['timestamp', 'channel'])
        med = d.groupby('channel').rms.median(); dead = (d.rms < 0.1 * d.channel.map(med)).groupby(d.filename).all()
        bad = []
        for fn in sorted(dead.index, reverse=True):
            if dead[fn]: bad.append(fn)
            else: break
        d = d[~d.filename.isin(bad)]
        ts = pd.Series(sorted(d.timestamp.unique())); dt = ts.diff().dt.total_seconds().fillna(0) / 3600
        clock = pd.DataFrame({'timestamp': ts, 't_h': np.where(dt > 0.25, 10 / 60, dt).cumsum()}); d = d.merge(clock, on='timestamp')
        for b in bears:
            axes = []
            for ch, s in d[d.bearing == b].groupby('channel'):
                s = s.sort_values('t_h'); r = s.rms.values; axes.append(r / r[RUNIN:BASE_END].mean()); t = s.t_h.values
            life = t[-1]
            out.append(dict(dataset='IMS-nonfailed', unit=f'{run}-B{b}', t=t, life=life, ttg=(life - t) / life, axes=axes))
    return out

# ---------------- rules: ratio array -> boolean alarm (monitoring only)
def mon_mask(n): m = np.zeros(n, bool); m[BASE_END:] = True; return m

def rule(r, method, c=1.2, k=3, m=8, w=10, lam=0.1, off=1.1 / 1.2, H=5):
    n = len(r); mon = mon_mask(n)
    if method in ('M0', 'M5'):
        a = persist(r > c, k)
    elif method == 'M1':
        e = (r > c).astype(int); a = np.convolve(e, np.ones(10, int))[:n] >= m
    elif method == 'M2':
        f = pd.Series(r).rolling(w, min_periods=1).median().values; a = persist(f > c, k)
    elif method == 'M3':
        z = np.ones(n); zz = 1.0
        for i in range(BASE_END, n): zz = lam * r[i] + (1 - lam) * zz; z[i] = zz
        a = persist(z > c, k)
    elif method in ('M4', 'M6'):
        x = pd.Series(r).rolling(w, min_periods=1).median().values if method == 'M6' else r
        on, offc = persist(x > c, k), persist(x < c * off, k)
        a = np.zeros(n, bool); st = False
        for i in range(BASE_END, n):
            if not st and on[i]: st = True
            elif st and offc[i]: st = False
            a[i] = st
    elif method == 'M7':
        K = np.log(c * 1.1 / 1.2); h = H * np.log(1.2); S = np.zeros(n); s = 0.0
        lr = np.log(np.maximum(r, 1e-9))
        for i in range(BASE_END, n): s = max(0.0, s + lr[i] - K); S[i] = s
        a = persist(S > h, k)
    return a & mon

def bearing_signal(b, method, **kw):
    al = [rule(r, method, **kw) for r in b['axes']]
    return np.logical_and.reduce(al) if method == 'M5' else np.logical_or.reduce(al)

def score(b, a, minlead=MINLEAD):
    t, life, n = b['t'], b['life'], len(a)
    rising = a & ~np.r_[False, a[:-1]]
    idx = np.where(a)[0]; sus_i = None
    if len(idx) and idx[-1] >= n - 3:
        j = len(idx) - 1
        while j > 0 and (t[idx[j]] - t[idx[j - 1]]) <= 0.01 * life + 1e-12: j -= 1
        sus_i = idx[j]
    sus_frac = (life - t[sus_i]) / life if sus_i is not None else 0.0
    det = sus_i is not None and sus_frac >= minlead
    edges = np.where(rising)[0]
    # protocol: every episode that ends before the end of the record (only an episode still on at the last snapshot is not nuisance)
    nuis = len(edges) - int(n > 0 and bool(a[-1]))
    ref = (np.arange(n) >= BASE_END) & (b['ttg'] > FINAL)
    first = float(t[edges[0]] / life) if len(edges) else None
    return dict(det=bool(det), sus_frac=float(sus_frac) if det else 0.0, nuis=nuis, hz_fa=int((rising & ref).sum()),
                first=first, early_alarm=bool(first is not None and first < 1 - FINAL))

def cost(rows, cfa=Cfa0, mode='pess'):
    num = sum((Cp if r['det'] else Cf) + r['nuis'] * cfa for r in rows)
    den = sum(EL * (1 - (r['sus_frac'] if mode == 'pess' else MINLEAD) * r['det']) for r in rows)
    return 1000 * num / den

def summarize(rows):
    n = len(rows); d = sum(r['det'] for r in rows); nu = [r['nuis'] for r in rows]
    leads = [r['sus_frac'] for r in rows if r['det']]
    return dict(n=n, detected=d, ci=cp_ci(d, n), nuis_mean=float(np.mean(nu)), nuis_median=float(np.median(nu)),
                nuis_bearings=int(sum(x > 0 for x in nu)), nuis_total=int(sum(nu)), hz_fa_bearings=int(sum(r['hz_fa'] > 0 for r in rows)),
                early_alarm_bearings=int(sum(r['early_alarm'] for r in rows)),
                median_lead=float(np.median(leads)) if leads else None, cost_pess=cost(rows), cost_opt=cost(rows, mode='opt'))

def evaluate(B, method, minlead=MINLEAD, **kw):
    return [dict(dataset=b['dataset'], unit=b['unit'], **score(b, bearing_signal(b, method, **kw), minlead)) for b in B]

