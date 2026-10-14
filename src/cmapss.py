"""Paper 4: confirmation of the bearing alarm findings on NASA C-MAPSS (protocol in PROTOCOL_CMAPSS.md).
Usage: python cmapss.py <C-MAPSS dir with train_FD00x.txt> <out dir>"""
import sys, json, numpy as np
from scipy.stats import f as Fdist
from cmapss_lib import *
DATA, OUT = sys.argv[1:3]
E = load(DATA)
print('engines', len(E), {s: sum(e['dataset'] == s for e in E) for s in SUBSETS}, 'min life', min(e['life'] for e in E), 'sensors used', sorted(set(e['p'] for e in E)), flush=True)
C0 = 2.0
METHODS = ['M0', 'M1', 'M2', 'M3', 'M4', 'M6', 'M7']
R = dict(protocol='PROTOCOL_CMAPSS.md', n=len(E), costs=dict(Cp=Cp, Cf=Cf, Cfa=Cfa0, EL=EL, age_opt=AGE, rtf=1000 * Cf / EL))
healthy = np.concatenate([e['axes'][0][BASE_END:int(0.5 * e['life'])] for e in E])
R['healthy_D'] = dict(median=float(np.median(healthy)), p95=float(np.percentile(healthy, 95)), p99=float(np.percentile(healthy, 99)))
main, RW = {}, {}
for m in METHODS:
    rows = evaluate(E, m, c=C0); RW[m] = rows; s = summarize(rows)
    s['by_subset'] = {sub: summarize([r for r in rows if r['dataset'] == sub]) for sub in SUBSETS}
    if m != 'M0': s['paired'] = paired(RW['M0'], rows)
    s['units'] = rows; main[m] = s
    print(m, s['detected'], [round(x, 3) for x in s['ci']], 'nuis', s['nuis_total'], s['nuis_units'], 'early', s['early_alarm_units'], 'lead', s['median_lead_cycles'],
          round(s['cost_pess']), round(s['cost_opt']), s.get('paired'), flush=True)
R['main'] = main
R['hypotheses'] = dict(
    H1=dict(nuis_M0=main['M0']['nuis_total'], nuis_M4=main['M4']['nuis_total'], p=main['M4']['paired']['wilcoxon_p'], lost=main['M4']['paired']['lost'],
            confirmed=bool(main['M4']['nuis_total'] < main['M0']['nuis_total'] and main['M4']['paired']['wilcoxon_p'] < 0.05 and main['M4']['paired']['lost'] == 0)),
    H2={m: dict(detected=main[m]['detected'], M0=main['M0']['detected'], confirmed=bool(main[m]['detected'] < main['M0']['detected'])) for m in ['M1', 'M2', 'M3', 'M6']})

# sensitivity
sens = {}
sens['abs5'] = {m: summarize(evaluate(E, m, abs_lead=5, c=C0)) for m in METHODS}
sens['abs10'] = {m: summarize(evaluate(E, m, abs_lead=10, c=C0)) for m in METHODS}
sens['k5'] = {m: summarize(evaluate(E, m, k=5, c=C0)) for m in METHODS}
p, n = 14, 50; ucl = p * (n + 1) * (n - 1) / (n * (n - p)) * Fdist.ppf(1 - 1e-4, p, n - p)
sens['T2'] = dict(ucl=ucl, M0=summarize(evaluate(E, 'M0', sig_key='T2', c=ucl)), M4=summarize(evaluate(E, 'M4', sig_key='T2', c=ucl, off=1.1 / 1.2)))
sens['T2']['paired'] = paired(evaluate(E, 'M0', sig_key='T2', c=ucl), evaluate(E, 'M4', sig_key='T2', c=ucl, off=1.1 / 1.2))
for cfa in (1000.0, 5000.0): sens[f'cfa{int(cfa)}'] = {m: summarize(RW[m], cfa) for m in METHODS}
R['sensitivity'] = sens
for lab, v in sens.items():
    if lab != 'T2': print(lab, {m: (x['detected'], x['nuis_total'], round(x['cost_pess'])) for m, x in v.items()}, flush=True)
print('T2', round(ucl, 2), {m: (sens['T2'][m]['detected'], sens['T2'][m]['nuis_total']) for m in ('M0', 'M4')}, sens['T2']['paired'], flush=True)

# grid + leave-one-subset-out (H3)
CG = [1.25 + 0.25 * i for i in range(12)]
grid = []
for c in CG:
    grid += [('M0', dict(c=c))] + [('M1', dict(c=c, m=v)) for v in (6, 8, 10)] + [('M2', dict(c=c, w=v)) for v in (5, 10, 20)]
    grid += [('M3', dict(c=c, lam=v)) for v in (0.05, 0.1, 0.2)] + [('M4', dict(c=c, off=v)) for v in (0.85, 0.90, 0.95)]
    grid += [('M6', dict(c=c, off=v)) for v in (0.85, 0.90, 0.95)] + [('M7', dict(c=c, H=v)) for v in (3, 5, 10)]
G = [dict(method=m, kw=kw, rows=evaluate(E, m, **kw)) for m, kw in grid]
print('grid', len(G), flush=True)
def sub(rows, ss): return [r for r in rows if r['dataset'] in ss]
lodo = {}; wins = 0
for hold in SUBSETS:
    tr = [s for s in SUBSETS if s != hold]; best = min(G, key=lambda g: cost(sub(g['rows'], tr)))
    t = summarize(sub(best['rows'], [hold])); d4 = summarize(sub(RW['M4'], [hold]))
    beat = t['cost_pess'] < d4['cost_pess']; wins += int(beat)
    lodo[hold] = dict(method=best['method'], kw=best['kw'], test=t, untuned_M4=d4, untuned_M0=summarize(sub(RW['M0'], [hold])), tuned_beats_M4=bool(beat))
R['lodo'] = lodo; R['hypotheses']['H3'] = dict(tuned_wins=wins, folds=4, confirmed=bool(wins <= 2))
R['grid_summary'] = [dict(method=g['method'], **g['kw'], **{k: v for k, v in summarize(g['rows']).items() if k != 'ci'}) for g in G]
print('LODO', {h: (v['method'], v['kw'], v['test']['detected'], v['test']['n'], round(v['test']['cost_pess']), round(v['untuned_M4']['cost_pess']), v['tuned_beats_M4']) for h, v in lodo.items()}, flush=True)

# watch-and-act (H4)
def act_rows(E, act, c_act=2.0 * 2.0 / 1.2, alpha=0.01, w=30, lead=0.01):
    out = []
    for e in E:
        D = e['axes'][0]; n = len(D); W = rule(D, 'M4', c=C0)
        if act == 'A0': A = W.copy()
        elif act == 'A1': A = persist(D > c_act, 3) & W
        else: A = persist(np.nan_to_num(mk_p(np.log(np.maximum(D, 1e-9)), w), nan=1.0) < alpha, 3) & W
        A[:BASE_END] = False; rise = np.where(W & ~np.r_[False, W[:-1]])[0]; ia = np.where(A)[0]
        if len(ia):
            ta = e['t'][ia[0]]; it = bool(ta <= e['life'] * (1 - lead) + 1e-12); cur = rise[rise <= ia[0]]
            out.append(dict(dataset=e['dataset'], intime=it, wasted=float((e['life'] - ta) / e['life']) if it else 0.0, nins=(len(cur) - 1) if it else len(rise)))
        else: out.append(dict(dataset=e['dataset'], intime=False, wasted=0.0, nins=len(rise)))
    return out
def act_sum(rows):
    k = sum(r['intime'] for r in rows); wl = [r['wasted'] for r in rows if r['intime']]
    num = sum((Cp if r['intime'] else Cf) + r['nins'] * Cfa0 for r in rows); den = sum(EL * (1 - r['wasted'] * r['intime']) for r in rows)
    return dict(n=len(rows), intime=k, ci=cp_ci(k, len(rows)), median_wasted=float(np.median(wl)) if wl else None, inspections=int(sum(r['nins'] for r in rows)), cost=1000 * num / den)
wa = {a: act_sum(act_rows(E, a)) for a in ('A0', 'A1', 'A2')}
R['watch_act'] = wa; R['hypotheses']['H4'] = dict(A0=wa['A0']['cost'], A1=wa['A1']['cost'], confirmed=bool(wa['A0']['cost'] < wa['A1']['cost']))
print('watch-act', {a: (v['intime'], round(v['median_wasted'] or 0, 3), v['inspections'], round(v['cost'])) for a, v in wa.items()}, flush=True)
print('HYPOTHESES', json.dumps({h: v.get('confirmed') if 'confirmed' in v else {m: x['confirmed'] for m, x in v.items()} for h, v in R['hypotheses'].items()}))
print('healthy D', R['healthy_D'])
json.dump(R, open(f'{OUT}/cmapss.json', 'w'), indent=1, default=float)
