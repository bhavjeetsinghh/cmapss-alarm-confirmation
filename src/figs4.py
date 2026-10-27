"""Figures for paper 4. Usage: python figs4.py <C-MAPSS dir> <cmapss.json> <paper3.json> <out dir>"""
import sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cmapss_lib as cm
DATA, RJ, P3J, OUT = sys.argv[1:5]
E = cm.load(DATA); R = json.load(open(RJ)); P3 = json.load(open(P3J))
BLUE, ORANGE, AQUA, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#1baf7a', '#222222', '#6b6b6b', '#e4e4e0'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'Liberation Serif', 'DejaVu Serif'], 'font.size': 8, 'axes.titlesize': 8,
                     'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.edgecolor': MUTED,
                     'axes.linewidth': 0.6, 'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.5, 'savefig.dpi': 600, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
W1, W2 = 3.4, 7.0
def save(fig, name): fig.savefig(f'{OUT}/{name}.png', bbox_inches='tight'); fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight'); plt.close(fig)

# Fig 1: health indicator D for 25 engines per subset
fig, axs = plt.subplots(1, 4, figsize=(W2, 2.2), sharey=True)
rng = np.random.default_rng(0)
for ax, sub in zip(axs, cm.SUBSETS):
    es = [e for e in E if e['dataset'] == sub]; pick = rng.choice(len(es), 25, replace=False)
    for i in pick:
        e = es[i]; ax.plot(e['t'] / e['life'], e['axes'][0], lw=0.4, color=BLUE, alpha=0.5)
    ax.axhline(2.0, color=INK, ls=':', lw=0.7); ax.axhline(2.0 * 1.1 / 1.2, color=INK, ls='--', lw=0.6)
    ax.set_yscale('log'); ax.set_yticks([0.5, 1, 2, 5, 10]); ax.set_yticklabels(['0.5', '1', '2', '5', '10']); ax.minorticks_off(); ax.set_ylim(0.4, 15)
    ax.set_title(f'{sub} ({len(es)} engines)', loc='left'); ax.set_xlabel('Fraction of life')
axs[0].set_ylabel('Health indicator D')
save(fig, 'fig1_indicator')

# Fig 2: cross-system comparison of alarm logic (defaults)
ms = ['M1', 'M2', 'M3', 'M4', 'M6', 'M7']; lab = ['m-of-n', 'Median', 'EWMA', 'Hysteresis', 'Median+hyst.', 'CUSUM']
b0 = P3['main']['M0']['nuis_total']; e0 = R['main']['M0']['nuis_total']
bn = [100 * P3['main'][m]['nuis_total'] / b0 for m in ms]; en = [100 * R['main'][m]['nuis_total'] / e0 for m in ms]
bl = [100 * (P3['main']['M0']['detected'] - P3['main'][m]['detected']) / 33 for m in ms]; el = [100 * (R['main']['M0']['detected'] - R['main'][m]['detected']) / 709 for m in ms]
fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 2.2))
x = np.arange(len(ms)); wd = 0.38
a.bar(x - wd / 2, bn, wd, color=ORANGE, label='bearings (33)'); a.bar(x + wd / 2, en, wd, color=BLUE, label='engines (709)')
a.set_xticks(x); a.set_xticklabels(lab, rotation=30, ha='right'); a.set_ylabel('Nuisance episodes (% of reference)'); a.legend(); a.grid(axis='x', visible=False); a.set_title('(a) Nuisance alarms', loc='left')
b.bar(x - wd / 2, bl, wd, color=ORANGE, label='bearings (33)'); b.bar(x + wd / 2, el, wd, color=BLUE, label='engines (709)')
b.set_xticks(x); b.set_xticklabels(lab, rotation=30, ha='right'); b.set_ylabel('Net failures missed (% of failures)'); b.legend(); b.grid(axis='x', visible=False); b.set_title('(b) Detection lost vs reference', loc='left')
b.set_ylim(-1, 14); b.text(0.97, 0.06, 'engines: no rule missed a failure', transform=b.transAxes, ha='right', fontsize=7, color=BLUE)
save(fig, 'fig2_cross_system')

# Fig 3: detection vs nuisance, grid configurations (engines)
g = pd.DataFrame(R['grid_summary']); M = R['main']
fig, ax = plt.subplots(figsize=(W1, 2.6))
ax.scatter(g.nuis_total + 1, g.detected, s=5, color='#bdbdb8', lw=0, label='grid configurations')
for m in ['M0', 'M1', 'M2', 'M3', 'M4', 'M6', 'M7']:
    s = M[m]; col = ORANGE if m == 'M0' else (BLUE if m == 'M4' else INK)
    ax.scatter(s['nuis_total'] + 1, s['detected'], s=20, color=col, zorder=3)
    ax.annotate(m, (s['nuis_total'] + 1, s['detected']), xytext=(-4, 6 if m in ('M0', 'M1', 'M4', 'M7') else -11), textcoords='offset points', fontsize=7, color=col)
ax.set_xscale('log'); ax.set_xlabel('Nuisance episodes, 709 engines (+1, log)'); ax.set_ylabel('Failures detected (of 709)'); ax.legend(loc='lower left')
save(fig, 'fig3_tradeoff')
print('ok')
