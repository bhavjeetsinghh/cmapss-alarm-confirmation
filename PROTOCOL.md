# Paper 4 protocol: independent confirmation on C-MAPSS turbofan run-to-failure data
Fixed 2026-10-14, before any C-MAPSS result was computed (only file row and unit counts were checked).
Hypotheses come from the bearing studies (papers 2-3); this is a confirmatory test on a different system.

## Data
NASA C-MAPSS training sets FD001-FD004 (Saxena, Goebel, Simon & Eklund, 2008): 709 simulated turbofan engines run to failure
(FD001 100, FD002 260, FD003 100, FD004 249). Test sets are not used (they end before failure). One snapshot per cycle; life = last cycle.
Sensors: the 14 non-constant sensors commonly used in the literature: 2, 3, 4, 7, 8, 9, 11, 12, 13, 14, 15, 17, 20, 21.

## Health indicator (label-free)
1. Operating regime: FD001/FD003 one regime; FD002/FD004 regime = operational setting 1 rounded to the nearest integer (six regimes).
2. Regime normalization (fleet history, early life only): for each subset, regime and sensor, mean and sd over cycles 11-60 of all units; each reading is standardized with these.
3. Per unit, run-in = cycles 1-10, baseline = cycles 11-60 (n = 50). z_jt = (x_jt - m_j) / s_j with m_j, s_j the unit's baseline mean and sd of normalized sensor j.
4. D_t = sqrt(mean_j z_jt^2) (healthy D is about 1). Monitoring starts at cycle 61. All units have at least 128 cycles, so all are evaluable.

## Rules (same structure as the bearing studies; thresholds in D units)
Default trigger c = 2.0, with ratios carried over from the bearing rules.
- M0 reference: D > 2.0 for k = 3 consecutive cycles.
- M1 m-of-n: >= 8 of the last 10 cycles with D > 2.0.
- M2 median: causal median of D over the last 10 cycles > 2.0, k = 3.
- M3 EWMA: z_t = 0.1 D_t + 0.9 z_{t-1} (z = 1 at baseline end) > 2.0, k = 3.
- M4 hysteresis: on after 3 cycles with D > 2.0; off after 3 cycles with D < 2.0 x 1.1/1.2 (= 1.833).
- M6 median + hysteresis: M2 filter, then M4 logic.
- M7 CUSUM: S_t = max(0, S_{t-1} + ln D_t - ln(2.0 x 1.1/1.2)) > 5 ln 1.2, k = 3.
(M5 axis voting has no analogue: one multivariate indicator per engine.)

## Outcomes (as paper 3)
Detection: sustained alarm block reaching the last cycle (gaps <= 1% of life) of at least 1% of life; k/n with Clopper-Pearson 95% CI.
Nuisance episodes: alarm episodes that end before the last cycle (total, engines with >= 1).
Early alarm: first alarm before 85% of life. Median sustained warning in cycles and as a fraction of life.
Paired vs M0: Wilcoxon signed-rank (nuisance per engine), exact McNemar (detection).
Empirical cost (pessimistic: replace at sustained-warning start; optimistic: just in time), costs and E[L] scaling as papers 2-3.

## Confirmatory hypotheses (from the bearing studies)
H1 Hysteresis (M4) reduces nuisance episodes vs M0 (Wilcoxon p < 0.05) and loses no detection (McNemar lost = 0 by construction; reported).
H2 Filters (M1, M2, M3, M6) each detect fewer failures than M0.
H3 Leave-one-subset-out: the configuration (method + parameters, lowest pessimistic cost on three subsets) does not beat untuned M4 on the held-out subset in more than 2 of 4 folds.
H4 Watch-and-act (paper 3 secondary analysis): replacing at the first M4 alarm (A0) costs less than waiting for escalation D > 2.0 x 2.0/1.2 (= 3.33) (A1).
A hypothesis is "confirmed" if its stated direction holds; all results are reported regardless.

## Grid (for H3 and the operating curve)
c in {1.25, 1.50, ..., 4.00}; method parameters as in paper 3: m in {6, 8, 10}; window in {5, 10, 20}; lambda in {0.05, 0.1, 0.2}; off/on in {0.85, 0.90, 0.95}; H in {3, 5, 10} x ln 1.2.

## Watch-and-act (H4)
WATCH = M4. ACT latches. A0 first WATCH cycle; A1 D > 3.33 for 3 cycles; A2 one-sided Mann-Kendall on ln D over the last 30 cycles, p < 0.01, 3 cycles.
In time = ACT at least 1% of life before the last cycle. Cost as paper 3 secondary analysis.

## Sensitivity
Planning lead/detection threshold 5 and 10 cycles (absolute) instead of 1% of life; k = 5; health indicator Hotelling T^2 on the 14 z-scores (baseline covariance, pseudo-inverse) in place of D with c chosen as the F-based UCL at alpha = 1e-4 and hysteresis off at 0.917 x UCL; Cfa = 1,000 and 5,000.
