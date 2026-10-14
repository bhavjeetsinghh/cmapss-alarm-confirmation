# Do bearing alarm findings transfer to turbofan engines?

Code and data for the paper *Which Alarm-Logic Findings Transfer from Bearings to Turbofan Engines? A Confirmatory, Label-Free Test on 709 C-MAPSS Run-to-Failure Trajectories* (Bhavjeet Singh Hora, Florida International University).

Four findings from label-free studies of 33 run-to-failure bearings are tested on the 709 simulated turbofan engines of the NASA C-MAPSS training sets, with the health indicator, rules, outcomes and hypotheses fixed in advance ([`PROTOCOL.md`](PROTOCOL.md)).

| Hypothesis (from the bearing studies) | Result on 709 engines |
|---|---|
| H1 Hysteresis removes nuisance alarms without missing failures | Confirmed: 874 → 38 nuisance episodes, all 709 failures detected |
| H2 Smoothing filters miss failures | Not confirmed: no filter missed a failure (gradual degradation) |
| H3 Tuned thresholds do not transfer | Not confirmed: tuned settings beat the default in 4 of 4 folds (one engine model) |
| H4 Replace at the first alarm rather than wait for escalation | Confirmed: waiting let 58 engines fail and cost 21% more |

Companion repositories: [ims-bearing-label-free-early-warning](https://github.com/bhavjeetsinghh/ims-bearing-label-free-early-warning), [bearing-early-warning-transfer](https://github.com/bhavjeetsinghh/bearing-early-warning-transfer), [bearing-nuisance-alarms](https://github.com/bhavjeetsinghh/bearing-nuisance-alarms).

## Reproduce
```bash
pip install numpy pandas scipy matplotlib
bash run_all.sh
```
Outputs: `results/cmapss.json` (all reported numbers) and `figures/fig1–3` (PNG and PDF).

## Data
- `data/CMAPSSData/` — C-MAPSS training sets FD001–FD004 and readme (Saxena, Goebel, Simon & Eklund, 2008; NASA Ames Prognostics Data Repository). The test sets are not used.
- `data/bearing_paper3_results.json` — results of the bearing nuisance-alarm study, used only for the cross-system comparison in Figure 2.

## Files
- `src/cmapss_lib.py` — loading, regime normalization, health indicator *D*, Hotelling *T*², scoring
- `src/p3lib.py` — alarm rules M0–M7 and scoring shared with the bearing study
- `src/cmapss.py` — main analysis, hypotheses H1–H4, sensitivity, leave-one-subset-out selection
- `src/figs4.py` — figures

## License
MIT for the code. The C-MAPSS data remain under NASA's terms.
