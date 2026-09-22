import json, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import *

HOR = [(-60,-50),(-30,-25),(-15,-10),(-10,-5),(-5,-2),(-2,-1),(-1,0)]

res = {}
tot = dict(events_in_cache=0, onset_resolved=0, onset_miss=0)
for aid in animals():
    D = order_by_time(load(aid))
    if len(D["ev"]) == 0:
        res[aid] = dict(n_events=0, n_epochs=int(len(D["t"])), note="no annotated sessions")
        continue
    ons, miss = event_onsets(D)
    tot["events_in_cache"] += len(D["ev"]); tot["onset_resolved"] += len(ons); tot["onset_miss"] += miss
    ids, lead4 = lead_mask(ons, 4*HOUR)
    _, lead2 = lead_mask(ons, 2*HOUR)
    _, lead1 = lead_mask(ons, 1*HOUR)
    pool = interictal_pool(D, ons)
    # onset-epoch blanking
    nblank_onset = sum(1 for k in ons if not ons[k]["onset_epoch_valid"])
    lab = {}
    for k in ons:
        lab[ons[k]["label"]] = lab.get(ons[k]["label"], 0) + 1
    res[aid] = dict(
        n_epochs=int(len(D["t"])), n_valid=int(D["valid"].sum()),
        hours=round(float(D["dur"].sum())/3600, 1),
        n_events=len(D["ev"]), n_onset_resolved=len(ons), n_onset_miss=miss,
        n_lead_1h=int(sum(lead1.values())), n_lead_2h=int(sum(lead2.values())),
        n_lead_4h=int(sum(lead4.values())),
        onset_epoch_blanked=nblank_onset,
        pool_epochs=int(pool.sum()), pool_hours=round(float(D["dur"][pool].sum())/3600, 1),
        pool_hours_present=int(len(np.unique(hod(D["t"][pool])))),
        labels=lab,
    )
    print(aid, res[aid], flush=True)
json.dump(dict(per_animal=res, totals=tot), open("cohort.json","w"), indent=1)
print("TOTALS", tot)
print("events total", sum(v.get("n_events",0) for v in res.values()))
for L in ("n_lead_1h","n_lead_2h","n_lead_4h"):
    print(L, sum(v.get(L,0) for v in res.values()))
