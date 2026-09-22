import json, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import *

HOR = [(-60,-50),(-30,-25),(-15,-10),(-10,-5),(-5,-2),(-2,-1),(-1,0)]
HNAME = {h: f"{-h[0]}..{-h[1]} min" for h in HOR}
POST = [(0,1),(1,2),(2,5),(5,10),(10,20),(20,40),(40,60),(60,120)]
MIN_COVER = 0.5
MIN_VALID = 3
MIN_CTRL_PER_HOUR = 30
LEAD_H = 4.0

def window_stats(D, t0, t1, sess):
    """epochs in [t0,t1) within session `sess`. returns n_present, n_valid, median F, ."""
    t = D["t"]
    lo = np.searchsorted(t, t0, "left"); hi = np.searchsorted(t, t1, "left")
    if hi <= lo:
        return 0, 0, None
    sl = slice(lo, hi)
    m = (D["sess"][sl] == sess)
    npres = int(m.sum())
    v = m & D["valid"][sl]
    nval = int(v.sum())
    if nval == 0:
        return npres, 0, None
    return npres, nval, np.median(D["F"][sl][v], axis=0)

out = {"per_animal": {}, "windows": [], "post": [], "diag": {}}
diag = dict(win_total=0, win_no_epochs=0, win_low_cover=0, win_no_valid=0,
            win_ok=0, win_epochs_present=0, win_epochs_valid=0,
            win_no_matched_ctrl=0, pool_hours_missing=0)

cohort = json.load(open("cohort.json"))["per_animal"]
PRIMARY = [a for a, v in cohort.items() if v.get("n_lead_4h", 0) >= 10]
print("primary animals:", len(PRIMARY), sorted(PRIMARY))

for aid in animals():
    if cohort.get(aid, {}).get("n_events", 0) == 0:
        continue
    D = order_by_time(load(aid))
    ons, _ = event_onsets(D)
    ids, lead = lead_mask(ons, LEAD_H*HOUR)
    pool = interictal_pool(D, ons)
    ph = hod(D["t"][pool]); PF = D["F"][pool]
    ref_h = {}
    for h in range(24):
        m = ph == h
        if m.sum() >= MIN_CTRL_PER_HOUR:
            ref_h[h] = np.median(PF[m], axis=0)
    ref_all = np.median(PF, axis=0) if pool.sum() >= MIN_CTRL_PER_HOUR else None
    mad = np.median(np.abs(PF - ref_all), axis=0)*1.4826 if ref_all is not None else None
    if mad is not None:
        mad = np.where(mad <= 0, np.nan, mad)
    out["per_animal"][aid] = dict(pool_n=int(pool.sum()),
                                  hours_with_ctrl=sorted(ref_h),
                                  sd=None if mad is None else [float(x) for x in mad],
                                  n_lead=int(sum(lead.values())))
    # ---------- pre-ictal windows
    for k in ids:
        if not lead[k]:
            continue
        e = ons[k]
        for (a, b) in HOR:
            t0 = e["onset"] + a*60.0; t1 = e["onset"] + b*60.0
            exp = (b-a)*6.0
            npres, nval, med = window_stats(D, t0, t1, e["session"])
            diag["win_total"] += 1
            diag["win_epochs_present"] += npres; diag["win_epochs_valid"] += nval
            if npres == 0: diag["win_no_epochs"] += 1; continue
            if npres/exp < MIN_COVER: diag["win_low_cover"] += 1; continue
            if nval < MIN_VALID: diag["win_no_valid"] += 1; continue
            hh = int(hod(np.array([(t0+t1)/2]))[0])
            rm = ref_h.get(hh)
            if rm is None: diag["win_no_matched_ctrl"] += 1
            diag["win_ok"] += 1
            out["windows"].append(dict(a=aid, ev=k, hor=f"{a},{b}", hod=hh, n=nval, npres=npres,
                                       stage=e["label"],
                                       v=[float(x) for x in med],
                                       rm=None if rm is None else [float(x) for x in rm],
                                       ru=None if ref_all is None else [float(x) for x in ref_all],
                                       sd=None if mad is None else [float(x) for x in mad]))
    # ---------- post-ictal recovery (events with a clean >=120 min tail)
    onsl = np.array([ons[k]["onset"] for k in ids]); st = np.array([ons[k]["start_true"] for k in ids])
    for i, k in enumerate(ids):
        e = ons[k]
        nxt = st[st > e["end"]]
        if nxt.size and nxt.min() - e["end"] < 120*60: continue
        for (a, b) in POST:
            t0 = e["end"] + a*60.0; t1 = e["end"] + b*60.0
            exp = (b-a)*6.0
            npres, nval, med = window_stats(D, t0, t1, e["session"])
            if npres == 0 or npres/exp < MIN_COVER or nval < MIN_VALID: continue
            hh = int(hod(np.array([(t0+t1)/2]))[0])
            rm = ref_h.get(hh)
            out["post"].append(dict(a=aid, hor=f"{a},{b}", hod=hh, n=nval,
                                    v=[float(x) for x in med],
                                    rm=None if rm is None else [float(x) for x in rm],
                                    sd=None if mad is None else [float(x) for x in mad]))
    print(aid, "windows", len(out["windows"]), "post", len(out["post"]), flush=True)

out["diag"] = diag
out["primary_animals"] = sorted(PRIMARY)
json.dump(out, open("windows.json", "w"))
print(diag)
