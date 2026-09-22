import json, os, sys, numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import ENAMES, ENDPOINTS

W = json.load(open("windows.json"))
PRIM = set(W["primary_animals"])
HOR = ["-60,-50","-30,-25","-15,-10","-10,-5","-5,-2","-2,-1","-1,0"]
HLAB = {"-60,-50":"-60..-50","-30,-25":"-30..-25","-15,-10":"-15..-10","-10,-5":"-10..-5",
        "-5,-2":"-5..-2","-2,-1":"-2..-1","-1,0":"-1..0 (SUSPECT)"}
MIN_EV = 5

def clustered(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    n = len(v)
    if n < 2: return dict(n=n, mean=float(v[0]) if n else float('nan'))
    m = v.mean(); sd = v.std(ddof=1); se = sd/np.sqrt(n); tc = stats.t.ppf(0.975, n-1)
    r = dict(n=n, mean=float(m), median=float(np.median(v)), lo=float(m-tc*se), hi=float(m+tc*se),
             n_pos=int((v>0).sum()), n_neg=int((v<0).sum()),
             p=float(stats.ttest_1samp(v, 0).pvalue))
    try: r["p_wil"]=float(stats.wilcoxon(v).pvalue)
    except Exception: r["p_wil"]=float('nan')
    return r

def bh(pv):
    p = np.asarray(pv, float); o = np.argsort(p); m = len(p)
    q = np.empty(m); run = 1.0
    for i in range(m-1, -1, -1):
        run = min(run, p[o[i]]*m/(i+1)); q[o[i]] = run
    return q

# per-animal, per-horizon, per-endpoint effect (median over events of (window - ctrl)/sd)
def build(mode):   # mode 'rm' matched, 'ru' unmatched
    per = {}
    for w in W["windows"]:
        if w["a"] not in PRIM: continue
        ref = w.get(mode)
        if ref is None or w["sd"] is None: continue
        z = (np.array(w["v"]) - np.array(ref)) / np.array(w["sd"])
        per.setdefault((w["a"], w["hor"]), []).append(z)
    eff = {}
    for (a, h), L in per.items():
        A = np.array(L)
        if len(A) < MIN_EV: continue
        eff[(a, h)] = (np.median(A, axis=0), len(A))
    return eff

res = {}
for mode in ("rm", "ru"):
    eff = build(mode)
    rows = []
    for hi_, h in enumerate(HOR):
        for ei, en in enumerate(ENAMES):
            vals = []; ani = []
            for (a, hh), (med, n) in eff.items():
                if hh != h: continue
                vals.append(med[ei]); ani.append(a)
            c = clustered(vals)
            c.update(horizon=h, hlabel=HLAB[h], endpoint=en,
                     family=dict(ENDPOINTS)[en],
                     n_events_total=int(sum(n for (a,hh),(m_,n) in eff.items() if hh==h)))
            rows.append(c)
    # BH within family, across horizons x endpoints
    for fam in ("primary","secondary"):
        idx = [i for i,r in enumerate(rows) if r["family"]==fam and np.isfinite(r.get("p",np.nan))]
        q = bh([rows[i]["p"] for i in idx])
        for i,qq in zip(idx,q): rows[i]["q_bh"]=float(qq)
    res[mode] = rows

json.dump(res, open("effects.json","w"), indent=1)

lab = {"rm":"HOUR-MATCHED control","ru":"UNMATCHED control"}
for mode in ("rm","ru"):
    print("="*100); print(lab[mode])
    print(f"{'horizon':>16} {'endpoint':>28} {'nA':>3} {'nEv':>5} {'d(mean)':>8} {'95% CI':>18} {'+/-':>7} {'p':>9} {'q':>8}")
    for r in res[mode]:
        if "lo" not in r: print(r["hlabel"], r["endpoint"], "n<2"); continue
        print(f"{r['hlabel']:>16} {r['endpoint']:>28} {r['n']:>3} {r['n_events_total']:>5} "
              f"{r['mean']:>8.3f} [{r['lo']:>7.3f},{r['hi']:>7.3f}] {r['n_pos']:>3}/{r['n_neg']:<3} "
              f"{r['p']:>9.2e} {r.get('q_bh',float('nan')):>8.3f}")
