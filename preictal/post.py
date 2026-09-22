import json,numpy as np,collections
from scipy import stats
W=json.load(open("windows.json")); PRIM=set(W["primary_animals"])
from core import ENAMES
POST=["0,1","1,2","2,5","5,10","10,20","20,40","40,60","60,120"]
per=collections.defaultdict(list)
for w in W["post"]:
    if w["a"] not in PRIM or w["rm"] is None or w["sd"] is None: continue
    per[(w["a"],w["hor"])].append((np.array(w["v"])-np.array(w["rm"]))/np.array(w["sd"]))
rows=[]
print(f"{'post win (min)':>14} {'nA':>3} {'nEv':>5} "+" ".join(f"{e[:14]:>14}" for e in ENAMES))
for h in POST:
    vals={e:[] for e in ENAMES}; nev=0; na=0
    for (a,hh),L in per.items():
        if hh!=h or len(L)<5: continue
        na+=1; nev+=len(L); m=np.median(np.array(L),axis=0)
        for i,e in enumerate(ENAMES): vals[e].append(m[i])
    r=dict(horizon=h,n_animals=na,n_events=nev)
    for e in ENAMES:
        v=np.array(vals[e]); n=len(v)
        if n<2: r[e]=(float('nan'),)*3; continue
        mm=v.mean(); se=v.std(ddof=1)/np.sqrt(n); tc=stats.t.ppf(.975,n-1)
        r[e]=(float(mm),float(mm-tc*se),float(mm+tc*se))
    rows.append(r)
    print(f"{h:>14} {na:>3} {nev:>5} "+" ".join(f"{r[e][0]:>14.3f}" for e in ENAMES))
json.dump(rows,open("post.json","w"),indent=1)
