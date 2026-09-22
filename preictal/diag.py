import json, numpy as np, collections
from scipy import stats
W=json.load(open("windows.json")); PRIM=set(W["primary_animals"])
HOR=["-60,-50","-30,-25","-15,-10","-10,-5","-5,-2","-2,-1","-1,0"]
# 1. which primary animals lose the matched control
cnt=collections.Counter(); cntm=collections.Counter(); nom=collections.Counter()
for w in W["windows"]:
    if w["a"] not in PRIM: continue
    cnt[w["a"]]+=1
    if w["rm"] is not None: cntm[w["a"]]+=1
    else: nom[w["a"]]+=1
print("per-animal windows total / with matched ctrl / without:")
for a in sorted(PRIM): print(f"  {a} {cnt[a]:5d} {cntm[a]:5d} {nom[a]:5d}")
print("total windows(primary)",sum(cnt.values()),"matched",sum(cntm.values()),"unmatched-fail",sum(nom.values()))
# per-horizon usable animal counts
for h in HOR:
    per=collections.Counter(w["a"] for w in W["windows"] if w["a"] in PRIM and w["hor"]==h and w["rm"] is not None)
    print(h,"animals with >=5 matched windows:",sum(1 for v in per.values() if v>=5), sorted(a for a,v in per.items() if v<5))
# 2. hour-of-day contrast: pre-ictal windows vs pool
ph=collections.Counter(w["hod"] for w in W["windows"] if w["a"] in PRIM and w["hor"]=="-30,-25")
print("\npre-ictal window hour histogram (-30..-25):",[ph[h] for h in range(24)])
light=sum(ph[h] for h in range(6,18)); dark=sum(ph.values())-light
print("light 06-18:",light,"dark:",dark,"frac light", round(light/max(1,light+dark),3))
