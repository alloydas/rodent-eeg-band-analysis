import json,sys,os,numpy as np,collections
from scipy import stats
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from core import *
HOR=[(-60,-50),(-30,-25),(-15,-10),(-10,-5),(-5,-2),(-2,-1),(-1,0)]
cohort=json.load(open("cohort.json"))["per_animal"]
MINLEAD=10; MIN_EV=5; MIN_CTRL=30
def window_med(D,t0,t1,sess):
    t=D["t"]; lo=np.searchsorted(t,t0,"left"); hi=np.searchsorted(t,t1,"left")
    if hi<=lo: return None,0,0
    sl=slice(lo,hi); m=(D["sess"][sl]==sess); npres=int(m.sum())
    v=m&D["valid"][sl]
    if v.sum()<3: return None,npres,int(v.sum())
    return np.median(D["F"][sl][v],axis=0),npres,int(v.sum())
OUT={}
CACHE={}
for aid in animals():
    if cohort.get(aid,{}).get("n_events",0)==0: continue
    CACHE[aid]=order_by_time(load(aid))
for LH in (1.0,2.0,4.0):
    eff=collections.defaultdict(list); nlead={}
    for aid,D in CACHE.items():
        ons,_=event_onsets(D); ids,lead=lead_mask(ons,LH*HOUR)
        nlead[aid]=int(sum(lead.values()))
        if nlead[aid]<MINLEAD: continue
        pool=interictal_pool(D,ons); ph=hod(D["t"][pool]); PF=D["F"][pool]
        if pool.sum()<MIN_CTRL: continue
        ref_all=np.median(PF,axis=0); sd=np.median(np.abs(PF-ref_all),axis=0)*1.4826
        sd=np.where(sd<=0,np.nan,sd)
        ref_h={h:np.median(PF[ph==h],axis=0) for h in range(24) if (ph==h).sum()>=MIN_CTRL}
        for k in ids:
            if not lead[k]: continue
            e=ons[k]
            for (a,b) in HOR:
                t0=e["onset"]+a*60; t1=e["onset"]+b*60
                med,npres,nval=window_med(D,t0,t1,e["session"])
                if med is None or npres/((b-a)*6.)<0.5: continue
                hh=int(hod(np.array([(t0+t1)/2]))[0])
                if hh not in ref_h: continue
                eff[(aid,f"{a},{b}")].append((med-ref_h[hh])/sd)
    rows=[]
    for (a,b) in HOR:
        h=f"{a},{b}"
        per=[np.median(np.array(v),axis=0) for (aa,hh),v in eff.items() if hh==h and len(v)>=MIN_EV]
        if len(per)<2: continue
        A=np.array(per); n=len(A)
        for i,e in enumerate(ENAMES):
            v=A[:,i]; v=v[np.isfinite(v)]; nn=len(v)
            m=v.mean(); s=v.std(ddof=1); se=s/np.sqrt(nn); tc=stats.t.ppf(.975,nn-1)
            mde=(stats.t.ppf(.975,nn-1)+stats.t.ppf(.80,nn-1))*se
            rows.append(dict(lead_h=LH,horizon=h,endpoint=e,n_animals=nn,
                             n_events=int(sum(len(v2) for (aa,hh),v2 in eff.items() if hh==h and len(v2)>=MIN_EV)),
                             mean=float(m),lo=float(m-tc*se),hi=float(m+tc*se),
                             p=float(stats.ttest_1samp(v,0).pvalue),n_pos=int((v>0).sum()),
                             mde80=float(mde)))
    OUT[str(LH)]=dict(rows=rows,n_lead=nlead,n_animals_ge10=sum(1 for v in nlead.values() if v>=MINLEAD))
    print(f"--- lead={LH}h  animals>=10 lead: {OUT[str(LH)]['n_animals_ge10']}  total lead events {sum(nlead.values())}")
    for r in rows:
        if r["endpoint"] in ("log10_total_power","tilt_log10_gamma_over_delta"):
            print(f"   {r['horizon']:>10} {r['endpoint']:>28} nA={r['n_animals']:2d} nEv={r['n_events']:4d} d={r['mean']:+.3f} [{r['lo']:+.3f},{r['hi']:+.3f}] p={r['p']:.3g} MDE80=+-{r['mde80']:.3f}")
json.dump(OUT,open("leadsens.json","w"),indent=1)
