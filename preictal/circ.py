import json,sys,os,numpy as np,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from core import *
cohort=json.load(open("cohort.json"))["per_animal"]
PRIM=[a for a,v in cohort.items() if v.get("n_lead_4h",0)>=10]
rows=[]; tot_all=np.zeros(24); tot_lead=np.zeros(24); tot_rec=np.zeros(24); tot_pool=np.zeros(24)
for aid in animals():
    if cohort.get(aid,{}).get("n_events",0)==0: continue
    D=order_by_time(load(aid)); ons,_=event_onsets(D); ids,lead=lead_mask(ons,4*HOUR)
    h_all=hod(np.array([ons[k]["onset"] for k in ids]))
    h_lead=hod(np.array([ons[k]["onset"] for k in ids if lead[k]]))
    h_rec=hod(D["t"])
    pool=interictal_pool(D,ons); h_pool=hod(D["t"][pool])
    A=np.bincount(h_all,minlength=24); L=np.bincount(h_lead,minlength=24)
    R=np.bincount(h_rec,minlength=24)*10/3600.; P=np.bincount(h_pool,minlength=24)*10/3600.
    if aid in PRIM: tot_all+=A; tot_lead+=L; tot_rec+=R; tot_pool+=P
    lf=lambda v: v[6:18].sum()/max(1e-9,v.sum())
    rows.append(dict(a=aid,n_all=int(A.sum()),n_lead=int(L.sum()),
                     light_all=float(lf(A)),light_lead=float(lf(L)),light_rec=float(lf(R)),
                     rate_light_all=float(A[6:18].sum()/max(1e-9,R[6:18].sum())*100),
                     rate_dark_all=float((A.sum()-A[6:18].sum())/max(1e-9,R.sum()-R[6:18].sum())*100),
                     rate_light_lead=float(L[6:18].sum()/max(1e-9,R[6:18].sum())*100),
                     rate_dark_lead=float((L.sum()-L[6:18].sum())/max(1e-9,R.sum()-R[6:18].sum())*100),
                     primary=aid in PRIM))
for r in rows: print(f"{r['a']} all={r['n_all']:5d} lead={r['n_lead']:4d} lightfrac all={r['light_all']:.3f} lead={r['light_lead']:.3f} rec={r['light_rec']:.3f}  rate/100h all L{r['rate_light_all']:6.1f} D{r['rate_dark_all']:6.1f} | lead L{r['rate_light_lead']:5.1f} D{r['rate_dark_lead']:5.1f}")
def lf(v): return v[6:18].sum()/v.sum()
print("\nPRIMARY-cohort pooled: all onsets lightfrac %.3f  lead lightfrac %.3f  recording lightfrac %.3f  pool lightfrac %.3f"%(lf(tot_all),lf(tot_lead),lf(tot_rec),lf(tot_pool)))
print("rate/100h all: light %.1f dark %.1f ratio %.2f"%(tot_all[6:18].sum()/tot_rec[6:18].sum()*100,(tot_all.sum()-tot_all[6:18].sum())/(tot_rec.sum()-tot_rec[6:18].sum())*100,(tot_all[6:18].sum()/tot_rec[6:18].sum())/((tot_all.sum()-tot_all[6:18].sum())/(tot_rec.sum()-tot_rec[6:18].sum()))))
print("rate/100h lead: light %.1f dark %.1f ratio %.2f"%(tot_lead[6:18].sum()/tot_rec[6:18].sum()*100,(tot_lead.sum()-tot_lead[6:18].sum())/(tot_rec.sum()-tot_rec[6:18].sum())*100,(tot_lead[6:18].sum()/tot_rec[6:18].sum())/((tot_lead.sum()-tot_lead[6:18].sum())/(tot_rec.sum()-tot_rec[6:18].sum()))))
# animal-clustered ratio for lead
import scipy.stats as st
rl=[np.log(r['rate_light_lead']/r['rate_dark_lead']) for r in rows if r['primary'] and r['rate_dark_lead']>0 and r['rate_light_lead']>0]
ra=[np.log(r['rate_light_all']/r['rate_dark_all']) for r in rows if r['primary'] and r['rate_dark_all']>0 and r['rate_light_all']>0]
for nm,v in (("all",ra),("lead",rl)):
    v=np.array(v); n=len(v); m=v.mean(); se=v.std(ddof=1)/np.sqrt(n); tc=st.t.ppf(.975,n-1)
    print(f"{nm}: n={n} light/dark rate ratio {np.exp(m):.2f} CI [{np.exp(m-tc*se):.2f},{np.exp(m+tc*se):.2f}] p={st.ttest_1samp(v,0).pvalue:.3g} animals>1: {(v>0).sum()}/{n}")
json.dump(dict(rows=rows,tot_all=tot_all.tolist(),tot_lead=tot_lead.tolist(),tot_rec=tot_rec.tolist(),tot_pool=tot_pool.tolist()),open("circ.json","w"),indent=1)
