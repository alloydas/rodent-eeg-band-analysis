"""What each trap would have produced if it had NOT been controlled.
Pre-ictal band = epochs in [-30,-5] min before onset. Five nested analyses."""
import json,sys,os,numpy as np,collections
from scipy import stats
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from core import *
cohort=json.load(open("cohort.json"))["per_animal"]
PRIM=[a for a,v in cohort.items() if v.get("n_lead_4h",0)>=10]
EI={e:i for i,e in enumerate(ENAMES)}
ALLPRE=[]; ALLCTL=[]; peranimal={}
for aid in animals():
    if cohort.get(aid,{}).get("n_events",0)==0: continue
    D=order_by_time(load(aid)); ons,_=event_onsets(D); ids,lead=lead_mask(ons,4*HOUR)
    t=D["t"]; pre_all=np.zeros(len(t),bool); pre_lead=np.zeros(len(t),bool)
    for k in ids:
        e=ons[k]; lo=np.searchsorted(t,e["onset"]-30*60,"left"); hi=np.searchsorted(t,e["onset"]-5*60,"left")
        m=slice(lo,hi)
        sel=(D["sess"][m]==e["session"])
        idx=np.arange(lo,hi)[sel]
        pre_all[idx]=True
        if lead[k]: pre_lead[idx]=True
    ict_near=np.zeros(len(t),bool)
    for e in ons.values():
        lo=np.searchsorted(t,e["start_true"],"left"); hi=np.searchsorted(t,e["end"],"right"); ict_near[lo:hi]=True
    pre_all&=D["valid"]&~ict_near; pre_lead&=D["valid"]&~ict_near
    ctl_naive=D["valid"]&(D["ovf"]==0)&~pre_all&~ict_near          # ALL non-ictal, no 4h margin
    ctl_strict=interictal_pool(D,ons)
    peranimal[aid]=dict(F=D["F"],pre_all=pre_all,pre_lead=pre_lead,ctl_naive=ctl_naive,
                        ctl_strict=ctl_strict,hod=hod(t))
    if aid in PRIM:
        ALLPRE.append(D["F"][pre_all]); ALLCTL.append(D["F"][ctl_naive])
P=np.vstack(ALLPRE); C=np.vstack(ALLCTL)
print("A) NAIVE epoch-pooled, all events, no 4h control margin, no hour match, no animal clustering")
print(f"   n_preictal_epochs={len(P):,}  n_control_epochs={len(C):,}")
for e in ENAMES:
    i=EI[e]; a=P[:,i]; b=C[:,i]; a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
    sd=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    tt=stats.ttest_ind(a,b,equal_var=False)
    print(f"   {e:>28} d={(a.mean()-b.mean())/sd:+.3f}  t={tt.statistic:+9.1f}  p={tt.pvalue:.3g}")

def clus(vals):
    v=np.array([x for x in vals if np.isfinite(x)]); n=len(v)
    m=v.mean(); se=v.std(ddof=1)/np.sqrt(n); tc=stats.t.ppf(.975,n-1)
    return n,m,m-tc*se,m+tc*se,stats.ttest_1samp(v,0).pvalue,int((v>0).sum())

def run(tag,prekey,ctlkey,match):
    out={}
    print(f"\n{tag}")
    for e in ENAMES:
        i=EI[e]; per=[]
        for aid in PRIM:
            d=peranimal[aid]; F=d["F"][:,i]; pre=d[prekey]; ctl=d[ctlkey]
            if pre.sum()<30 or ctl.sum()<30: continue
            ref=np.median(F[ctl]); sd=np.median(np.abs(F[ctl]-ref))*1.4826
            if not np.isfinite(sd) or sd<=0: continue
            if match:
                num=[]; H=d["hod"]
                for h in range(24):
                    cm=ctl&(H==h); pm=pre&(H==h)
                    if cm.sum()<30 or pm.sum()<3: continue
                    num.append(np.median(F[pm])-np.median(F[cm]))
                if not num: continue
                per.append(float(np.median(num))/sd)
            else:
                per.append(float(np.median(F[pre])-ref)/sd)
        n,m,lo,hi,p,npos=clus(per)
        out[e]=dict(n=n,mean=float(m),lo=float(lo),hi=float(hi),p=float(p),n_pos=npos)
        print(f"   {e:>28} nA={n:2d} d={m:+.3f} [{lo:+.3f},{hi:+.3f}] p={p:.3g} sign {npos}/{n}")
    return out

R={}
R["A_naive_pooled"]={e:None for e in ENAMES}
R["B_cluster_allev_naivectl"]=run("B) + animal clustering only (all events, naive controls, no hour match)","pre_all","ctl_naive",False)
R["C_cluster_allev_strictctl"]=run("C) + 4h-clean controls (all events, no hour match)","pre_all","ctl_strict",False)
R["D_cluster_lead_strictctl"]=run("D) + LEAD-seizure restriction (no hour match)","pre_lead","ctl_strict",False)
R["E_cluster_lead_strict_matched"]=run("E) + hour-of-day matching  [FULL CONTROL]","pre_lead","ctl_strict",True)
json.dump(R,open("naive.json","w"),indent=1)
