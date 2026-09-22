"""(a) 1-min descriptive time course -10..0 ; (b) window-level in-sample AUROC vs
hour-matched interictal windows (characterisation, NOT prospective prediction)."""
import json,sys,os,numpy as np,collections
from scipy import stats
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from core import *
cohort=json.load(open("cohort.json"))["per_animal"]
PRIM=[a for a,v in cohort.items() if v.get("n_lead_4h",0)>=10]
MIN_CTRL=30
def clus(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]; n=len(v)
    m=v.mean(); se=v.std(ddof=1)/np.sqrt(n); tc=stats.t.ppf(.975,n-1)
    return dict(n=n,mean=float(m),lo=float(m-tc*se),hi=float(m+tc*se),
                p=float(stats.ttest_1samp(v,0).pvalue),n_pos=int((v>0).sum()))
def auroc(x,y):
    x=np.asarray(x,float); y=np.asarray(y,bool); n1=y.sum(); n0=(~y).sum()
    if n1==0 or n0==0: return np.nan
    r=stats.rankdata(x); return float((r[y].sum()-n1*(n1+1)/2)/(n1*n0))
FINE=[(-m-1,-m) for m in range(10)][::-1]           # -10..-9 ... -1..0
AU_H=[(-30,-25),(-5,-2)]
fine=collections.defaultdict(list); au=collections.defaultdict(list); aun=collections.defaultdict(list)
for aid in animals():
    if aid not in PRIM: continue
    D=order_by_time(load(aid)); ons,_=event_onsets(D); ids,lead=lead_mask(ons,4*HOUR)
    pool=interictal_pool(D,ons); ph=hod(D["t"][pool]); PF=D["F"][pool]; PT=D["t"][pool]
    if pool.sum()<MIN_CTRL: continue
    ref_all=np.median(PF,axis=0); sd=np.median(np.abs(PF-ref_all),axis=0)*1.4826
    sd=np.where(sd<=0,np.nan,sd)
    ref_h={h:np.median(PF[ph==h],axis=0) for h in range(24) if (ph==h).sum()>=MIN_CTRL}
    t=D["t"]
    def wmed(t0,t1,sess):
        lo=np.searchsorted(t,t0,"left"); hi=np.searchsorted(t,t1,"left")
        if hi<=lo: return None
        sl=slice(lo,hi); m=(D["sess"][sl]==sess); v=m&D["valid"][sl]
        if v.sum()<3 or m.sum()/((t1-t0)/10.)<0.5: return None
        return np.median(D["F"][sl][v],axis=0)
    for k in ids:
        if not lead[k]: continue
        e=ons[k]
        for (a,b) in FINE:
            med=wmed(e["onset"]+a*60,e["onset"]+b*60,e["session"])
            if med is None: continue
            hh=int(hod(np.array([e["onset"]+(a+b)*30]))[0])
            if hh not in ref_h: continue
            fine[(aid,f"{a},{b}")].append((med-ref_h[hh])/sd)
    # ---- control windows: consecutive blocks inside the pool, same length, same hour
    for (a,b) in AU_H:
        L=int((b-a)*6)                                  # epochs per window
        preW=[];preH=[]
        for k in ids:
            if not lead[k]: continue
            e=ons[k]; med=wmed(e["onset"]+a*60,e["onset"]+b*60,e["session"])
            if med is None: continue
            preW.append(med); preH.append(int(hod(np.array([e["onset"]+(a+b)*30]))[0]))
        # build control blocks from contiguous pool runs
        ctlW=[];ctlH=[]
        brk=np.flatnonzero(np.diff(PT)>15.0)+1
        for s,eend in zip(np.r_[0,brk],np.r_[brk,len(PT)]):
            for j in range(s,eend-L+1,L):
                blk=PF[j:j+L]
                if not np.isfinite(blk).all(): continue
                ctlW.append(np.median(blk,axis=0)); ctlH.append(int(hod(np.array([PT[j+L//2]]))[0]))
        if len(preW)<10 or len(ctlW)<20: continue
        preW=np.array(preW); ctlW=np.array(ctlW); preH=np.array(preH); ctlH=np.array(ctlH)
        # hour-adjust both by the animal's per-hour control median
        def adj(Wm,Hh):
            out=np.full_like(Wm,np.nan)
            for h in np.unique(Hh):
                if h not in ref_h: continue
                out[Hh==h]=(Wm[Hh==h]-ref_h[h])/sd
            return out
        A1=adj(preW,preH); A0=adj(ctlW,ctlH)
        ok1=np.isfinite(A1).all(1); ok0=np.isfinite(A0).all(1)
        A1=A1[ok1]; A0=A0[ok0]
        if len(A1)<10 or len(A0)<20: continue
        for i,en in enumerate(ENAMES):
            au[(f"{a},{b}",en)].append(auroc(np.r_[A0[:,i],A1[:,i]],
                                             np.r_[np.zeros(len(A0),bool),np.ones(len(A1),bool)]))
        aun[f"{a},{b}"].append((aid,len(A1),len(A0)))
out={"fine":[], "auroc":[]}
print("1-min descriptive time course (hour-matched, robust-z, mean over animals):")
print(f"{'win':>9} {'nA':>3} "+" ".join(f"{e[:13]:>13}" for e in ENAMES))
for (a,b) in FINE:
    h=f"{a},{b}"
    per=[np.median(np.array(v),axis=0) for (aa,hh),v in fine.items() if hh==h and len(v)>=5]
    if len(per)<2: continue
    A=np.array(per); r=dict(horizon=h,n_animals=len(A),
                            n_events=int(sum(len(v) for (aa,hh),v in fine.items() if hh==h and len(v)>=5)))
    for i,e in enumerate(ENAMES): r[e]=clus(A[:,i])
    out["fine"].append(r)
    print(f"{h:>9} {len(A):>3} "+" ".join(f"{r[e]['mean']:>13.3f}" for e in ENAMES))
print("\nWindow-level in-sample AUROC (pre-ictal window vs hour-matched interictal window), per animal then clustered:")
for h in ["-30,-25","-5,-2"]:
    print(f"  horizon {h}: animals={len(aun[h])} pre-windows={sum(x[1] for x in aun[h])} ctrl-windows={sum(x[2] for x in aun[h])}")
    for e in ENAMES:
        v=np.array(au[(h,e)]); c=clus(v-0.5)
        out["auroc"].append(dict(horizon=h,endpoint=e,n_animals=c["n"],mean_auroc=0.5+c["mean"],
                                 lo=0.5+c["lo"],hi=0.5+c["hi"],p=c["p"],
                                 n_pre=int(sum(x[1] for x in aun[h])),n_ctl=int(sum(x[2] for x in aun[h]))))
        print(f"      {e:>28} AUROC={0.5+c['mean']:.3f} [{0.5+c['lo']:.3f},{0.5+c['hi']:.3f}] p={c['p']:.3g}")
json.dump(out,open("fine.json","w"),indent=1)
