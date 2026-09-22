"""Pre-ictal cohort construction + endpoint extraction.

Source: /work/mech-ai-scratch/alloy/EEG/cache/epochs/*.npz   (READ ONLY)
Contract already enforced there: channel_index==1 (dedupe on recording_id,epoch_index),
RoomD/RN245 dropped, .interrupted quarantined, pandas csv reader (never awk -F,),
event stage from SOLE-id rows, event bounds from every touching epoch.
Here: numpy/scipy only.
"""
import json, glob, os
import numpy as np

CACHE = "/work/mech-ai-scratch/alloy/EEG/cache/epochs"
MAN = json.load(open(os.path.join(CACHE, "manifest.json")))
FEATURES = MAN["features"]
FI = {c: i for i, c in enumerate(FEATURES)}
BANDS = ["delta", "theta", "alpha", "beta", "gamma"]

# endpoint definitions, all in log10 units
ENDPOINTS = [
    ("log10_total_power", "primary"),
    ("tilt_log10_gamma_over_delta", "primary"),
    ("log10_delta", "secondary"), ("log10_theta", "secondary"),
    ("log10_alpha", "secondary"), ("log10_beta", "secondary"),
    ("log10_gamma", "secondary"), ("log10_line_length", "secondary"),
]
ENAMES = [e for e, _ in ENDPOINTS]

HOUR = 3600.0


def animals():
    return sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(CACHE, "*.npz"))
                  if "_" not in os.path.basename(p)[:-4])


def load(aid):
    d = np.load(os.path.join(CACHE, f"{aid}.npz"), allow_pickle=True)
    X = d["X"]; t = d["t"].astype(np.float64); dur = d["dur"].astype(np.float64)
    valid = d["valid"]; ovf = d["ovf"].astype(np.float64); sess = d["session"]
    # endpoint matrix
    tp = np.log10(np.maximum(X[:, FI["total_power_0p5_80_hz_uV2"]].astype(np.float64), 1e-12))
    lg = {b: X[:, FI[f"{b}_log10_power_uV2"]].astype(np.float64) for b in BANDS}
    ll = np.log10(np.maximum(X[:, FI["line_length_uV_per_s"]].astype(np.float64), 1e-12))
    F = np.column_stack([tp, lg["gamma"] - lg["delta"],
                         lg["delta"], lg["theta"], lg["alpha"], lg["beta"], lg["gamma"], ll])
    F[~valid] = np.nan
    ev = json.load(open(os.path.join(CACHE, f"{aid}_events.json")))
    return dict(t=t, dur=dur, valid=valid, ovf=ovf, sess=sess, F=F, ev=ev, aid=aid)


def order_by_time(D):
    o = np.argsort(D["t"], kind="stable")
    for k in ("t", "dur", "valid", "ovf", "sess"):
        D[k] = D[k][o]
    D["F"] = D["F"][o]
    return D


def event_onsets(D):
    """Onset = start of the EARLIEST epoch with annotation_overlap_fraction > 0
    that overlaps the event, in the event's own session."""
    t, dur, ovf, sess = D["t"], D["dur"], D["ovf"], D["sess"]
    touch = ovf > 0
    out = {}
    miss = 0
    # index epochs by session for the touching set
    ti = t[touch]; di = dur[touch]; si = sess[touch]; vi = D["valid"][touch]
    order = np.argsort(ti, kind="stable")
    ti, di, si, vi = ti[order], di[order], si[order], vi[order]
    for eid, e in D["ev"].items():
        lo = np.searchsorted(ti, e["start"] - 60.0, "left")
        hi = np.searchsorted(ti, e["end"] + 60.0, "right")
        if hi <= lo:
            miss += 1; continue
        m = (si[lo:hi] == e["session"]) & (ti[lo:hi] < e["end"]) & (ti[lo:hi] + di[lo:hi] > e["start"])
        idx = np.flatnonzero(m)
        if idx.size == 0:
            miss += 1; continue
        j = lo + idx[0]
        out[eid] = dict(onset=float(ti[j]), end=float(e["end"]), start_true=float(e["start"]),
                        session=e["session"], label=e.get("label", "<unresolved>"),
                        n_touch=int(idx.size), onset_epoch_valid=bool(vi[j]))
    return out, miss


def lead_mask(ons, lookback_s):
    """event is LEAD if no OTHER event overlaps [onset-lookback, onset)."""
    ids = sorted(ons, key=lambda k: ons[k]["onset"])
    starts = np.array([ons[k]["start_true"] for k in ids])
    ends = np.array([ons[k]["end"] for k in ids])
    onsets = np.array([ons[k]["onset"] for k in ids])
    lead = {}
    for i, k in enumerate(ids):
        w0 = onsets[i] - lookback_s
        other = np.ones(len(ids), bool); other[i] = False
        hit = other & (ends > w0) & (starts < onsets[i])
        lead[k] = not hit.any()
    return ids, lead


def interictal_pool(D, ons, margin_s=4 * HOUR):
    """valid, ovf==0, >= margin from every event interval (both sides)."""
    t, dur = D["t"], D["dur"]
    near = np.zeros(len(t), bool)
    for e in ons.values():
        lo = np.searchsorted(t, e["start_true"] - margin_s, "left")
        hi = np.searchsorted(t, e["end"] + margin_s, "right")
        near[lo:hi] = True
    pool = D["valid"] & (D["ovf"] == 0) & (~near)
    return pool


def hod(t):
    return np.floor((t % 86400.0) / 3600.0).astype(int)
