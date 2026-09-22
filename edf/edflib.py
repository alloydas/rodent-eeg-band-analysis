"""Minimal dependency-free EDF / EDF+ reader (numpy only).

Written because this box has numpy+scipy but no mne / edfio.
Returns physical values exactly as the EDF spec defines them:
    phys = (digital - dig_min) * (phys_max - phys_min)/(dig_max - dig_min) + phys_min
"""
import datetime as _dt
import numpy as np


class EdfHeader:
    __slots__ = ("path", "version", "patient", "recording", "start", "n_header_bytes",
                 "reserved", "n_records", "record_dur", "n_signals", "labels",
                 "transducers", "units", "phys_min", "phys_max", "dig_min", "dig_max",
                 "prefilt", "n_samp_rec", "sfreqs", "duration_s")

    def __repr__(self):
        return (f"<EDF {self.path.split('/')[-1]} start={self.start} "
                f"n_sig={self.n_signals} dur={self.duration_s:.3f}s "
                f"labels={self.labels}>")


def _s(b):
    return b.decode("latin-1").strip()


def read_header(path):
    h = EdfHeader()
    h.path = path
    with open(path, "rb") as f:
        raw = f.read(256)
        if len(raw) < 256:
            raise ValueError(f"{path}: truncated EDF header")
        h.version = _s(raw[0:8])
        h.patient = _s(raw[8:88])
        h.recording = _s(raw[88:168])
        d, t = _s(raw[168:176]), _s(raw[176:184])
        h.n_header_bytes = int(_s(raw[184:192]))
        h.reserved = _s(raw[192:236])
        h.n_records = int(_s(raw[236:244]))
        h.record_dur = float(_s(raw[244:252]))
        h.n_signals = int(_s(raw[252:256]))

        # EDF clipping rule: 2-digit year 85..99 -> 19xx, 00..84 -> 20xx.
        # EDF+ puts the unclipped date in the recording field ("Startdate dd-MMM-yyyy").
        dd, mm, yy = (int(x) for x in d.split("."))
        year = 1900 + yy if yy >= 85 else 2000 + yy
        if h.recording.startswith("Startdate"):
            parts = h.recording.split()
            if len(parts) > 1 and parts[1] != "X":
                try:
                    year = _dt.datetime.strptime(parts[1], "%d-%b-%Y").year
                except ValueError:
                    pass
        hh, mi, ss = (int(x) for x in t.split("."))
        h.start = _dt.datetime(year, mm, dd, hh, mi, ss)

        ns = h.n_signals
        buf = f.read(h.n_header_bytes - 256)

    def field(off, width):
        base = off * ns
        return [_s(buf[base + i * width: base + (i + 1) * width]) for i in range(ns)]

    h.labels = field(0, 16)
    o = 16 * ns
    h.transducers = [_s(buf[o + i * 80: o + (i + 1) * 80]) for i in range(ns)]
    o += 80 * ns
    h.units = [_s(buf[o + i * 8: o + (i + 1) * 8]) for i in range(ns)]
    o += 8 * ns
    def nums(width, cast):
        nonlocal o
        vals = [cast(_s(buf[o + i * width: o + (i + 1) * width])) for i in range(ns)]
        o += width * ns
        return vals
    h.phys_min = nums(8, float)
    h.phys_max = nums(8, float)
    h.dig_min = nums(8, float)
    h.dig_max = nums(8, float)
    h.prefilt = [_s(buf[o + i * 80: o + (i + 1) * 80]) for i in range(ns)]
    o += 80 * ns
    h.n_samp_rec = nums(8, int)
    h.sfreqs = [n / h.record_dur for n in h.n_samp_rec]
    h.duration_s = h.n_records * h.record_dur
    return h


def read_signal(path, sig, start_s=0.0, dur_s=None, header=None):
    """Read one signal's physical values over [start_s, start_s+dur_s).

    Reads only the data records it needs. Returns (values_float64, true_start_s).
    """
    h = header or read_header(path)
    n_per_rec = h.n_samp_rec[sig]
    rec_samples = sum(h.n_samp_rec)
    rec_bytes = rec_samples * 2
    off_in_rec = sum(h.n_samp_rec[:sig]) * 2
    sf = h.sfreqs[sig]

    if dur_s is None:
        dur_s = h.duration_s - start_s
    r0 = max(0, int(np.floor(start_s / h.record_dur)))
    r1 = min(h.n_records, int(np.ceil((start_s + dur_s) / h.record_dur)))
    if r1 <= r0:
        return np.empty(0), r0 * h.record_dur

    out = np.empty((r1 - r0) * n_per_rec, dtype=np.int16)
    with open(path, "rb") as f:
        for k, r in enumerate(range(r0, r1)):
            f.seek(h.n_header_bytes + r * rec_bytes + off_in_rec)
            chunk = f.read(n_per_rec * 2)
            out[k * n_per_rec:(k + 1) * n_per_rec] = np.frombuffer(chunk, dtype="<i2")

    dmin, dmax = h.dig_min[sig], h.dig_max[sig]
    pmin, pmax = h.phys_min[sig], h.phys_max[sig]
    gain = (pmax - pmin) / (dmax - dmin)
    phys = (out.astype(np.float64) - dmin) * gain + pmin

    block_start_s = r0 * h.record_dur
    i0 = int(round((start_s - block_start_s) * sf))
    i1 = i0 + int(round(dur_s * sf))
    i0 = max(0, i0)
    i1 = min(phys.size, i1)
    return phys[i0:i1], block_start_s + i0 / sf
