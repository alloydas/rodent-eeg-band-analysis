import os
"""Verify that raw_epoch_features/*.csv powerbands were computed from the same EDFs
our clips were cut from.

Method
------
Every clip EDF under data_full/Data_<ANIMAL>_cropped/<day>/<clip>/eeg.edf keeps the
SOURCE recording's start time in its EDF header (mne's raw.crop leaves meas_date
alone), and info.txt records the clip's wall-clock span. So a CSV epoch at
`epoch_start_s` seconds into the source maps to `epoch_start_s - (clip_start -
source_start)` seconds into the clip. Epochs fully inside the clip can be recomputed
sample-for-sample and compared column by column.

The clip EDFs were requantized to 16 bit on write, so a perfect match floors at
~1e-5 relative, not 0.

Discovered feature recipe (reproduces all six power columns to 2.7e-5):
    Welch PSD, nperseg=4000 (4 s), Hann, 50% overlap, mean-detrend, density scaling,
    band power = trapezoid over [lo, hi] inclusive.
    Bands: delta 0.5-4, theta 4-8, alpha 8-13, beta 13-30, gamma 30-80 Hz.
"""
import csv, datetime as dt, glob, json, os, sys
import numpy as np
from scipy import signal as sg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import edflib

import os
_DEFAULT_REPO = os.environ.get('EEG_ROOT', '/work/mech-ai-scratch/alloy/EEG')
REPO = os.environ.get('EEG_ROOT', _DEFAULT_ROOT)  # override for another machine
FEAT = os.path.join(REPO, 'raw_epoch_features')
BANDS = dict(delta=(0.5, 4.0), theta=(4.0, 8.0), alpha=(8.0, 13.0),
             beta=(13.0, 30.0), gamma=(30.0, 80.0))
csv.field_size_limit(10 ** 8)


def room_of(animal):
    for room in ('RoomC', 'RoomD'):
        for pre in ('', 'DONE_'):
            p = os.path.join(FEAT, room, pre + animal)
            if os.path.isdir(p):
                yield room, p


def csv_for(animal, day):
    """All export CSVs for an animal/day, including '(2)' continuations."""
    out = []
    for _room, d in room_of(animal):
        for c in sorted(glob.glob(os.path.join(d, '*.csv'))):
            stem = os.path.basename(c)[:-4]
            if stem == day or stem.startswith(day + '('):
                out.append(c)
    return out


def parse_info(clip_dir):
    m = {}
    with open(os.path.join(clip_dir, 'info.txt')) as f:
        for line in f:
            if ':' in line:
                k, v = line.split(':', 1)
                m[k.strip()] = v.strip()
    def ts(k):
        v = m.get(k)
        if not v:
            return None
        for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S'):
            try:
                return dt.datetime.strptime(v, fmt)
            except ValueError:
                pass
        return None
    return dict(raw=m, clip_start=ts('Clip start'), clip_end=ts('Clip end'),
                sz_start=ts('Seizure start'), sz_end=ts('Seizure end'),
                edf=m.get('EDF file'), label=m.get('Label'),
                sfreq=float(m['EEG sfreq (Hz)']) if m.get('EEG sfreq (Hz)') else None)


def biopotential_index(header):
    """Index of the channel the export used: the RAW biopotential, not the FIR-HP copy."""
    cands = [i for i, l in enumerate(header.labels)
             if l.upper().startswith(('ECG', 'EEG')) and 'FIR' not in l.upper()]
    return cands[0] if cands else None


def features(x, fs):
    n = x.size
    m = float(x.mean()); v = float(x.var(ddof=0)); sd = float(np.sqrt(v))
    d1 = np.diff(x)
    z = (x - m) / sd if sd > 0 else np.zeros_like(x)
    mob = float(np.std(d1) / sd) if sd > 0 else 0.0
    d2 = np.diff(d1)
    mob2 = float(np.std(d2) / np.std(d1)) if np.std(d1) > 0 else 0.0
    f, P = sg.welch(x, fs=fs, nperseg=min(4000, n), noverlap=min(4000, n) // 2,
                    window='hann', detrend='constant', scaling='density')

    def bp(lo, hi):
        sel = (f >= lo) & (f <= hi)
        return float(np.trapezoid(P[sel], f[sel]))

    tot = bp(0.5, 80.0)
    out = dict(mean_uV=m, std_uV=sd, variance_uV2=v, rms_uV=float(np.sqrt((x ** 2).mean())),
               min_uV=float(x.min()), max_uV=float(x.max()), peak_to_peak_uV=float(np.ptp(x)),
               skewness=float((z ** 3).mean()), excess_kurtosis=float((z ** 4).mean() - 3.0),
               line_length_uV_per_s=float(np.abs(d1).sum() / (n / fs)),
               hjorth_activity_uV2=v, hjorth_mobility_per_sample=mob,
               hjorth_complexity=(mob2 / mob if mob > 0 else 0.0),
               total_power_0p5_80_hz_uV2=tot)
    for b, (lo, hi) in BANDS.items():
        p = bp(lo, hi)
        out[f'{b}_power_uV2'] = p
        out[f'{b}_relative_power'] = p / tot if tot > 0 else 0.0
        out[f'{b}_log10_power_uV2'] = float(np.log10(p)) if p > 0 else float('nan')
    return out


EXACT = ['mean_uV', 'std_uV', 'variance_uV2', 'rms_uV', 'min_uV', 'max_uV',
         'peak_to_peak_uV', 'skewness', 'excess_kurtosis', 'line_length_uV_per_s',
         'hjorth_activity_uV2', 'hjorth_mobility_per_sample', 'hjorth_complexity']
SPECTRAL = (['total_power_0p5_80_hz_uV2'] +
            [f'{b}_{s}' for b in BANDS for s in ('power_uV2', 'relative_power', 'log10_power_uV2')])


def check_clip(clip_dir, animal, day, max_epochs=None):
    """Recompute every export epoch fully contained in this clip. Returns a dict."""
    res = dict(clip=os.path.relpath(clip_dir, REPO), animal=animal, day=day,
               n_epochs=0, status='ok', errors={}, notes=[])
    try:
        meta = parse_info(clip_dir)
    except Exception as e:
        res['status'] = f'no_info:{type(e).__name__}'; return res
    edfp = os.path.join(clip_dir, 'eeg.edf')
    if not os.path.exists(edfp):
        res['status'] = 'no_edf'; return res
    try:
        h = edflib.read_header(edfp)
    except Exception as e:
        res['status'] = f'bad_edf:{type(e).__name__}'; return res

    ch = biopotential_index(h)
    if ch is None:
        res['status'] = 'no_biopotential_channel'; return res
    res['channel_index_used'] = ch
    res['channel_label'] = h.labels[ch]
    res['source_start_from_edf_header'] = h.start.isoformat()
    res['clip_start_from_info'] = meta['clip_start'].isoformat() if meta['clip_start'] else None
    res['info_edf'] = meta['edf']
    fs = h.sfreqs[ch]
    res['clip_sfreq'] = fs

    files = csv_for(animal, day)
    if not files:
        res['status'] = 'no_export_csv'; return res

    off = (meta['clip_start'] - h.start).total_seconds()
    worst = {}
    matched_src, matched_sr, matched_unit, matched_patient = set(), set(), set(), set()
    ann_checks = []
    n = 0
    for cf in files:
        with open(cf, newline='') as f:
            for r in csv.DictReader(f):
                try:
                    s0, s1 = float(r['epoch_start_s']), float(r['epoch_end_s'])
                except (TypeError, ValueError):
                    continue
                if not (off <= s0 and s1 <= off + h.duration_s):
                    continue
                # the CSV's own wall clock must agree with the source EDF header
                nat = r['epoch_start_native']
                if nat:
                    want = h.start + dt.timedelta(seconds=s0)
                    got = dt.datetime.fromisoformat(nat)
                    res.setdefault('native_clock_delta_s', []).append((got - want).total_seconds())
                matched_src.add(r['source_file']); matched_sr.add(r['sample_rate_hz'])
                matched_unit.add(r['source_unit']); matched_patient.add(r['source_patient'])
                x, _ = edflib.read_signal(edfp, ch, s0 - off, s1 - s0, header=h)
                if x.size != int(round((s1 - s0) * fs)):
                    res['notes'].append(f'short read at epoch {r["epoch_index"]}'); continue
                got = features(x, fs)
                for k, gv in got.items():
                    try:
                        rv = float(r[k])
                    except (TypeError, ValueError, KeyError):
                        continue
                    scale = abs(got['std_uV']) if k in ('mean_uV',) else abs(rv)
                    denom = max(scale, 1e-12)
                    rel = abs(gv - rv) / denom
                    if k not in worst or rel > worst[k][0]:
                        worst[k] = (rel, rv, gv, r['epoch_index'])
                # annotation-overlap timing check against our own info.txt
                if meta['sz_start'] and meta['sz_end']:
                    a = max(h.start + dt.timedelta(seconds=s0), meta['sz_start'])
                    b = min(h.start + dt.timedelta(seconds=s1), meta['sz_end'])
                    exp = max(0.0, (b - a).total_seconds()) / (s1 - s0)
                    try:
                        ann_checks.append(abs(float(r['annotation_overlap_fraction']) - exp))
                    except (TypeError, ValueError):
                        pass
                n += 1
                if max_epochs and n >= max_epochs:
                    break
    res['n_epochs'] = n
    res['errors'] = {k: dict(rel=v[0], csv=v[1], recomputed=v[2], epoch=v[3]) for k, v in worst.items()}
    res['max_rel_exact'] = max([worst[k][0] for k in EXACT if k in worst], default=None)
    res['max_rel_spectral'] = max([worst[k][0] for k in SPECTRAL if k in worst], default=None)
    res['source_file'] = sorted(matched_src)
    res['sample_rate_hz'] = sorted(matched_sr)
    res['source_unit'] = sorted(matched_unit)
    res['source_patient'] = sorted(matched_patient)
    res['edf_patient'] = h.patient
    res['max_native_clock_delta_s'] = (max(abs(d) for d in res['native_clock_delta_s'])
                                       if res.get('native_clock_delta_s') else None)
    res.pop('native_clock_delta_s', None)
    res['max_annotation_overlap_err'] = max(ann_checks) if ann_checks else None
    if n == 0:
        res['status'] = 'no_contained_epochs'
    return res


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--animal', required=True)
    ap.add_argument('--day')
    ap.add_argument('--max-clips', type=int, default=3)
    ap.add_argument('--clips-per-day', type=int, default=None,
                    help='sample this many clips from EVERY day instead of --max-clips overall')
    ap.add_argument('--max-epochs', type=int, default=None)
    a = ap.parse_args()
    root = os.path.join(REPO, 'data_full', f'Data_{a.animal}_cropped')
    days = [a.day] if a.day else sorted(os.listdir(root))
    out = []
    for day in days:
        clips = sorted(glob.glob(os.path.join(root, day, 'seizure_*')))
        k = a.clips_per_day if a.clips_per_day else a.max_clips
        for c in clips[:k]:
            out.append(check_clip(c, a.animal, day, a.max_epochs))
    json.dump(out, sys.stdout, indent=1, default=str)


if __name__ == '__main__':
    main()
