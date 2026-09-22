import os, subprocess, sys
import os
_DEFAULT_REPO = os.environ.get('EEG_ROOT', '/work/mech-ai-scratch/alloy/EEG')
REPO = os.environ.get('EEG_ROOT', _DEFAULT_ROOT)  # override for another machine
report=[]
for z in ('RoomC.zip','RoomD.zip'):
    out = subprocess.run(['unzip','-l',os.path.join(REPO,z)],capture_output=True,text=True).stdout
    entries={}
    for line in out.splitlines():
        p=line.split(None,3)
        if len(p)==4 and p[0].isdigit() and '-' in p[1]:
            entries[p[3]]=int(p[0])
    missing=[];mismatch=[];extra=[]
    for name,size in entries.items():
        if name.endswith('/'): continue
        disk=os.path.join(REPO,'raw_epoch_features',name)
        if not os.path.exists(disk):
            # tolerate the DONE_ prefix the user added to finished animals
            parts=name.split('/')
            alt=None
            if len(parts)>=3:
                alt=os.path.join(REPO,'raw_epoch_features',parts[0],'DONE_'+parts[1],*parts[2:])
            if alt and os.path.exists(alt): disk=alt
            else: missing.append(name); continue
        ds=os.path.getsize(disk)
        if ds!=size: mismatch.append((name,size,ds))
    report.append((z,len(entries),len(missing),len(mismatch),missing[:8],mismatch[:8]))
for z,n,nm,nx,miss,mis in report:
    print(f'{z}: {n} entries | missing on disk: {nm} | size mismatch: {nx}')
    for m in miss: print('   MISSING ',m)
    for m in mis: print('   SIZE    ',m)
