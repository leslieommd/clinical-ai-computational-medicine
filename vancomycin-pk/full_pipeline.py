#!/usr/bin/env python3
"""Educational two-compartment PK workflow for synthetic vancomycin courses.

Usage: python full_pipeline.py --input data --output results [--limit 5]
Dependencies: numpy, scipy. Fits four parameters and projects seven regimens.
"""
import argparse, csv, json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares


def read_csv(folder, name):
    with (folder/name).open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
def stamp(s): return datetime.fromisoformat(s.replace('Z','+00:00')).timestamp()/3600

def simulate(times, infusion, renal, pars, initial=None, auc_start=None):
    """Exact exponential transition on every discontinuity; optional AUC state."""
    cl,vc,vp,q = pars
    t0=min(times)
    ends=max(times)
    boundaries=sorted(set([t0,ends,*times]+[t for t,_ in infusion if t0<t<ends]+[t for t,_ in renal if t0<t<ends]))
    x=np.zeros(3) if initial is None else np.array([initial[0],initial[1],0.],float)
    answers={t0:x.copy()}
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        if b<=a:continue
        rate=sum(0 for _ in [])  # actual aggregate rate already supplied as piecewise changes
        for t,r in infusion:
            if t<=a+1e-10:rate=r
            else:break
        scr=renal[0][1]
        for t,c in renal:
            if t<=a+1e-10:scr=c
            else:break
        clearance=cl*(1/scr)**0.65
        # [central amount, peripheral amount, accumulated integral of central C]
        H=np.array([[-(clearance+q)/vc,q/vp,0], [q/vc,-q/vp,0], [1/vc if auc_start is None or a>=auc_start-1e-8 else 0,0,0.]])
        augmented=np.zeros((4,4));augmented[:3,:3]=H;augmented[0,3]=rate
        trans=expm(augmented*(b-a))
        x=(trans@np.r_[x,1])[:3]
        answers[b]=x.copy()
    return answers

def at_time(rows,t):
    return next(v for s,v in reversed(rows) if s<=t+1e-10)

def aggregate_events(events, anchor):
    discarded={e['revises_id'] for e in events if e['revises_id']}
    active=sorted((e for e in events if e['event_id'] not in discarded),key=lambda e:(e['effective_utc'],e['recorded_utc'],e['event_id']))
    inf_rates={}; total=0.; timeline=[(0.,0.)]
    for e in active:
        t=stamp(e['effective_utc'])-anchor
        k=e['infusion_id']; nr=float(e['rate_mg_h']) if e['action']=='SET' else 0.
        total+=nr-inf_rates.get(k,0.);inf_rates[k]=nr
        timeline.append((t,total))
    return sorted(timeline),len(discarded)

def prep_levels(records,anchor):
    versions={}
    for e in records:
        key=e['specimen_id']
        if key not in versions or int(e['revision'])>int(versions[key]['revision']):versions[key]=e
    valid=[]
    for e in versions.values():
        t=stamp(e['collection_utc'])-anchor
        if e['status']=='FINAL' and float(e['concentration_mg_l'])>0 and 0<t<=48:
            valid.append((t,float(e['concentration_mg_l'])))
    return sorted(valid),len(records)-len(valid)

def prior(weight): return np.array([3.8*(weight/75)**.75,.35*weight,.40*weight,4.2*(weight/75)**.75])
def bounds(weight): return np.array([1.,.20*weight,.20*weight,1.2]),np.array([11.,.65*weight,.90*weight,10.5])

def fit(levels,pump,renal,wt,starts=2):
    times=[0]+[t for t,_ in levels]
    goal=np.array([v for _,v in levels]);pr=prior(wt);lo,hi=bounds(wt)
    def residual(x):
        pars=np.exp(x)
        states=simulate(times,pump,renal,pars)
        pred=np.array([states[t][0]/pars[1] for t,_ in levels])
        return np.r_[(pred-goal)/.7,np.log(pars/pr)/[.6,.5,.35,.45]]
    starts_list=[np.clip(pr,lo,hi)]
    if starts>1:starts_list.append(np.sqrt(lo*hi))
    best=None
    for guess in starts_list:
        res=least_squares(residual,np.log(guess),bounds=(np.log(lo),np.log(hi)),max_nfev=160,ftol=1e-9,xtol=1e-9,gtol=1e-9)
        obj=float(np.dot(res.fun,res.fun))
        if best is None or obj<best[0]:best=(obj,np.exp(res.x))
    return best[1]

def regimen_schedule(option):
    # Project prescribed intervals as a piecewise-constant aggregate rate.
    changes=defaultdict(float)
    t=48.
    while t<72-1e-8:
        duration=float(option['infusion_h']);rate=float(option['dose_mg'])/duration
        changes[t]+=rate;changes[min(72.,t+duration)]-=rate
        t+=float(option['interval_h'])
    level=0.;out=[(48.,0.)]
    for t,delta in sorted(changes.items()):level+=delta;out.append((t,level))
    return out

def rank_options(scored,protocol):
    target=protocol['ranking'];lo,hi=target['target_auc_range'];mid=target['target_midpoint']
    def key(rec):
        auc=rec['auc'];o=rec['opt'];dose=0.;t=48.
        while t<72-1e-8:
            dose+=float(o['dose_mg'])*min(float(o['infusion_h']),72-t)/float(o['infusion_h'])
            t+=float(o['interval_h'])
        return (not(lo<=auc<=hi),abs(auc-mid),dose,o['option_id'])
    return sorted(scored,key=key)

def run(source,dest,limit=None):
    courses=read_csv(source,'courses.csv');pump=read_csv(source,'pump_events.csv');lev=read_csv(source,'drug_levels.csv');renal=read_csv(source,'renal_labs.csv')
    protocol=json.loads((source/'protocol.json').read_text())
    by_p=defaultdict(list);by_l=defaultdict(list);by_r=defaultdict(list)
    for e in pump:by_p[e['course_id']].append(e)
    for e in lev:by_l[e['course_id']].append(e)
    for e in renal:by_r[e['course_id']].append(e)
    exposures=[];regimens=[];replaced=excluded=0
    for c in courses[:limit]:
        course=c['course_id'];anchor=stamp(c['anchor_utc']);wt=float(c['weight_kg'])
        p,n=aggregate_events(by_p[course],anchor);replaced+=n
        levels,n=prep_levels(by_l[course],anchor);excluded+=n
        scr=sorted((stamp(r['collection_utc'])-anchor,float(r['creatinine_mg_dl'])) for r in by_r[course]);nlev=len(levels)
        if nlev<3:
            exposures.append([course,'','','',nlev,'insufficient_levels']);continue
        pars=fit(levels,p,scr,wt)
        # Evaluate pre-48 amount and AUC in a single exact state propagation.
        model=simulate([0.,24.,48.],p,scr,pars,auc_start=24.)
        obs=float(model[48.][2]);state=model[48.][:2]
        exposures.append([course,f'{pars[0]:.7f}',f'{pars[1]:.7f}',f'{obs:.7f}',nlev,'estimated'])
        estimates=[]
        for o in protocol['regimen_options']:
            scheduled=regimen_schedule(o)
            projected=simulate([48.,72.],scheduled,scr,pars,initial=state,auc_start=48.)
            estimates.append({'opt':o,'auc':float(projected[72.][2])})
        for index,row in enumerate(rank_options(estimates,protocol),1):
            regimens.append([course,row['opt']['option_id'],f'{row["auc"]:.7f}',index])
    dest.mkdir(parents=True,exist_ok=True)
    def write(name,header,data):
        with (dest/name).open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f,lineterminator='\n');w.writerow(header);w.writerows(data)
    write('exposure.csv',['course_id','cl_l_h','v_l','auc24_mg_h_l','n_levels','status'],exposures)
    write('regimen_scores.csv',['course_id','option_id','predicted_auc24_mg_h_l','rank'],regimens)
    (dest/'qc_report.json').write_text(json.dumps({'courses_processed':len(exposures),'records_replaced':replaced,'levels_excluded':excluded},indent=2)+'\n')
    print(f'Processed {len(exposures)} courses; {len(regimens)} regimen projections')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',default='data');parser.add_argument('--output',default='results');parser.add_argument('--limit',type=int)
    a=parser.parse_args();run(Path(a.input),Path(a.output),a.limit)
