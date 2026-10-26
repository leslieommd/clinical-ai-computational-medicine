#!/usr/bin/env python3
"""Educational implementation of autosomal-recessive haplotype reconstruction.
Usage: python solution.py --input /workspace/input --output /workspace/output
"""
import argparse
import collections
import csv
import itertools
import json
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument('--input', default='/workspace/input')
parser.add_argument('--output', default='/workspace/output')
args=parser.parse_args()
input_dir=Path(args.input)
output_dir=Path(args.output)
output_dir.mkdir(parents=True, exist_ok=True)

rows={}
for name in ('pedigree.csv','variants.csv','cnvs.csv','phase.csv','genes.csv'):
    with (input_dir/name).open(newline='',encoding='utf-8') as f:
        rows[name]=list(csv.DictReader(f))
for q in rows['cnvs.csv']:
    q['copy_number']=int(q['copy_number'])
if any(q['inheritance'] != 'AR' for q in rows['genes.csv']):
    raise ValueError('All supplied genes must be autosomal recessive')

by_sample=collections.defaultdict(dict)
for r in rows['pedigree.csv']: by_sample[r['sample_id']]=r
var=collections.defaultdict(list); cnv=collections.defaultdict(list); ph=collections.defaultdict(list)
for r in rows['variants.csv']: var[(r['sample_id'],r['gene'])].append(r)
for r in rows['cnvs.csv']: cnv[(r['sample_id'],r['gene'])].append(r)
for r in rows['phase.csv']: ph[r['sample_id']].append(r)

def overlap(a,b):
 st=max(int(a['start']),int(b['start']));en=min(int(a['end']),int(b['end']))
 inter=max(0,en-st+1)
 return inter/(int(a['end'])-int(a['start'])+1)>=.8 and inter/(int(b['end'])-int(b['start'])+1)>=.8

def merged_cnv_for(sample,gene):
 raw=cnv[(sample,gene)]
 parents=list(range(len(raw)))
 def rt(x):
  while parents[x]!=x: x=parents[x]
  return x
 for i,j in itertools.combinations(range(len(raw)),2):
  if raw[i]['copy_number']==raw[j]['copy_number'] and raw[i]['classification']==raw[j]['classification'] and overlap(raw[i],raw[j]):
   parents[rt(j)]=rt(i)
 groups=collections.defaultdict(list)
 for i in range(len(raw)): groups[rt(i)].append(raw[i])
 # Keep a canonical row (without '-DUP' if possible)
 ans={}
 for rr in groups.values():
  orig=sorted(rr,key=lambda r: (r['event_id'].endswith('-DUP'), r['event_id']))[0]
  ans[orig['event_id']] = (orig,rr)
 return ans

out_assign=[];out_diagnosis=[]; conflicts=0; cnv_merges=0; total_ambig=0; checked=0
for r in rows['pedigree.csv']:
 C=r['sample_id']
 if not (r['mother_id'] and r['father_id']):continue
 M,F=r['mother_id'],r['father_id']
 has_genes=sorted(set([g for s,g in var if s==C] +[g for s,g in cnv if s==C]))
 for gene in has_genes:
  cv=merged_cnv_for(C,gene); mcv=merged_cnv_for(M,gene);fcv=merged_cnv_for(F,gene)
  cnv_merges+=sum(len(rr)-1 for _,rr in cv.values())
  child_v=var[(C,gene)]
  mom_v={q['event_id']:q for q in var[(M,gene)]}
  dad_v={q['event_id']:q for q in var[(F,gene)]}
  candidates={};both=[];klass={}; evid={};overlaps=[]
  for q in child_v:
   eid=q['event_id'];klass[eid]=q['classification'];evid[eid]=[q]
   if q['genotype']=='1/1': both.append(eid);continue
   mo=eid in mom_v; fa=eid in dad_v
   candidates[eid]=set((['maternal'] if mo else [])+(['paternal'] if fa else []))
  for eid,(q,rr) in cv.items():
   evid[eid]=rr;klass[eid]=q['classification']
   if q['copy_number']==0: both.append(eid);continue
   mo=any((x['copy_number']==1 or x['copy_number']==0) and overlap(q,x) for x,grp in mcv.values())
   fa=any((x['copy_number']==1 or x['copy_number']==0) and overlap(q,x) for x,grp in fcv.values())
   candidates[eid]=set((['maternal'] if mo else [])+(['paternal'] if fa else []))
  # CNV and SNP overlapping constraints imply opposite inherited alleles.
  for eid,(q,rr) in cv.items():
   if eid not in candidates:continue
   for v in child_v:
    if v['event_id'] in candidates and int(q['start'])<=int(v['position'])<=int(q['end']): overlaps.append((eid,v['event_id'],'trans'))
  edges=[(p['event_a'],p['event_b'],p['relationship']) for p in ph[C] if p['event_a'] in candidates and p['event_b'] in candidates]
  edges+=overlaps
  # Enumerate connected phase components with parity constraints, honoring parental information.
  adj=collections.defaultdict(list)
  for a,b,rel in edges:
   parity=0 if rel=='cis' else 1
   adj[a].append((b,parity));adj[b].append((a,parity))
  color={}; components=[];bad=False
  for node in candidates:
   if node in color: continue
   color[node]=0; stack=[node]; component=[]
   while stack:
    at=stack.pop();component.append(at)
    for nxt,p in adj[at]:
     x=color[at]^p
     if nxt not in color:color[nxt]=x;stack.append(nxt)
     elif color[nxt]!=x:bad=True
   components.append(component)
  choices=[]
  for comp in components:
   flips=[]
   for flip in [0,1]:
    if all(('maternal' if color[x]^flip==0 else 'paternal') in candidates[x] for x in comp): flips.append(flip)
   if not flips:bad=True
   choices.append(flips)
  if bad:
   conflicts+=1
   assignments={k:'unresolved' for k in candidates}
   status='unresolved'
  else:
   assignments={}
   for comp,flips in zip(components,choices):
    for x in comp:
     origins={('maternal' if color[x]^flip==0 else 'paternal') for flip in flips}
     assignments[x]=next(iter(origins)) if len(origins)==1 else 'unresolved'
   # Diagnosis invariant over all logically consistent phase assignments.
   outcome=set()
   for flip_tuple in itertools.product(*choices):
    disrupted={'maternal':False,'paternal':False}
    if any(klass[eid] in ('P','LP') for eid in both):disrupted={'maternal':True,'paternal':True}
    for comp,flip in zip(components,flip_tuple):
     for eid in comp:
      if klass[eid] in ('P','LP'):
       disrupted['maternal' if color[eid]^flip==0 else 'paternal']=True
    outcome.add(disrupted['maternal'] and disrupted['paternal'])
   status='unresolved' if len(outcome)>1 else ('diagnosed' if True in outcome else 'negative')
  for q in child_v:
   eid=q['event_id']
   out_assign.append(dict(sample_id=C,gene=gene,event_id=eid,assignment='both' if eid in both else assignments[eid]))
  for eid,(q,rr) in cv.items():
   assignment='both' if eid in both else assignments[eid]
   for raw in rr:out_assign.append(dict(sample_id=C,gene=gene,event_id=raw['event_id'],assignment=assignment))
  out_diagnosis.append(dict(sample_id=C,gene=gene,status=status))
  total_ambig+=(status=='unresolved');checked+=1

out_assign.sort(key=lambda x:(x['sample_id'],x['gene'],x['event_id']))
out_diagnosis.sort(key=lambda x:(x['sample_id'],x['gene']))

for name,data,headers in (
    ('resolved_alleles.csv',out_assign,['sample_id','gene','event_id','assignment']),
    ('family_diagnoses.csv',out_diagnosis,['sample_id','gene','status']),
):
    with (output_dir/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=headers,lineterminator='\n')
        w.writeheader()
        w.writerows(data)
with (output_dir/'qc_report.json').open('w',encoding='utf-8') as f:
    json.dump({'phase_conflicts':conflicts,'cnv_merges':cnv_merges,'unresolved_sample_genes':total_ambig},f,indent=2)
print('allele records:',len(out_assign),'patient-gene records:',len(out_diagnosis))
print('phase conflicts:',conflicts,'cnv merges:',cnv_merges,'unresolved patient-genes:',total_ambig)
