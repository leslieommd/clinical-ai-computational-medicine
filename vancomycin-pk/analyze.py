"""Audit synthetic therapeutic drug-monitoring event records.
Run: python analyze.py
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv, json
DATA = Path(__file__).parent / 'data'
def read(name):
    with (DATA/name).open(encoding='utf-8', newline='') as f: return list(csv.DictReader(f))
C, P, L, R = [read(n) for n in ('courses.csv','pump_events.csv','drug_levels.csv','renal_labs.csv')]
anchors={r['course_id']:datetime.fromisoformat(r['anchor_utc'].replace('Z','+00:00')) for r in C}
ids={r['event_id'] for r in P}
replaced={r['revises_id'] for r in P if r['revises_id']}
assert replaced <= ids, 'Unresolved pump-event revision reference'
# Apply latest revision per specimen, as determined by integer revision numbers.
specimens={}
for r in L:
    key=(r['course_id'],r['specimen_id'])
    if key not in specimens or int(r['revision'])>int(specimens[key]['revision']): specimens[key]=r
accepted=[]
for r in specimens.values():
    collection=datetime.fromisoformat(r['collection_utc'].replace('Z','+00:00'))
    hours=(collection-anchors[r['course_id']]).total_seconds()/3600
    if r['status']=='FINAL' and float(r['concentration_mg_l'])>0 and 0<hours<=48:
        accepted.append(r)
counts=Counter(r['course_id'] for r in accepted)
summary={
 'courses':len(C),
 'pump_events':len(P),
 'superseded_pump_events':len(replaced),
 'drug_level_records':len(L),
 'accepted_specimens':len(accepted),
 'excluded_level_records':len(L)-len(accepted),
 'renal_measurements':len(R),
 'insufficient_level_courses':sum(counts[r['course_id']]<3 for r in C),
}
print(json.dumps(summary,indent=2))
