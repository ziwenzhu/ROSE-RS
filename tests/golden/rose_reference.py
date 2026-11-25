"""Independent re-implementation of ROSE 0.1 stitching and per-region signal rules
(ROSE_main.py, ROSE_utils.py, ROSE_bamToGFF.py; floor 1, extension 200, rpm), used to make
expected_rose_region_map.tsv. Usage: python3 rose_reference.py peaks.gff chip.bam input.bam annotation tss stitch"""
import sys, pysam, re
from collections import defaultdict
gff, bam_fn, ctrl_fn, annot, tssw, stitch = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
peaks=[]
for l in open(gff):
    c=l.rstrip('\n').split('\t'); peaks.append((c[0],int(c[3]),int(c[4]),c[1]))
# refseq: first row per ID (importRefseq/makeStartDict)
start={}
rows=[l.rstrip('\n').split('\t') for l in open(annot)][1:]
for r in rows:
    if r[1] in start: continue
    tss=int(r[4]) if r[3]=='+' else int(r[5])
    start[r[1]]=(r[2],tss,r[12])
def tss_in(c,lo,hi,w):  # TSS loci [t-w,t+w] overlapping closed [lo,hi]
    return [v for v in start.values() if v[0]==c and v[1]-w<=hi and v[1]+w>=lo]
bound=[]
for p in peaks:
    c,lo,hi,i=p
    if tssw and any(t[1]-tssw<=lo and hi<=t[1]+tssw for t in start.values() if t[0]==c): continue
    bound.append(p)
# stitching: closed-interval overlap of [lo-w,hi+w] -> connected components
bychr=defaultdict(list)
for p in bound: bychr[p[0]].append(p)
regions=[]
for c,ps in bychr.items():
    ps.sort(key=lambda x:(x[1],x[2]))
    cl=[]; cur_hi=None
    for p in ps:
        if cl and p[1] <= cur_hi + stitch: cl.append(p); cur_hi=max(cur_hi,p[2])
        else:
            if cl: regions.append((c,cl))
            cl=[p]; cur_hi=p[2]
    if cl: regions.append((c,cl))
final=[]
for c,cl in regions:
    lo=min(p[1] for p in cl); hi=max(p[2] for p in cl)
    names={t[2] for t in tss_in(c,lo,hi,50)} if tssw else set()
    if len(names)>2:
        for p in cl: final.append((c,p[1],p[2],p[3]))
    else:
        final.append((c,lo,hi,f"{len(cl)}_{cl[0][3]}_lociStitched"))
final=[r for r in final if r[0]!='chrY']
def mmr(fn):
    s=pysam.samtools.flagstat(fn)
    for l in s.splitlines():
        if 'mapped (' in l: return round(int(l.split(' ')[0])/1e6,4)
def density(fn):
    m=mmr(fn); b=pysam.AlignmentFile(fn); out={}
    for (c,lo,hi,i) in final:
        h=defaultdict(int)
        for r in b.fetch(c, lo-200-1, hi+200):  # samtools view c:lo-200-hi+200 (1-based closed)
            if 'N' in r.cigarstring: continue
            P=r.reference_start+1; L=len(r.query_sequence)
            s,e=(P-200,P+L) if r.is_reverse else (P,P+L+200)
            for x in range(s,e+1): h[x]+=1
        keys=[x for x in h if h[x]>1 and lo<x<hi]
        bs=hi-lo
        out[i]= 0.0 if bs==0 else round(sum(h[x] for x in keys)/bs/m,4)*bs
    return out
t=density(bam_fn); k=density(ctrl_fn)
print("REGION_ID\tCHROM\tSTART\tSTOP\tchip\tinput")
for (c,lo,hi,i) in final: print(f"{i}\t{c}\t{lo}\t{hi}\t{t[i]}\t{k[i]}")
