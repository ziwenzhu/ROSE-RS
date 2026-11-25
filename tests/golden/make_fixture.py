"""Regenerates the synthetic golden fixture: python3 make_fixture.py <outdir> (needs pysam)."""
import random, pysam, os
random.seed(7)
import sys
OUT=sys.argv[1]
chroms={'chr1':400_000,'chr2':250_000,'chrY':40_000}
# genes
genes=[]; gi=0
for c,L in chroms.items():
    pos=20000
    while pos < L-60000:
        gi+=1
        strand=random.choice('+-'); ln=random.randint(5000,40000)
        genes.append((f"NM_{gi:05d}",c,strand,pos,pos+ln,f"GENE{gi}"))
        if random.random()<0.3:  # alt transcript same gene different TSS
            gi+=1; genes.append((f"NM_{gi:05d}",c,strand,pos+3000,pos+ln,f"GENE{gi-1}"))
        pos+=ln+random.randint(5000,30000)
# a dense 4-gene cluster on chr1 near 1.5Mb
for k in range(4):
    gi+=1; genes.append((f"NM_{gi:05d}",'chr1','+',200_000+k*4000,200_000+k*4000+2000,f"DENSE{k}"))
with open(OUT+'/hg19_refseq.ucsc','w') as f:
    f.write("#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames\n")
    for (n,c,s,a,b,g) in genes:
        f.write(f"0\t{n}\t{c}\t{s}\t{a}\t{b}\t{a}\t{b}\t1\t{a},\t{b},\t0\t{g}\tcmpl\tcmpl\t0,\n")
# peaks
peaks=[]; pid=0
def addpeak(c,s,w):
    global pid; pid+=1; peaks.append((c,s,s+w,f"peak_{pid}"))
for c,L in chroms.items():
    for _ in range(L//25000):
        s=random.randint(1000,L-5000); addpeak(c,s,random.randint(200,2500))
    for _ in range(max(1,L//100000)):  # clusters
        base=random.randint(10000,max(10001,L-80000))
        for k in range(random.randint(3,9)):
            addpeak(c,base,random.randint(500,3000)); base+=random.randint(1000,14000)
# promoter peaks
for (n,c,s,a,b,g) in random.sample(genes,6):
    t=a if s=='+' else b
    addpeak(c,t-300,600)
# cluster across DENSE genes
for k in range(5): addpeak('chr1',199_000+k*3500,1500)
# exact-gap edge cases for stitch boundary: gaps of 12499,12500,12501 (closed-coord) 
for g in (12499,12500,12501):
    base=random.randint(20000,200000)
    addpeak('chr2',base,1000); addpeak('chr2',base+1000+g,1000)
peaks.sort()
with open(OUT+'/peaks.gff','w') as f:
    for c,s,e,i in peaks:
        f.write(f"{c}\t{i}\t\t{s}\t{e}\t\t.\t\t{i}\n")
with open(OUT+'/peaks.bed','w') as f:
    for c,s,e,i in peaks: f.write(f"{c}\t{s-1}\t{e}\t{i}\n")
hdr={'HD':{'VN':'1.6','SO':'unsorted'},'SQ':[{'SN':c,'LN':L} for c,L in chroms.items()]}
def write_bam(fn,n_bg,enrich):
    tmp=fn+'.tmp.bam'
    with pysam.AlignmentFile(tmp,'wb',header=hdr) as out:
        k=0
        def emit(c,p,rev,dup=False):
            nonlocal k
            k+=1; a=pysam.AlignedSegment(); a.query_name=f"r{k}"; a.reference_id=list(chroms).index(c)
            a.reference_start=max(0,p); a.mapping_quality=random.choice([0,30,60,60,60]); a.cigarstring='50M'
            a.query_sequence='A'*50; a.query_qualities=pysam.qualitystring_to_array('I'*50)
            a.flag=(16 if rev else 0)|(1024 if dup else 0); out.write(a)
        for c,L in chroms.items():
            for _ in range(int(n_bg*L/5.5e6)):
                emit(c,random.randint(0,L-60),random.random()<.5)
        if enrich:
            for c,s,e,i in peaks:
                strength=random.lognormvariate(2,1.2)
                for _ in range(int(strength*(e-s)/200)):
                    emit(c,random.randint(s-150,e),random.random()<.5,dup=random.random()<0.1)
    pysam.sort('-o',fn,tmp); os.remove(tmp); pysam.index(fn)
write_bam(OUT+'/chip.bam',60000,True)
write_bam(OUT+'/input.bam',60000,False)
with open(OUT+'/ref.fa','w') as f:
    for c,L in chroms.items():
        f.write('>'+c+'\n')
        for i in range(0,L,60): f.write('ACGT'[(i//60)%4]*min(60,L-i)+'\n')
pysam.faidx(OUT+'/ref.fa')
pysam.view('-C','-T',OUT+'/ref.fa','-o',OUT+'/chip.cram',OUT+'/chip.bam',catch_stdout=False)
pysam.index(OUT+'/chip.cram')
print(len(genes),'genes',len(peaks),'peaks')
