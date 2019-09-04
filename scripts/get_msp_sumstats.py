import msprime as msp,pyslim,os,numpy as np,allel,itertools
from scipy import stats
from tqdm import tqdm
import argparse

#command line args
parser = argparse.ArgumentParser('get summary stats for msPrime grid simulations')
parser.add_argument('--basedir')
parser.add_argument('--n',type=int)
parser.add_argument('--outfile')
args=parser.parse_args()

## debug
# args = argparse.Namespace(basedir="/Users/cj/Dropbox/speciation_cyclical_migration/simulations/msprime/periodic/",
#                           n=0,
#                           outfile="/Users/cj/Desktop/test.txt")

#pairwise IBS tract distribution function for summary stats
def getPairwiseIbsTractLengths(x,y,positions,maxlen,min_len_to_keep=0):
    '''
    input:
    x: haplotype as 1D array
    y: haplotype as 1D array
    positions: a 1D array of SNP positions
    maxlen: length of chromosome/contig

    Returns:
    1d array listing distances between adjacent SNPs in a pair of sequences.
    Assumes all sites are accessible.
    '''
    snps=~np.equal(x,y)
    snp_positions=positions[snps]
    l=len(snp_positions)
    ibs_tracts=[]
    if(l==0):
        ibs_tracts=[maxlen]
    else:
        if(l>1):
            ibs_tracts=snp_positions[np.arange(1,l-1,1)]-snp_positions[np.arange(0,l-2,1)] #middle blocks
        np.append(ibs_tracts,snp_positions[0]+1)        #first block
        np.append(ibs_tracts,maxlen-snp_positions[l-1]) #last block
        con=[x>=min_len_to_keep for x in ibs_tracts] #drop blocks < min_len_to_keep
        ibs_tracts=np.extract(con,ibs_tracts)
    return ibs_tracts

#directory to summarize
os.chdir(args.basedir)

#get the target tree
trees=[x for x in os.listdir() if x.endswith(".trees") and not x.startswith(".")]
t=trees[args.n]
sumstats=[]

#add mutations
ts=msp.load(t)
ts1=msp.mutate(ts,2.3e-9)

#dump a vcf
with open(t+".vcf", "w") as vcf_file:
    ts1.write_vcf(vcf_file, 2)

###############summary stats##############
##data structures for stats (note assumes 20 samples per population)
pops={'0':[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],
      '1':[20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]}
haps=np.array(ts1.genotype_matrix())
positions=np.array([s.position for s in ts1.sites()])
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
allele_counts=genotypes.count_alleles()
subpop_allele_counts=genotypes.count_alleles_subpops(subpops=pops)
genotype_allele_counts=genotypes.to_allele_counts()

##SNP stats
segsites=np.shape(genotypes)[0]
pi=allel.sequence_diversity(positions,allele_counts,start=1,stop=1e7)
tajD=allel.tajima_d(ac=allele_counts,start=1,stop=1e7)
thetaW=allel.watterson_theta(pos=positions,ac=allele_counts,start=1,stop=1e7)
het_o=np.mean(allel.heterozygosity_observed(genotypes))
fst=allel.stats.fst.average_weir_cockerham_fst(genotypes,
     [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],
      [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]],
     100)[0]
dxy=allel.stats.diversity.sequence_divergence(positions,subpop_allele_counts['0'],subpop_allele_counts['1'])

##Identical-by-state haplotype block length stats (i.e. LD-type stats)
#between population comparisons
bw_pairs=list(itertools.product(pops['0'],pops['1']))
ibs_bw=[]
for j in bw_pairs:
    ibspair=getPairwiseIbsTractLengths(haps[:,j[0]],
                                       haps[:,j[1]],
                                       positions,
                                       1e7,
                                       2)
    ibs_bw.append(ibspair)
ibs_bw=[x for y in ibs_bw for x in y]
ibs_bw_mean=np.mean(ibs_bw)
ibs_bw_skew=stats.skew(ibs_bw)
ibs_bw_blocks_over_1e4=len([x for x in ibs_bw if x > 1e4])

#within-population comparisons (must be a better way to get within-pop indicex pairs but this works)
wi_0_pairs=list(itertools.product(pops['0'],pops['0']))
wi_1_pairs=list(itertools.product(pops['1'],pops['1']))
wi_pairs=[wi_0_pairs,wi_1_pairs]
wi_pairs=[a for b in wi_pairs for a in b]
ibs_wi=[]
for j in wi_pairs:
    ibspair=getPairwiseIbsTractLengths(haps[:,j[0]],
                                       haps[:,j[1]],
                                       positions,
                                       1e7,
                                       2)
    ibs_wi.append(ibspair)
ibs_wi=[x for y in ibs_wi for x in y]
ibs_wi_mean=np.mean(ibs_wi)
ibs_wi_skew=stats.skew(ibs_wi)
ibs_wi_blocks_over_1e4=len([x for x in ibs_wi if x > 1e4])

#output summary stats
header="segsites pi thetaW tajD het_o fst dxy ibs_bw_mean ibs_bw_skew ibs_bw_blocks_over_1e4 ibs_wi_mean ibs_wi_skew ibs_wi_blocks_over_1e4\n"
if not os.path.exists(args.outfile):
    out=open(args.outfile,'a')
    out.write(header)
    out.write(str(segsites)+ " "+str(pi)+" "+str(thetaW)+" "+str(tajD)+" "+
              str(het_o)+" "+str(fst)+" "+str(dxy)+" "+str(ibs_bw_mean)+" "+
              str(ibs_bw_skew)+" "+str(ibs_bw_blocks_over_1e4)+" "+
              str(ibs_wi_mean)+" "+str(ibs_wi_skew)+" "+
              str(ibs_wi_blocks_over_1e4)+"\n")
    out.close()
else:
    out=open(args.outfile,'a')
    out.write(str(segsites)+ " "+str(pi)+" "+str(thetaW)+" "+str(tajD)+" "+
              str(het_o)+" "+str(fst)+" "+str(dxy)+" "+str(ibs_bw_mean)+" "+
              str(ibs_bw_skew)+" "+str(ibs_bw_blocks_over_1e4)+" "+
              str(ibs_wi_mean)+" "+str(ibs_wi_skew)+" "+
              str(ibs_wi_blocks_over_1e4)+"\n")
    out.close()
