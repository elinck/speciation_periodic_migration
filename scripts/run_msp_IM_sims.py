#msprime simulations for constant v cyclical gene flow
import os, sys, msprime as msp, numpy as np, re, allel
os.chdir("/home/cbattey2/speciation_cyclical_migration")
#os.chdir("/Users/cj/Dropbox/speciation_cyclical_migration")

#run with -> parallel -j 50 -N0 python scripts/run_msp_IM_sims.py ::: {1..100}

###############################################
############# periodic migration ##############
nsamples=40 #n haploid samples per population
Ne=1e5 #size of each subpopulation

#get years where CO2 passes 250ppm
vostok=np.loadtxt("Vostok.txt")
CO2cutoff=250
startrate=0;switchyears=[];rates=[]
for i in vostok:
    if i[3] < 250:
        newrate=1e-4
    else:
        newrate=0
    if not startrate==newrate:
        switchyears.append(i[2])
        rates.append(newrate)
    startrate=newrate

#two starting populations
pop_configs=[msp.PopulationConfiguration(sample_size=nsamples,
                                         initial_size=Ne),
             msp.PopulationConfiguration(sample_size=nsamples,
                                         initial_size=Ne)]

#switch migration to 10 migrants per generation when CO2 is below 250ppm
dem_events=[]
for i in range(len(switchyears)):
    dem_events.append(msp.MigrationRateChange(time=switchyears[i],
                                             rate=rates[i]))

dem_events.append(msp.MassMigration(time=6e5,source=1,destination=0,proportion=1)) #populations merge 600kybp

#simulate
print("\nsimulating periodic migration")
trees=msp.simulate(population_configurations=pop_configs,
                   demographic_events=dem_events,
                   migration_matrix=[[0,0],[0,0]],
                   mutation_rate=2.3e-9,
                   recombination_rate=1e-8,
                   length=1e7)
simID=int(np.random.choice(np.arange(1e7,9e7,1),1)[0])
#trees.dump("simulations/msprime/periodic/vostok_sim_"+str(simID)+".trees")
with open("simulations/msprime/periodic/vostok_sim_"+str(simID)+".vcf","w") as out:
    trees.write_vcf(out,ploidy=2)
#calculate mean fst to set migration for constant-migration simulation
haps=np.array(trees.genotype_matrix())
positions=np.array([s.position for s in trees.sites()])
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
fst=allel.stats.fst.average_weir_cockerham_fst(genotypes,[list(range(0,int(nsamples/2))),list(range(int(nsamples/2),nsamples))],100)[0]
print("\nfst="+str(fst))

###########################################
########### constant migration ############
m=((1/fst)-1)/(8*Ne) #corrected for 2 populations -- from Latter 1973 but see Wang 2004 Cons Biol for this form

pop_configs=[msp.PopulationConfiguration(sample_size=nsamples,
                                         initial_size=Ne),
             msp.PopulationConfiguration(sample_size=nsamples,
                                         initial_size=Ne)]

dem_events=[msp.MassMigration(time=6e5,source=1,destination=0,proportion=1)]

#simulate
print("\nsimulating constant migration")
trees=msp.simulate(population_configurations=pop_configs,
                   demographic_events=dem_events,
                   migration_matrix=[[0,m],
                                     [m,0]],
                   mutation_rate=2.3e-9,
                   recombination_rate=1e-8,
                   length=1e8)
#simID=int(np.random.choice(np.arange(1e7,9e7,1),1)[0])
#trees.dump("simulations/msprime/constant/constant_sim_"+str(simID)+".trees")
with open("simulations/msprime/constant/vostok_sim_"+str(simID)+".vcf","w") as out:
    trees.write_vcf(out,ploidy=2)
haps=np.array(trees.genotype_matrix())
positions=np.array([s.position for s in trees.sites()])
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
fst=allel.stats.fst.average_weir_cockerham_fst(genotypes,[list(range(0,int(nsamples/2))),list(range(int(nsamples/2),nsamples))],100)[0]
print("\nfst="+str(fst))
