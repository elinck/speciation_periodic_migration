#msprime simulations for constant v cyclical gene flow
import os, sys, msprime as msp, numpy as np, re
os.chdir("~/Dropbox/speciation_periodic_migration")

###############################################
############# periodic migration ##############
#get years where CO2 passes 250ppm
vostok=np.loadtxt("data/Vostok.txt")
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
dem_events=[]

#two starting populations
pop_configs=[msp.PopulationConfiguration(sample_size=20,
                                         initial_size=1e5),
             msp.PopulationConfiguration(sample_size=20,
                                         initial_size=1e5)]

#switch migration to 10 migrants per generation (m=1e-4) whenever CO2 is below 250ppm
for i in range(len(switchyears)):
    dem_events.append(msp.MigrationRateChange(time=switchyears[i],
                                             rate=rates[i]))
dem_events.append(msp.MassMigration(time=6e5,source=1,destination=0,proportion=1)) #populations merge 600k ybp

#simulate
trees=msp.simulate(population_configurations=pop_configs,
                   demographic_events=dem_events,
                   migration_matrix=[[0,0],[0,0]],
                   mutation_rate=0,
                   recombination_rate=1e-8,
                   length=5e7)
simID=int(np.random.choice(np.arange(1e7,9e7,1),1)[0])
trees.dump("simulations/msprime/periodic/vostok_sim_"+str(simID)+".trees")


###########################################
########### constant migration ############
#final Fst in period sims is around 0.06, so setting m to 3.574587e-05/2 to produce similar equilibrium Fst in constant-migration sims
pop_configs=[msp.PopulationConfiguration(sample_size=20,
                                         initial_size=1e5),
             msp.PopulationConfiguration(sample_size=20,
                                         initial_size=1e5)]

dem_events=[msp.MassMigration(time=6e5,source=1,destination=0,proportion=1)]

#simulate
trees=msp.simulate(population_configurations=pop_configs,
                   demographic_events=dem_events,
                   migration_matrix=[[0,3.574587e-05/2],[3.574587e-05/2,0]],
                   mutation_rate=0,
                   recombination_rate=1e-8,
                   length=5e7)
simID=int(np.random.choice(np.arange(1e7,9e7,1),1)[0])
trees.dump("simulations/msprime/constant/constant_sim_"+str(simID)+".trees")
