# Speciation with periodic gene flow  
  
Scripts and processed data for Linck and Battey 2019.
  
### Scripts:  
  
(in `/scripts` subdirectory)  
  
`IM_constant.slim:` Isolation-with-migration SLiM simulation; not used in ms v1.  
    
`IM_flux.slim:` Periodic migration SLiM simulation; not used in ms v1.  
  
`IM_moments.py:` `moments` isolation-with-migration demographic model. 
  
`PM_moments.py:` `moments` periodic gene flow demographic model.  
  
`allopatry_BDMI_model.slim:` Simple BDMI model validation, allopatric speciation.  
  
`allopatry_model_loop.sh:` Loop to run simple BDMI model validation, allopatric speciation.  
  
`flux01_BDMI_model.slim:` Simple BDMI model validation, speciation with periodic gene flow, 1% divergence time in allopatry.  
   
`flux01_model_loop.sh:` Loop to run simple BDMI model validation, speciation with periodic gene flow, 1% divergence time in allopatry.  
  
`flux10_BDMI_model.slim:` Simple BDMI model validation, speciation with periodic gene flow, 10% divergence time in allopatry.
   
`flux10_model_loop.sh:` Loop to run simple BDMI model validation, speciation with periodic gene flow, 10% divergence time in allopatry. 

`flux50_BDMI_model.slim:` Simple BDMI model validation, speciation with periodic gene flow, 50% divergence time in allopatry.  
   
`flux50_model_loop.sh:` Loop to run simple BDMI model validation, speciation with periodic gene flow, 50% divergence time in allopatry.  

`get_msp_sumstats.py:` Get summary statistics from `msprime` simulations.   
  
`parapatry_BDMI_model.slim:` Simple BDMI model validation, parapatric speciation.   
  
`parapatry_model_loop.sh:`  Loop to run simple BDMI model validation, parapatric speciation.
  
`run_msp_IM_sims.py:` Run `msprime` isolation-with-migration simulations. 
  
`run_neutral_sims.sh:`  Run neutral SliM simulations; not used in ms v1. 
  
`sim_figs.R:`  Generate figures from simulations.
  
`simulation_waiting_time_speciation.R:` Process data and generate figures from simple model validation.   
  
`waiting_time_speciation.R:`  Plot simple model validation.   

### Data:  
    
(in `/data/simple_model_output/` subdirectory)   

`allopatry.txt:`  SLiM script output for simple BDMI model validation, allopatric speciation.  

`flux01.txt:`  SLiM script output for simple BDMI model validation, speciation with periodic gene flow, 1% divergence time in allopatry.  

`flux10.txt:`  SLiM script output for simple BDMI model validation, speciation with periodic gene flow, 10% divergence time in allopatry.  

`flux50.txt:`  SLiM script output for simple BDMI model validation, speciation with periodic gene flow, 50% divergence time in allopatry.  

`parapatry.txt:`  SLiM script output for simple BDMI model validation, parapatric speciation.  
  
(in `/data/sumstats/` subdirectory)   

`constant.txt:`  Summary statistics from isolation-with-migration `msprime` simulations.  
  
`periodic.txt:`  Summary statistics from periodic gene flow `msprime` simulations.  
  
(in `/data/vcf/` subdirectory)   

`constant_sim_10005626.trees.vcf:`  Randomly selected `.vcf` file from `msprime` isolation-with-migration simulations used for `moments` demographic inference.  

`vostok_sim_10388971.trees.vcf:`  Randomly selected `.vcf` file from `msprime` periodic gene flow simulations used for `moments` demographic inference.   
  
(in `/data/` subdirectory)   
    
`IM_params.txt:`  Parameters from isolation-with-migration `moments` model.    
  
`PM_params.txt:` Parameters from periodic gene flow `moments` model.  
  
`Vostok.txt:`  [Vostok Ice Core data](https://cdiac.ess-dive.lbl.gov/ftp/trends/co2/vostok.icecore.co2).  
  
`ID_pop_moments.txt:`  Population IDs for `moments.`
    
`samples:`  Plain text file of sample IDs.   
  
`sim_cycles_01.txt:`  Data for simulations w/ 1% total time in allopatry.   
   
`sim_cycles_10.txt:` Data for simulations w/ 10% total time in allopatry.    
    
`sim_cycles_50.txt:` Data for simulations w/ 50% total time in allopatry.    

  