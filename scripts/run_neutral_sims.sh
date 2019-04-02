#!/bin/bash

cd ~/speciation_cyclical_migration

for i in {1..50}
	do 
	#~/SLiM_build/slim slim_recipes/IM_flux_vostok.slim &
	~/SLiM_build/slim slim_recipes/IM_constant.slim &
done
