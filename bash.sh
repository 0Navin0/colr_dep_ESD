#!/bin/bash
#config1=(~/configs_n_signals2/fullpofz/run0_random_rbin5_fullpofz/signal_dr72safe/signal_dr72safe*)
#config2=(~/configs_n_signals3/fullpofz/run1_noRandom_rbin5_fullpofz/signal_dr72safe/signal_dr72safe*)
#for ((i=0;i<8;i++)) ; do diff -s "${config1[i]}" "${config2[i]}"; done
config1="find /home/navin/git/colr_dep_ESD/Obesrved_signal/*/*/*/*png"
config2="find /home/navin/git//Obesrved_signal/*/*/*/*png"
for ((i=0;i<287;i++)) ; do diff -s "${config1[i]}" "${config2[i]}"; done



##some bash commands
#awk 'NR!=1{print $6,$8,$12}'
