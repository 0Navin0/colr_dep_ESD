#!/bin/bash
#config1=(~/configs_n_signals2/fullpofz/run0_random_rbin5_fullpofz/signal_dr72safe/signal_dr72safe*)
#config2=(~/configs_n_signals3/fullpofz/run1_noRandom_rbin5_fullpofz/signal_dr72safe/signal_dr72safe*)
#for ((i=0;i<8;i++)) ; do diff -s "${config1[i]}" "${config2[i]}"; done
config1="find /home/navin/git/colr_dep_ESD/Obesrved_signal/*/*/*/*png"
config2="find /home/navin/git//Obesrved_signal/*/*/*/*png"
for ((i=0;i<287;i++)) ; do diff -s "${config1[i]}" "${config2[i]}"; done



##some bash commands
#awk:
#awk 'NR!=1{print $6,$8,$12}'

##brace expansion
#touch z{1..10}.c
#touch z{a..z}.c
#echo *{1..10}*
#echo {001..009}
#echo {A..Z}{0..9}
#echo {{A..Z},{a..z}}
#echo -e \\n{{A..Z},{a..z}}
## ues cases:
#wget http://docs.example.com/documentation/slides_part{1..6}.html
#mkdir /home/bash/test/{foo,bar,baz,cat,dog}   ## command to remove only the directories: rm -r */ OR rm -rf `ls -d */`
#echo img{00{1..9},0{10..99},{100..999}}.png
#printf "%s\n" img{00{1..9},0{10..99},{100..999}}.png
