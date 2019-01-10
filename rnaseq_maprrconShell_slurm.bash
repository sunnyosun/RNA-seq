#!/bin/bash                                                                                                                
#$ -S /bin/bash                                                                                                                        
#$ -cwd

for i in $(cat ../name.txt); do qsub rnaseq_maprrcon_star.bash $i 4; done