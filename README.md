# theWeights-examples
Example code implementing the parental weighting scheme outlined the manuscript "TheWeight: A simple and flexible algorithm for simulating non-ideal populations"

# the R directory 
This directory contains R code implementing three example weighting schemes: 
  - 1) The same weights for all parents, equivalent to Wright-Fisher reproduction.
  - 2) Sequential weights from 1, 2, ... N for N individuals
  - 3) Inverse weights from 1/1, 1/2, ... 1/N for N individuals
  
  To run this code use Rscript, Rstudio, or simply paste the contents of sim.R into an active R session. 
  
 # the SLiM directory
 This directory contains code implementing the sequential weights (1, 2, ... N for N individuals) in SLiM.
 
To run this code, first run some replicate simulations in SLiM: 
```
# (at a bash shell )
for REP in {1..100};
do 
 slim -d N=100 -d L=100 -d REP=$REP fromfile.test.slim
done
```
This will simulate 100 replicates of 100 loci (L=100) in 100 individuals (N=100) for 20 generations, and put the results in res.txt. 


Then use the supplied R code to summarize and plot the results:
```
# (at a bash shell )

Rscript summary_and_plot.R

```
