# deepDILS
```
nIndividuals=40
nReplicates=100
nSNPs=2000
rho=300
length=100000

window_width=0.01
window_step=0.005

# simulated 100 replicates of a sweep
msms $nIndividuals $nReplicates -s $nSNPs -r $rho $length -SAA 200 -SaA 100 -SF 1e-2 -N 100000 -Sp 0.5 >output.ms

# calculated statistics along chromosomes
python3 msmscalc_onePop.py output.ms $nIndividuals 1 $length $window_width $window_step $nReplicates outputStats
```

