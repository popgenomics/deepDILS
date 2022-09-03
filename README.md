# deepDILS
```
# simulated 100 replicates of a sweep
msms 40 100 -s 2000 -r 300 100000 -SAA 200 -SaA 100 -SF 1e-2 -N 100000 -Sp 0.5 >100LD.ms

# calculated statistics along chromosomes
python3 msmscalc_onePop.py 100LD.ms 40 1 100000 0.01 0.005 100 output1
```

