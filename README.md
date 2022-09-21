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
python3 ../scripts/msmscalc_onePop.py infile=output.ms nIndiv=$nIndividuals nCombParam=1 regionSize=$length width=$window_width step=$window_step nRep=$nReplicates root=outputStats
```
Here, a sweep was simulated 100 times in the middle of a chromosome of 100-kb, sequenced in 40 sampled gametes.  
![Alt text](pictures/simulated_sweep.png "simulated sweep")


With two simulated datasets using the Guillaume Lan-Fong's SLiM pipeline producing the following files:
- 101_neutral_positions.txt  
- 101_neutral_sumStats.txt  
- 100_neutral_positions.txt  
- 100_neutral_sumStats.txt  
- 101_neutral.ms  
- 101_neutral_recap_mut.trees  
- 101_neutral_recap.trees  
- 100_neutral.ms  
- 100_neutral_recap_mut.trees  
- 100_neutral_recap.trees  
- 101_neutral_parameters.txt  
- 101_neutral.trees  
- 100_neutral_parameters.txt  
- 100_neutral.trees  
- 101_sweep_positions.txt  
- 101_sweep_sumStats.txt  
- 100_sweep_positions.txt  
- 100_sweep_sumStats.txt  
- 101_sweep.ms  
- 101_sweep_recap_mut.trees  
- 101_sweep_recap.trees  
- 100_sweep.ms  
- 100_sweep_recap_mut.trees  
- 100_sweep_recap.trees  
- 100_sweep_parameters.txt  
- 100_sweep.trees  
- 101_sweep_parameters.txt  
- 101_sweep.trees  

## calculating statistics in sliding windows  
In the _example_ subdirectory:  
```
cd example
for iteration in 100 101; do
	for model in neutral sweep; do
		echo ${iteration}_${model}
		python3 ../scripts/msmscalc_onePop.py infile=${iteration}_${model}.ms nIndiv=40 nCombParam=1 regionSize=100000 width=0.01 step=0.005 nRep=1 root=${iteration}_${model}
	done
done
```
  
## producing jpg for learning  
In the _example_ subdirectory:  
```
for iteration in 100 101; do
	python3 ../scripts/sim2box_single_YOLOv5.py dpi=300 datapath=$PWD simulation=${iteration} object=posSelection theta=0 phasing=1 plotStats=0
done
```
For phased data (phasing=1): statistics relative to LD are computed (nHaplotypes, H1, H2, H12, H2 over H1, D, r2).   
For only plotting jpg files used for machine learning (and not one additional plot per individual statistics): plotStats=0. If you want multiple jpg files: plotStats=1.  
datapath: datapath of a directory with all *neutral* and *sweep* pairs of simulations.  
dpi: resolution of the jpg files.  
theta: specify the way we define the bounding box. If theta=1, then the bounding box is delimited by 4.N.mu. If theta=0, then the bounding box is delimited by the average pi calculated from the *neutral* simulation.  
  
## main outputs  
- 100_neutral_rawData.txt
- 100_neutral_rawData.jpg
- 100_sweep_rawData.txt
- 100_sweep_rawData.jpg
- 101_neutral_rawData.txt
- 101_neutral_rawData.jpg
- 101_sweep_rawData.txt
- 101_sweep_rawData.jpg
- 100_neutral_globalPic.txt
- 100_neutral_globalPic.jpg
- 100_sweep_globalPic.txt
- 100_sweep_globalPic.jpg
- 101_neutral_globalPic.txt
- 101_neutral_globalPic.jpg
- 101_sweep_globalPic.txt
- 101_sweep_globalPic.jpg
  
jpg files correspond to simulated genomes.  
txt files are associated to jpg files and correspond to coordinates of the box for YOLOv5 : object x_center y_center width height.  
globalPic: summary statistics along windows.  
rawData: simulated haplotypes.  

