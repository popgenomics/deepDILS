#!/bin/bash
module unload python
module load java-jdk/11.0.9.1
module load python/3.9

# number of iterations
nIterations=10000

# Vérifie s'il y a exactement deux arguments
if [ "$#" -ne 2 ]; then
	echo -e "\nUsage: $0 <model_value> <m_value>"
	echo -e "\nmodel equal to CST, EXP or BTL"
	echo "m is 4.N.m"
	exit 1
fi

# Assignation des arguments aux variables 'model' et 'm'
model=$1
m=$2

# example 
binpath=/shared/home/croux/softwares/deepDILS/scripts/dev

# number of sequenced gametes
n=40

# chromosome length
L=100000

# number of polymorphic sites
S=2000

# get the model_name (BTL / CST / EXP / MGB / MIG / MGX)
if [[ $model == "BTL" ]]; then
    if [[ $(awk 'BEGIN { print ('"$m"' == 0 || '"$m"' == 0.0) }') -ne 0 ]]; then
        model_name="BTL"
    else
        model_name="MGB"
    fi
elif [[ $model == "CST" ]]; then
    if [[ $(awk 'BEGIN { print ('"$m"' == 0 || '"$m"' == 0.0) }') -ne 0 ]]; then
        model_name="CST"
    else
        model_name="MIG"
    fi
elif [[ $model == "EXP" ]]; then
    if [[ $(awk 'BEGIN { print ('"$m"' == 0 || '"$m"' == 0.0) }') -ne 0 ]]; then
        model_name="EXP"
    else
        model_name="MGX"
    fi
else
    echo "Aucune correspondance trouvée"
fi

timeStamp=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1)
number_date=$(date '+%d%m%y')

folderName=${model_name}-${number_date}-${timeStamp}

mkdir ${folderName}
rootName=${model_name}-${number_date}-${timeStamp}
outfile=${folderName}/${rootName}

#python3 ${binpath}/launch_pipeline_msms.py --outfile ${outfile} --model ${model} -n ${n} -S ${S} -L ${L} -r 2e-6 -fA 0.01 -m ${m}


# produce the bash script to submit to slurm using sbatch:
bashScript=${model_name}-${number_date}-${timeStamp}.sh
echo "#!/bin/bash" > ${bashScript}
echo "#SBATCH --array=1-${nIterations}%25  # Définit l'array de 1 à ${nIterations} (paquets de 25 jobs)" >> ${bashScript}
echo "#SBATCH --job-name=${model_name}  # Nom du job" >> ${bashScript}
echo -e "#SBATCH --ntasks=1  # Nombre de tâches pour chaque job\n" >> ${bashScript}
#echo "#SBATCH --partition=grey_zone_in_green_world" >> ${bashScript}

echo "# Génère un identifiant unique pour chaque tâche à partir de l'ARRAY" >> ${bashScript}
echo "# Utilise \$SLURM_ARRAY_TASK_ID pour différencier les tâches" >> ${bashScript}
echo -e "task_id=\$((\$SLURM_ARRAY_TASK_ID))  # Exemple : 101 à 110 pour différencier\n" >> ${bashScript}

echo "# binpath" >> ${bashScript}
echo -e "binpath=${binpath}\n" >> ${bashScript}

echo "# number of sequenced gametes" >> ${bashScript}
echo -e "n=${n}\n" >> ${bashScript}

echo "# chromosome length" >> ${bashScript}
echo -e "L=${L}\n" >> ${bashScript}

echo "# number of polymorphic sites" >> ${bashScript}
echo -e "S=${S}\n" >> ${bashScript}

echo "# model" >> ${bashScript}
echo -e "model=${model}\n" >> ${bashScript}

echo "# migration" >> ${bashScript}
echo -e "m=${m}\n" >> ${bashScript}

echo "# outfile" >> ${bashScript}
echo -e "outfile=${outfile}\n" >> ${bashScript}

echo "# Exécute le job en utilisant l'identifiant unique" >> ${bashScript}
echo -e "python3 ${binpath}/launch_pipeline_msms.py --outfile \${outfile}-\${task_id} --model \${model} -n \${n} -S \${S} -L \${L} -r 1e-6 -fA 0.01 -m \${m}\n" >> ${bashScript}

echo -e "\n\tSubmission of :\n\t\tsbatch ${bashScript}"

sbatch ${bashScript}

