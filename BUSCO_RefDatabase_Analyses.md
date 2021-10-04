/data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref

/data/putnamlab/erin_chille/mcap2019/data/symbiont/Trinity_psymTrans/1_Trinity/config.ini


Adig.fasta
Ahya.fasta
Amil.fasta
Apal.fasta
Apoc.fasta
Aque.fasta
Emue.fasta
Gasp.fasta
Gfas.fasta
Mcap.fasta
Mlei.fasta
Nvec.fasta
Pbac.fasta
Pdam.fasta
Plut.fasta
Prus.fasta
Pver.fasta
Spis.fasta


# BUSCO Config

```
mkdir busco_downloads

wget https://busco-data.ezlab.org/v4/data/lineages/metazoa_odb10.2020-09-10.tar.gz
mv metazoa_odb10.2020-09-10.tar.gz busco_downloads/


wget https://busco-data.ezlab.org/v4/data/file_versions.tsv
mv file_versions.tsv busco_downloads/
```

# A. digitifera BUSCO
```
nano /data/putnamlab/erin_chille/Fellous_Epi_Machinery/scripts/BUSCO.sh
```

```
#!/bin/bash
#SBATCH --job-name="BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref
#SBATCH --mem=50GB
#SBATCH -q putnamlab

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4

#run BUSCO
busco --config config.ini -m transcriptome -i Adig.fasta -o Adig_BUSCO -l /data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref/busco_downloads/metazoa_odb10 --offline

echo "BUSCO Mission complete!" $(date)
```
```
sbatch /data/putnamlab/erin_chille/Fellous_Epi_Machinery/scripts/BUSCO.sh
```

# A. hyacinthus BUSCO Test Run 
```
nano /data/putnamlab/erin_chille/Fellous_Epi_Machinery/scripts/BUSCO.sh
```

```
#!/bin/bash
#SBATCH --job-name="BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref
#SBATCH --mem=50GB
#SBATCH -q putnamlab

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4

#run BUSCO
busco --config config.ini -m transcriptome -i Ahya.fasta -o Ahya_BUSCO -l /data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref/busco_downloads/metazoa_odb10 --offline

echo "BUSCO Mission complete!" $(date)
```
```
sbatch /data/putnamlab/erin_chille/Fellous_Epi_Machinery/scripts/BUSCO.sh
```

### A. digitifera BUSCO Results
```C:78.6%[S:77.7%,D:0.9%],F:8.3%,M:13.1%,n:954
750    Complete BUSCOs (C)                       
741    Complete and single-copy BUSCOs (S)       
9      Complete and duplicated BUSCOs (D)        
79     Fragmented BUSCOs (F)                     
125    Missing BUSCOs (M)                        
954    Total BUSCO groups searched
```


```
nano /data/putnamlab/erin_chille/Fellous_Epi_Machinery/scripts/Loop_BUSCO.sh
```

```
#!/bin/bash
#SBATCH --job-name="BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref
#SBATCH --mem=100GB
#SBATCH -q putnamlab

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4

#run BUSCO

for i in "Adig" "Ahya" "Amil" "Mcap" "Plut" "Prus" "Apoc" "Gasp" "Gfas" "Pver" "Pdam" "Spis" "Nvec" "Apal" "Mlei" "Pbac" "Aque" "Emue"; do
busco --config config.ini -m transcriptome -i ${i}.fasta -o ${i}_BUSCO -l /data/putnamlab/erin_chille/Fellous_Epi_Machinery/ref/busco_downloads/metazoa_odb10 --offline
done
```

```
sbatch /data/putnamlab/erin_chille/Fellous_Epi_Machinery/scripts/Loop_BUSCO.sh
```