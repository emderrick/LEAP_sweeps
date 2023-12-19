
Annotation of my candidate MAGs with bakta since I had previously done it with prokka. I installed bakta in a virtual environment in the venvs directory in home. This is bakta version 1.8.1.

```bash
module load python/3.10
virtualenv --no-download bakta
source bakta/bin/activate
pip install --no-index --upgrade pip
pip install bakta

wget www.drive5.com/pilercr/pilercr1.06.tar.gz # install pilercr
tar -zxvf pilercr1.06.tar.gz
cd pilercr1.06
make
cd ..

git clone https://github.com/ncbi/amr.git
cd amr
git checkout master
make
make clean
make DEFAULT_DB_DIR=/home/ederrick/scratch/dbs/bakta/amrfinderplus-db/
make install INSTALL_DIR=/home/ederrick/venvs/
```

Then I downloaded the database with

```bash
wget https://zenodo.org/record/7669534
tar -xzf db.tar.gz
rm db.tar.gz
```

Then I updated the database to include the amrfinder db with bakta's internal command.
```bash
export PATH=/home/ederrick/venvs/amr:/home/ederrick/venvs/pilercr1.06:$PATH
amrfinder_update --force_update --database db/amrfinderplus-db
```

Then I ran bakta on each MAG.
```bash
#!/usr/bin/bash
#SBATCH --account=
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=16

source /home/ederrick/venvs/bakta/bin/activate

module load StdEnv/2020  gcc/9.3.0 trnascan-se/2.0.12
module load aragorn/1.2.41
module load infernal/1.1.4
module load hmmer/3.3.2
module load diamond/2.0.15
module load blast+/2.13.0
module load circos/0.69-9

export PATH=/home/ederrick/venvs/amr:/home/ederrick/venvs/pilercr1.06:$PATH

for file in candidate_MAGs/*.fa
do
out="${file//.fa/_full_bakta_output}"
bakta $file --db /home/ederrick/scratch/dbs/bakta  --output $out --threads 16
done
deactivate
```
