#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -l walltime=744:00:00

cat $PBS_NODEFILE > hostlist
NPROCS=`wc -l < $PBS_NODEFILE`

cd /home/sangiova/prog/weak-coupling/check_moment_bethe/mu-0.33

export DIR=/home/sangiova/prog/weak-coupling

declare -i niter
declare -i n
declare -i nm1
let niter=20

if [ !  -e Iter_1 ]; then
mkdir Iter_1
cd Iter_1
cp ../Start/* .
nice -19 $DIR/Self/self.out
nice -19 $DIR/CTQMC/CTQMC.out
nice -19 $DIR/Analysis/jackv4.out
nice -19 $DIR/Self/self.out
cd ..
fi

let n=2
while [ $n -le  $niter ];
do
    
    let nm1=n-1
    export NewDir="Iter_"$n
    export OldDir="Iter_"$nm1
    if [ ! -e  $NewDir ]; then
       mkdir $NewDir
       cd  $NewDir
       cp ../$OldDir/* .
       bash restart.c
       nice -19 $DIR/CTQMC/CTQMC.out
       nice -19 $DIR/Analysis/jackv4.out
       nice -19 $DIR/Self/self.out
       cd ..
    fi
    let n=n+1
done 


