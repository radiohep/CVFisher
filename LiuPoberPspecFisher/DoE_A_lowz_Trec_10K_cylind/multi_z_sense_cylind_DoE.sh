#!/bin/bash

zStart=6.0
deltaZ=-0.5
numZs=12
restFreq=1.420 # Rest frequency of spectral line 
fgOption='opt'
fgAmp=2000.0 # in K, at 100 MHz
fgOffset=2.7 # in K
fgSpecIndex=2.4
zHalfWidth=`echo "scale=20; $deltaZ / 2.0" | bc -l`
base_bg_path='/Users/acliu/Research/radiohep/CVFisher/assumptions/background_zlow.dat'
base_pspec_path='/Users/acliu/Research/radiohep/CVFisher/assumptions/pk0.dat'
dish_size=10.0 # in meters
Trx='1e4' # receiver temp in mK
sqSidelength=32 # number of dishes on one side of the square
sqSpacing=10.0 # distance between dish centers in meters
sqSpacing_in_cm=`echo "scale=5; $sqSpacing * 100.0" | bc -l`
dish_size_in_lambda=`echo "scale=5; $dish_size / 2.0" | bc -l` # at 150 MHz
base_cal_fname='base_sq_array.py'
calFname='DoE_A'

# Make cal file for the array that we want to simulate
cat $base_cal_fname | sed s/"    'dish_size_in_lambda': ".*/"    'dish_size_in_lambda': $dish_size_in_lambda,"/ > asdf.py
cat asdf.py | sed s/"    'Trx': ".*/"    'Trx': $Trx"/ > asdf2.py
cat asdf2.py | sed s/"sqSidelength = ".*/"sqSidelength = $sqSidelength"/ > asdf3.py
cat asdf3.py | sed s/"dL = 1000.".*/"dL = $sqSpacing_in_cm"/ > ${calFname}.py
rm asdf.py asdf2.py asdf3.py


for (( i=0 ; $i<$numZs ; i=$i+1 )) ; do
    currentZ=`echo "scale=20; $zStart + $i * $deltaZ" | bc -l`
    currentFreq=`echo "scale=30; 1.420 / ( 1 + $currentZ )" | bc -l`

    upperZ=`echo "scale=20; $currentZ - $zHalfWidth" | bc -l`
    lowerZ=`echo "scale=20; $currentZ + $zHalfWidth" | bc -l`

    upperFreq=`echo "scale=30; 1.420 / ( 1 + $lowerZ )" | bc -l`
    lowerFreq=`echo "scale=30; 1.420 / ( 1 + $upperZ )" | bc -l`
    currentBandwidth=`echo "scale=20; $upperFreq - $lowerFreq" | bc -l`

    # echo $currentFreq, $currentZ, $upperZ, $lowerZ
    # echo $upperFreq, $lowerFreq, $currentBandwidth

    currentTsky=`echo "scale=30; $fgOffset + $fgAmp * e(-1.0 * $fgSpecIndex * l(10.0 * $currentFreq))" | bc -l`
    echo $currentFreq
    echo $currentTsky
    echo "----"

#   echo $currentZ

    # currentZ="`printf '%06.2f' $currentZ`"
    # #printf "%05d\n" $i
    # echo $currentZ,$currentFreq

    python mk_array_file_z.py -C $calFname -f $currentFreq
    current_array_fname="`cat temp.log | tail -1`"
    rm temp.log

    python scalePspec_noBias.py --pspecPath="$base_pspec_path" --bgPath="$base_bg_path" --chosenFreq=$currentFreq --restFreq=$restFreq --outPath="temp_Pk.npz"


    python calc_sense_z_cylind.py -m $fgOption -T $currentTsky -R $restFreq -f $currentFreq --bwidth=$currentBandwidth $current_array_fname --eor="temp_Pk.npz"
    rm temp_Pk.npz

done

### To do:
### ---Change input pspec --- DONE
### ---Change input sky temperature --- DONE
### ---Get Trx and dish_size_in_lambda to be inputted from bash script --- DONE
