#!/bin/bash

## This is a script to run novaCTF in completeness to generate a bin1 WBP tomogram
## and a bin4 SIRT-like filtered tomogram. It goes through each step of
## novaCTF, running one tomogram at a time through to completeness to avoid
## over-filling the hard drive storage with intermediate files. It is parallelized
## but without memory allocation/limits, so use with caution. Replace references of
## <ctfpath> with the path to your executable.

## Inputs: IMOD .rawtlt and .com files from alignment, ctffind4 or IMOD
## defocus determination file labeled as defocus_corrected.txt. Amplitude contrast,
## Cs and voltage should match ctffind4 inputs, and goldradius the radius used for 
## erasing of gold in etomo. IMOD defocus format is assumed; this can 
## be changed in the defocus.com file generation step.  Directory names are assumed
## to be two-digit integers, with one folder for each tomogram and the tomogram 
## name ##.mrc matching the directory name ##.

## This script is heavily based on guidance for running novaCTF, written by Beata
## Turonova and distributed at https://github.com/turonova/novaCTF/wiki. This
## wrapper was written by Lauren Ann Metskas on 1 November 2018.

## This work is licensed under a GPL v.3 license, in compliance with the license
## of the parent script and guidance at https://github.com/turonova/novaCTF/wiki.




## INPUTS
export maindir='<filename>'
export pixsize=<>    # pixel size in nm
defstep=15		# defocus step size to search (in nm)
nCores=1		# number of cores to run locally in parallel (CPU)
export goldradius=<>    # radius to erase for gold beads.  Suggested: 61 for 1.104 A pix size.



## OTHER PARAMETERS
amplitudeContrast=<>
Cs=<>
Voltage=<> 
tomos="12 14" #change this list definition to automate according to your data set-up
export interpolation=NearestNeighbor 	# imod interpolation desired
astig=1			# correct for astigmatism if information is present in file

################## CODE BELOW, DO NOT CHANGE BELOW THIS LINE ##################

a=1
for t in $tomos
do 
echo Starting on tomogram $t

if [ ${#t} -le 1 ]; then
	workingdir=$maindir/0$t
	a=0$t
else
	workingdir=$maindir/$t
	a=$t
fi

cd $workingdir/

header ${a}.mrc > temp.txt
export tomosizey=`more temp.txt | grep sections | cut -c45-48`
export tomosizex=`more temp.txt | grep sections | cut -c53-56`
export thickness=`more $workingdir/tilt.com | grep THICKNESS | cut -c11-16`

cd $workingdir/
newstack -input ${a}_dose-filt.st -output ${a}.ali -xform *_fid.xf -size $tomosizex,$tomosizey -origin -$interpolation

####################### ROUND 1: NOVA CTF DEFOCUS #############################

echo "Algorithm defocus" > setup_defocus.com
echo "InputProjections ${a}.ali" >> setup_defocus.com
echo "FULLIMAGE ${tomosizex} ${tomosizey}" >> setup_defocus.com 
echo "THICKNESS ${thickness}" >> setup_defocus.com
echo "TILTFILE ${a}.rawtlt" >> setup_defocus.com
echo "SHIFT 0.0 0.0" >> setup_defocus.com
echo "CorrectionType phaseflip" >> setup_defocus.com
echo "DefocusFileFormat ctffind4" >> setup_defocus.com
echo "DefocusFile defocus_corrected.txt" >> setup_defocus.com
echo "PixelSize ${pixsize}" >> setup_defocus.com
echo "DefocusStep ${defstep}" >> setup_defocus.com
echo "CorrectAstigmatism ${astig}" >> setup_defocus.com
if [ -e extra_defocus.txt ]; then
        echo "DefocusShiftFile extra_defocus.txt" >> setup_defocus.com
fi

<ctfpath> -param setup_defocus.com

####################### ROUND 2: NOVA CTF CORRECTION ##########################

nslices=`ls defocus_corrected.txt_* | wc -l`
lastSlice="$((nslices-1))"


for ((n=0;n<=lastSlice;n++))
do

echo "Algorithm ctfCorrection" > setup_ctfCorrection_$n.com
echo "InputProjections ${a}.ali" >> setup_ctfCorrection_$n.com
echo "DefocusFile defocus_corrected.txt_${n}" >> setup_ctfCorrection_$n.com
echo "OutputFile corrected_stack.st_${n}" >> setup_ctfCorrection_$n.com
echo "TILTFILE ${a}.rawtlt" >> setup_ctfCorrection_$n.com
echo "CorrectionType phaseflip" >> setup_ctfCorrection_$n.com
echo "DefocusFileFormat ctffind4" >> setup_ctfCorrection_$n.com
echo "PixelSize ${pixsize}" >> setup_ctfCorrection_$n.com
echo "AmplitudeContrast ${amplitudeContrast}" >> setup_ctfCorrection_$n.com
echo "Cs ${Cs}" >> setup_ctfCorrection_$n.com
echo "Volt ${Voltage}" >> setup_ctfCorrection_$n.com
echo "CorrectAstigmatism ${astig}" >> setup_ctfCorrection_$n.com

done

(
for ((k=0;k<=lastSlice;k++))
do
((i=i%nCores)); ((i++==0)) && wait   
<ctfpath> -param setup_ctfCorrection_$k.com &
done
)

####################### ROUND 3: NOVA FILTER AND IMOD PROCESSING #############

for ((l=0;l<=lastSlice;l++))
do

ccderaser -input corrected_stack.st_$l -output erased_stack.ali_$l -ModelFile ${a}_dose-filt_erase.fid -BetterRadius $goldradius -PolynomialOrder 0 -MergePatches -ExcludeAdjacent -CircleObjects /
echo "Algorithm filterProjections" > setup_filter_$l.com
echo "InputProjections corrected_stack_flipped.ali_${l}" >> setup_filter_$l.com
echo "OutputFile filtered_stack.ali_${l}" >> setup_filter_$l.com
echo "TILTFILE ${a}.rawtlt" >> setup_filter_$l.com
echo "StackOrientation xz" >> setup_filter_$l.com

done


(
for ((j=0;j<=lastSlice;j++))
do   
((i=i%nCores)); ((i++==0)) && wait

mrctaper -t 100 erased_stack.ali_$j
clip flipyz erased_stack.ali_$j corrected_stack_flipped.ali_$j
<ctfpath> -param setup_filter_$j.com

done
)

####################### ROUND 4: NOVA RECONSTRUCTION ######################

shift=`more $workingdir/tilt.com | grep SHIFT | cut -c11-${cutposition}`

echo "Algorithm 3dctf" > setup_reconstruction.com
echo "InputProjections filtered_stack.ali" >> setup_reconstruction.com
echo "OutputFile ${a}.rec" >> setup_reconstruction.com
echo "TILTFILE ${a}.rawtlt" >> setup_reconstruction.com
echo "THICKNESS ${thickness}" >> setup_reconstruction.com
echo "FULLIMAGE ${tomosizex} ${tomosizey}" >> setup_reconstruction.com
echo "SHIFT 0.0 ${shift}" >> setup_reconstruction.com
echo "PixelSize ${pixsize}" >> setup_reconstruction.com
echo "NumberOfInputStacks $nslices" >> setup_reconstruction.com
echo "Use3DCTF 1" >> setup_reconstruction.com

<ctfpath> -param setup_reconstruction.com

######################## STEP 5: POLISH AND CLEAN UP ####################

trimvol -rx ${a}.rec ${a}_bin1_phaseflip.rec

rm ${a}.rec
rm ${a}_output.txt_*

mv ${a}_bin1_phaseflip.rec ${a}_bin1.rec
binvol -bin 4 -antialias 5 ${a}_bin1.rec ${a}_bin4.rec
binvol -bin 2 -antialias 5 ${a}_bin4.rec ${a}_bin8.rec

####################### ROUND 3: NOVA FILTER AND IMOD PROCESSING #############

for ((l=0;l<=lastSlice;l++))
do

#ccderaser -input corrected_stack.st_$l -output erased_stack.ali_$l -ModelFile ${a}_dose-filt_erase.fid -BetterRadius $goldradius -PolynomialOrder 0 -MergePatches -ExcludeAdjacent -CircleObjects /
echo "Algorithm filterProjections" > setup_filter_$l.com
echo "InputProjections corrected_stack_flipped.ali_${l}" >> setup_filter_$l.com
echo "OutputFile filtered_stack.ali_${l}" >> setup_filter_$l.com
echo "TILTFILE ${a}.rawtlt" >> setup_filter_$l.com
echo "StackOrientation xz" >> setup_filter_$l.com
echo "FakeSIRTiterations 13" >> setup_filter_$l.com

done


(
for ((j=0;j<=lastSlice;j++))
do   
((i=i%nCores)); ((i++==0)) && wait

mrctaper -t 100 erased_stack.ali_$j
clip flipyz erased_stack.ali_$j corrected_stack_flipped.ali_$j
<ctfpath> -param setup_filter_$j.com

done
)

####################### ROUND 4: NOVA RECONSTRUCTION ######################

shift=`more $workingdir/tilt.com | grep SHIFT | cut -c11-${cutposition}`

echo "Algorithm 3dctf" > setup_reconstruction.com
echo "InputProjections filtered_stack.ali" >> setup_reconstruction.com
echo "OutputFile ${a}_SIRT.rec" >> setup_reconstruction.com
echo "TILTFILE ${a}.rawtlt" >> setup_reconstruction.com
echo "THICKNESS ${thickness}" >> setup_reconstruction.com
echo "FULLIMAGE ${tomosizex} ${tomosizey}" >> setup_reconstruction.com
echo "SHIFT 0.0 ${shift}" >> setup_reconstruction.com
echo "PixelSize ${pixsize}" >> setup_reconstruction.com
echo "NumberOfInputStacks $nslices" >> setup_reconstruction.com
echo "Use3DCTF 1" >> setup_reconstruction.com

<ctfpath> -param setup_reconstruction.com

######################## STEP 5: POLISH AND CLEAN UP ####################

trimvol -rx ${a}_SIRT.rec ${a}_SIRT_bin1.rec

rm ${a}_SIRT.rec
rm ${a}_output.txt_*
rm defocus_corrected.txt_*
rm setup_ctfCorrection_*.com
rm filtered_stack.ali_*
rm corrected_stack_flipped.ali_*
rm setup_filter_*.com
rm erased_stack.ali_*
rm corrected_stack.st_*
rm temp.txt

binvol -bin 4 -antialias 5 ${a}_SIRT_bin1.rec ${a}_SIRT_bin4.rec
rm ${a}_SIRT_bin1.rec

cd $maindir/

done

