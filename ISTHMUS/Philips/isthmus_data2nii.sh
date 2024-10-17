#!/bin/bash

# A function for converting a set of ISTHMUS .data/.list/.spar files into BIDS-compliant NIfTI files for each constituent part (short-TE / HERCULES / water-suo / ref).
#
# C.W. Davies-Jenkins, Johns Hopkins University 2024
#
# DEPENDENCIES:
# Python packages
# spec2nii (tested on v0.8.2)
#       https://github.com/wtclarke/spec2nii
# mrs_tools (tested on v0.1.7)
#       https://github.com/wtclarke/nifti_mrs_tools

# 
# Inputs: #############################################################################################
FileLoc=$1  # Path to the data file (assumes matching .list & .spar)
SubID=$2    # Subject ID (e.g. 001 for sub-001)
SesID=$3    # Session ID (e.g. 01 for ses-01)
VOI=$4      # Location of VOI (e.g. pcc for voi-pcc)
OutputDIR="$5"/sub-"$SubID"/ses-"$SesID"/mrs # Driectory to save the data (i.e. the BIDS data directory)
#########################################################################################################

# For Philips ".data" format, we need the locations of the list and spar. Assume the same filename for now, but we can put them as additional arguments.
ListLoc="${FileLoc%.*}.list"
SparLoc="${FileLoc%.*}.spar"

mkdir -p $OutputDIR

# Directory for temporary files:
Staging="$OutputDIR"/temp/
mkdir -p $Staging

# For data/list:
CMD="spec2nii philips_dl "$FileLoc" "$ListLoc" "$SparLoc
Ext=".data"

# Run initial spec2nii conversion
CMD="$CMD -o $Staging -f Staging"
eval $CMD

# Split the water-suppressed file and push it to the output directory:
MetFile="$Staging"Staging_STD_0.nii.gz
eval "mrs_tools split --file $MetFile --dim DIM_USER_0 --indices 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 --output $Staging"

# Move water-suppressed short-TE ############################################
eval spec2nii anon "$Staging"Staging_STD_0_selected.nii.gz 
eval spec2nii extract "$Staging"Staging_STD_0_selected.nii.gz 

# Correct JSON sidecar
PYCMD=$(cat <<EOF
import json
jsonHeaderFile = open("${Staging}Staging_STD_0_selected.json")
HeaderFileData = json.load(jsonHeaderFile)
jsonHeaderFile.close()
HeaderFileData["EchoTime"] = 0.035
HeaderFileData["WaterSuppressed"] = True
HeaderFileData["Manufacturer"] = "Philips"
HeaderFileData["dim_6"] = "DIM_DYN"
file = open("${Staging}Staging_STD_0_selected.json", 'w')
json.dump(HeaderFileData, file, indent=4)
file.close()
EOF
)
python3 -c "$PYCMD"
eval spec2nii insert "$Staging"Staging_STD_0_selected.nii.gz "$Staging"Staging_STD_0_selected.json -o "$Staging" -f Staging_STD_0_selected


cp "$Staging"Staging_STD_0_selected.nii.gz "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-press_te-35_svs.nii.gz
cp "$Staging"Staging_STD_0_selected.json "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-press_te-35_svs.json

# Move water-suppressed HERCULES ############################################
mrs_tools reshape --file "$Staging"Staging_STD_0_others.nii.gz --shape -1 56 4 --d6 DIM_DYN --d7 DIM_EDIT --output "$Staging" --filename Staging_STD_0_others.nii.gz 
eval spec2nii anon "$Staging"Staging_STD_0_others.nii.gz 
eval spec2nii extract "$Staging"Staging_STD_0_others.nii.gz

# Correct JSON sidecar
PYCMD=$(cat <<EOF
import json
jsonHeaderFile = open("${Staging}Staging_STD_0_others.json")
HeaderFileData = json.load(jsonHeaderFile)
jsonHeaderFile.close()
HeaderFileData["EchoTime"] = 0.080
HeaderFileData["WaterSuppressed"] = True
HeaderFileData["Manufacturer"] = "Philips"
file = open("${Staging}Staging_STD_0_others.json", 'w')
json.dump(HeaderFileData, file, indent=4)
file.close()
EOF
)
python3 -c "$PYCMD"
eval spec2nii insert "$Staging"Staging_STD_0_others.nii.gz "$Staging"Staging_STD_0_others.json -o "$Staging" -f Staging_STD_0_others

cp "$Staging"Staging_STD_0_others.nii.gz "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-hercules_te-80_svs.nii.gz
cp "$Staging"Staging_STD_0_others.json "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-hercules_te-80_svs.json

# Split the water reference file and push it to the output directory:
RefFile="$Staging"Staging_STD_1.nii.gz
eval "mrs_tools split --file $RefFile --dim DIM_USER_0 --indices 0 2 4 6 --output $Staging"

# Move water-reference short-TE ############################################
eval spec2nii anon "$Staging"Staging_STD_1_selected.nii.gz 
eval spec2nii extract "$Staging"Staging_STD_1_selected.nii.gz 

# Correct JSON sidecar
PYCMD=$(cat <<EOF
import json
jsonHeaderFile = open("${Staging}Staging_STD_1_selected.json")
HeaderFileData = json.load(jsonHeaderFile)
jsonHeaderFile.close()
HeaderFileData["EchoTime"] = 0.035
HeaderFileData["WaterSuppressed"] = False
HeaderFileData["Manufacturer"] = "Philips"
HeaderFileData["dim_6"] = "DIM_DYN"
file = open("${Staging}Staging_STD_1_selected.json", 'w')
json.dump(HeaderFileData, file, indent=4)
file.close()
EOF
)
python3 -c "$PYCMD"
eval spec2nii insert "$Staging"Staging_STD_1_selected.nii.gz "$Staging"Staging_STD_1_selected.json -o "$Staging" -f Staging_STD_1_selected

cp "$Staging"Staging_STD_1_selected.nii.gz "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-press_te-35_mrsref.nii.gz
cp "$Staging"Staging_STD_1_selected.json "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-press_te-35_mrsref.json

# Move water-reference HERCULES ############################################
mrs_tools reshape --file "$Staging"Staging_STD_1_others.nii.gz --shape -1 4 --d6 DIM_EDIT --output "$Staging" --filename Staging_STD_1_others.nii.gz 
eval spec2nii anon "$Staging"Staging_STD_1_others.nii.gz 
eval spec2nii extract "$Staging"Staging_STD_1_others.nii.gz 

# Correct JSON sidecar
PYCMD=$(cat <<EOF
import json
jsonHeaderFile = open("${Staging}Staging_STD_1_others.json")
HeaderFileData = json.load(jsonHeaderFile)
jsonHeaderFile.close()
HeaderFileData["EchoTime"] = 0.080
HeaderFileData["WaterSuppressed"] = False
HeaderFileData["Manufacturer"] = "Philips"
file = open("${Staging}Staging_STD_1_others.json", 'w')
json.dump(HeaderFileData, file, indent=4)
file.close()
EOF
)
python3 -c "$PYCMD"
eval spec2nii insert "$Staging"Staging_STD_1_others.nii.gz "$Staging"Staging_STD_1_others.json -o "$Staging" -f Staging_STD_1_others

cp "$Staging"Staging_STD_1_others.nii.gz "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-press_te-80_mrsref.nii.gz
cp "$Staging"Staging_STD_1_others.json "$OutputDIR"/sub-"$SubID"_ses-"$SesID"_voi-"$VOI"_acq-press_te-80_mrsref.json

# Remove temporary directory at the end
rm -r "$Staging"

echo "Conversion complete!"