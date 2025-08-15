#!/bin/bash

# A function for converting a set of ISTHMUS .sdat/.spar files into BIDS-compliant NIfTI files for each constituent part (short-TE / HERCULES / water-suo / ref).
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
FileLoc=$1  # Path to the data file (assumes matching .sdat & .spar)
SubID=$2    # Subject ID (e.g. 001 for sub-001)
SesID=$3    # Session ID (e.g. 01 for ses-01)
VOI=$4      # Location of VOI (e.g. pcc for voi-pcc)
OutputDIR="$5"/sub-"$SubID"/ses-"$SesID"/mrs # Driectory to save the data (i.e. the BIDS data directory)
#########################################################################################################

# For Philips ".sdat" format, we need the locations of the sdat and spar. Assume the same filename for now, and "_act" / "_ref" naming convention
SparLoc="${FileLoc%.*}.spar"
FileLoc_ref="${FileLoc%_act.*}_ref.sdat"
SparLoc_ref="${FileLoc_ref%.*}.spar"
mkdir -p $OutputDIR


Staging="$OutputDIR"/temp/
mkdir -p $Staging

# act temp files:
CMD="spec2nii philips "$FileLoc" "$SparLoc
CMD="$CMD -o $Staging -f Staging_act"
eval $CMD

# ref temp files:
CMD="spec2nii philips "$FileLoc_ref" "$SparLoc_ref
CMD="$CMD -o $Staging -f Staging_ref"
eval $CMD

