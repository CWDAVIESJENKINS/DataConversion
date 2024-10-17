function[] = isthmus_sdat2nii(InFile, SubID, SesID, VOI, OutFolder)
%% function[] = isthmus_sdat2nii(InFile, SubID, SesID, VOI, OutFolder)
%
% Description: Converts a pair of ISTHMUS sdat files (act & ref) into the 4
% NIfTI files that can be used for Osprey processing.
%
% Input:     InFile = location of the metabolite ("act") scan
%            SubID = subject identifier ['text']
%            SesID = session identifier ['text']
%            VOI = location of the voxel ['text']
%            OutFolder = Location of data output folder
%
% Example usage: 
% SplitIsthmus_sdat('data_act.sdat', 'JHU01', '01', 'pcc', 'project/data_raw')
%
% C.W. Davies-Jenkins, Johns Hopkins University 2024

OspVer=getCurrentVersion;
OspVer = OspVer.Version;



OutLoc = fullfile(OutFolder,['sub-',SubID],['ses-',SesID],'mrs');
if ~exist(OutLoc,"dir")
    mkdir(OutLoc);
end

Spec = io_loadspec_sdat(InFile,1);

%% Export short-TE metabolite data
Spec_ste = op_takeaverages(Spec,1:32);
Spec_ste.te = 35; 
Spec_ste.seq = 'PRESS';
OutLoc = fullfile(OutFolder,['sub-',SubID],['ses-',SesID],'mrs',['sub-',SubID,'_ses-',SesID,'_voi-',VOI,'_acq-press_te-35_svs.nii.gz']);
io_writeniimrs(Spec_ste, OutLoc, {}, OspVer);
%% Export HERCULES metabolite data
Spec_Herc = op_takeaverages(Spec,33:Spec.averages);
Spec_Herc.fids = reshape(Spec_Herc.fids,[2048,4,224/4]);
Spec_Herc.fids = permute(Spec_Herc.fids,[1,3,2]);
Spec_Herc.specs = reshape(Spec_Herc.specs,[2048,4,224/4]);
Spec_Herc.specs = permute(Spec_Herc.specs,[1,3,2]);
Spec_Herc.sz = size(Spec_Herc.specs);
Spec_Herc.seq = 'HERCULES';
Spec_Herc.dims.subSpecs=3;
Spec_Herc.averages = 56;
Spec_Herc.subspecs = 4;
OutLoc = fullfile(OutFolder,['sub-',SubID],['ses-',SesID],'mrs',['sub-',SubID,'_ses-',SesID,'_voi-',VOI,'_acq-hercules_te-80_svs.nii.gz']);
io_writeniimrs(Spec_Herc, OutLoc, {'DIM_EDIT'}, OspVer);

%% Try and find ref data
InFile = strrep(InFile,'_act.','_ref.');

if exist(InFile,"file")
    Spec = io_loadspec_sdat(InFile,1);

    Spec_ste = op_takeaverages(Spec,1:32);
    Spec_ste.te = 35; 
    Spec_ste.seq = 'PRESS';
    OutLoc = fullfile(OutFolder,['sub-',SubID],['ses-',SesID],'mrs',['sub-',SubID,'_ses-',SesID,'_voi-',VOI,'_acq-press_te-35_mrsref.nii.gz']);
    io_writeniimrs(Spec_ste, OutLoc, {}, OspVer);

    Spec_Herc = op_takeaverages(Spec,33:Spec.averages);
    Spec_Herc.fids = reshape(Spec_Herc.fids,[2048,4,224/4]);
    Spec_Herc.fids = permute(Spec_Herc.fids,[1,3,2]);
    Spec_Herc.specs = reshape(Spec_Herc.specs,[2048,4,224/4]);
    Spec_Herc.specs = permute(Spec_Herc.specs,[1,3,2]);
    Spec_Herc.sz = size(Spec_Herc.specs);
    Spec_Herc.seq = 'HERCULES';
    Spec_Herc.dims.subSpecs=3;
    Spec_Herc.averages = 56;
    Spec_Herc.subspecs = 4;
    OutLoc = fullfile(OutFolder,['sub-',SubID],['ses-',SesID],'mrs',['sub-',SubID,'_ses-',SesID,'_voi-',VOI,'_acq-hercules_te-80_mrsref.nii.gz']);
    io_writeniimrs(Spec_Herc, OutLoc, {'DIM_EDIT'}, OspVer);
end



end
