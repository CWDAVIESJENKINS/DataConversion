function[] = isthmus_twix2nii(InFile, OutFolder)
% Converts ISTHMUS twix files to NIfTI for Osprey processing.
% Acquisition order:
%   1) Long TE water
%   2) 32 × PRESS short TE
%   3) Short TE water
%   4) 32 × HERCULES
%   5) Long TE water
%   6) 32 × HERCULES
%   7) Short TE water
%   8) 32 × HERCULES
%   ... alternating long/short TE water with HERCULES blocks
%
% Inputs:
%   InFile   - Path to twix file
%   OutFolder - Output folder
%
% CW Davies-Jenkins / Updated by Gizeaddis Simegn, 2025
%
% V1: CWDJ initial script
% V2: GS generalized indices
% V3: Fix error in the delineation of HERCULES subspectra

%% Setup
OspVer = getCurrentVersion;
OspVer = OspVer.Version;
if ~exist(OutFolder, "dir")
    mkdir(OutFolder);
end
%% Load in the Twix file
twix_obj = mapVBVD(InFile);
if isstruct(twix_obj)
    disp('Single RAID file detected.');
elseif iscell(twix_obj)
    disp('Multi RAID file detected.');
    twix_obj = twix_obj{end}; % use last RAID block
end
Spec_IMG = TwixObj_2_FIDA(twix_obj, 'image'); Spec_IMG.OriginalFile = InFile;
Spec_RTF = TwixObj_2_FIDA(twix_obj, 'RTfeedback'); Spec_RTF.OriginalFile = InFile;
%% --- Dynamic indexing based on acquisition pattern ---
nTrans = size(Spec_IMG.fids, 3);
% Initialise index lists
longTE_water_idx = [];
shortTE_water_idx = [];
press_idx = [];
hercules_idx = [];
% First block: long TE water + 32 × PRESS
longTE_water_idx = [longTE_water_idx, 1];
press_idx        = [press_idx, 2:33];
% Remaining blocks: alternate short/long TE water + 32 × HERCULES
current_idx = 34;
water_toggle = 0; % 0 = short TE, 1 = long TE
while current_idx <= nTrans
    if water_toggle == 0
        shortTE_water_idx = [shortTE_water_idx, current_idx];
    else
        longTE_water_idx  = [longTE_water_idx, current_idx];
    end
    % HERCULES block
    if current_idx + 32 <= nTrans
        hercules_idx = [hercules_idx, current_idx+1 : current_idx+32];
    else
        hercules_idx = [hercules_idx, current_idx+1 : nTrans];
    end
    current_idx = current_idx + 33;
    water_toggle = 1 - water_toggle; % alternate water TE
end
%% --- Extract from correct datasets ---
% PRESS short TE (from Spec_IMG, excluding water)
SpecOut_STE = op_takeaverages(Spec_IMG, press_idx);
SpecOut_STE.te = 35;
io_writeniimrs(SpecOut_STE, fullfile(OutFolder, 'PRESS.nii.gz'), [], OspVer);
% HERCULES (from Spec_IMG, excluding water)
SpecOut_HERC = op_takeaverages(Spec_IMG, hercules_idx);
%SpecOut_HERC.fids = reshape(SpecOut_HERC.fids, SpecOut_HERC.sz(1), SpecOut_HERC.sz(2), [], 4);
%SpecOut_HERC.specs = reshape(SpecOut_HERC.specs, SpecOut_HERC.sz(1), SpecOut_HERC.sz(2), [], 4);
SpecOut_HERC.fids = permute(reshape(SpecOut_HERC.fids, SpecOut_HERC.sz(1), SpecOut_HERC.sz(2), 4, []), [1,2,4,3]);
SpecOut_HERC.specs = permute(reshape(SpecOut_HERC.specs, SpecOut_HERC.sz(1), SpecOut_HERC.sz(2), 4, []), [1,2,4,3]);

SpecOut_HERC.sz = size(SpecOut_HERC.specs);
SpecOut_HERC.dims.subSpecs = 4;
SpecOut_HERC.subspecs = 4;
SpecOut_HERC.rawSubspecs = 4;
SpecOut_HERC.averages = size(SpecOut_HERC.fids, 3) / 4;
SpecOut_HERC.rawAverages = SpecOut_HERC.averages;
SpecOut_HERC.te = 80;
io_writeniimrs(SpecOut_HERC, fullfile(OutFolder, 'HERCULES.nii.gz'), {'DIM_EDIT'}, OspVer);
% Water long TE (from Spec_RTF only)
water_longTE = op_takeaverages(Spec_RTF, longTE_water_idx);
water_longTE.te = 80;
io_writeniimrs(water_longTE, fullfile(OutFolder, 'Water_longTE.nii.gz'), [], OspVer);
% Water short TE (from Spec_RTF only)
water_shortTE = op_takeaverages(Spec_RTF, shortTE_water_idx);
water_shortTE.te = 35;
io_writeniimrs(water_shortTE, fullfile(OutFolder, 'Water_shortTE.nii.gz'), [], OspVer);
end

function[out] = TwixObj_2_FIDA(twix_obj,Struct)

dOut.data=twix_obj.(Struct)();
version=twix_obj.(Struct).softwareVersion;
sqzSize=twix_obj.(Struct).sqzSize;
sqzDims=twix_obj.(Struct).sqzDims;
leftshift = twix_obj.(Struct).freeParam(1);

%Check if the twix file is from a VE version
if contains(twix_obj.hdr.Dicom.SoftwareVersions, 'E11')
    version='ve';
end
if contains(twix_obj.hdr.Dicom.SoftwareVersions, 'XA')
    index = strfind(twix_obj.hdr.Dicom.SoftwareVersions,'XA');
    version=twix_obj.hdr.Dicom.SoftwareVersions(index:index+3);
end

%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging
%to find the correct sequence name.
sequence=twix_obj.hdr.Config.SequenceFileName;

seq = 'isthmus';

fids=double(squeeze(dOut.data));

sz = size(fids);
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;

%Get TxFrq
if strcmp(version,'ve')
    txfrq=twix_obj.hdr.Meas.lFrequency  ;
else
    try
        txfrq=twix_obj.hdr.Meas.Frequency;
    catch
        txfrq=twix_obj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;
    end
    if isempty(txfrq)
        txfrq=twix_obj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1, 1}.lFrequency;
    end
end

% Extract voxel dimensions
if (strcmp(version,'vd') || strcmp(version,'vb') || contains(version,'XA'))
    TwixHeader.VoI_RoFOV     = twix_obj.hdr.Config.VoI_RoFOV; % Voxel size in readout direction [mm]
    TwixHeader.VoI_PeFOV     = twix_obj.hdr.Config.VoI_PeFOV; % Voxel size in phase encoding direction [mm]
    TwixHeader.VoIThickness  = twix_obj.hdr.Config.VoI_SliceThickness; % Voxel size in slice selection direction [mm]
    TwixHeader.PosCor         = twix_obj.hdr.Config.VoI_Position_Cor; % Coronal coordinate of voxel [mm]
    TwixHeader.PosSag         = twix_obj.hdr.Config.VoI_Position_Sag; % Sagittal coordinate of voxel [mm]
    TwixHeader.PosTra         = twix_obj.hdr.Config.VoI_Position_Tra; % Transversal coordinate of voxel [mm]
    TwixHeader.VoI_InPlaneRot = twix_obj.hdr.Config.VoI_InPlaneRotAngle; % Voxel rotation in plane
    TwixHeader.NormCor        = twix_obj.hdr.Config.VoI_Normal_Cor; % Coronal component of normal vector of voxel
    TwixHeader.NormSag        = twix_obj.hdr.Config.VoI_Normal_Sag; % Sagittal component of normal vector of voxel
    TwixHeader.NormTra        = twix_obj.hdr.Config.VoI_Normal_Tra; % Transversal component of normal vector of voxel
else
    TwixHeader.VoI_RoFOV     = twix_obj.hdr.Spice.VoiReadoutFOV; % Voxel size in readout direction [mm]
    TwixHeader.VoI_PeFOV     = twix_obj.hdr.Spice.VoiPhaseFOV; % Voxel size in phase encoding direction [mm]
    TwixHeader.VoIThickness  = twix_obj.hdr.Spice.VoiThickness; % Voxel size in slice selection direction [mm]
    TwixHeader.PosCor         = twix_obj.hdr.Spice.VoiPositionCor; % Coronal coordinate of voxel [mm]
    TwixHeader.PosSag         = twix_obj.hdr.Spice.VoiPositionSag; % Sagittal coordinate of voxel [mm]
    TwixHeader.PosTra         = twix_obj.hdr.Spice.VoiPositionTra; % Transversal coordinate of voxel [mm]
    TwixHeader.VoI_InPlaneRot = twix_obj.hdr.Spice.VoiInPlaneRot; % Voxel rotation in plane
    TwixHeader.NormCor        = twix_obj.hdr.Spice.VoiNormalCor; % Coronal component of normal vector of voxel
    TwixHeader.NormSag        = twix_obj.hdr.Spice.VoiNormalSag; % Sagittal component of normal vector of voxel
    TwixHeader.NormTra        = twix_obj.hdr.Spice.VoiNormalTra; % Transversal component of normal vector of voxel
end
TwixHeader.TablePosSag    = twix_obj.hdr.Dicom.lGlobalTablePosSag; % Sagittal table position [mm]
TwixHeader.TablePosCor    = twix_obj.hdr.Dicom.lGlobalTablePosCor; % Coronal table position [mm]
TwixHeader.TablePosTra    = twix_obj.hdr.Dicom.lGlobalTablePosTra; % Transversal table position [mm]
% If a parameter is set to zero (e.g. if no voxel rotation is
% performed), the respective field is left empty in the TWIX file. This
% case needs to be intercepted. Setting to the minimum possible value.
VoI_Params = {'VoI_InPlaneRot','VoI_RoFOV','VoI_PeFOV','VoIThickness','NormCor','NormSag','NormTra', ...
              'PosCor','PosSag','PosTra','TablePosSag','TablePosCor','TablePosTra'};
for pp = 1:length(VoI_Params)
    if isempty(TwixHeader.(VoI_Params{pp}))
        TwixHeader.(VoI_Params{pp}) = realmin('double');
    end
end
geometry.size.VoI_RoFOV     = TwixHeader.VoI_RoFOV; % Voxel size in readout direction [mm]
geometry.size.VoI_PeFOV     = TwixHeader.VoI_PeFOV; % Voxel size in phase encoding direction [mm]
geometry.size.VoIThickness  = TwixHeader.VoIThickness; % Voxel size in slice selection direction [mm]
geometry.pos.PosCor         = TwixHeader.PosCor; % Coronal coordinate of voxel [mm]
geometry.pos.PosSag         = TwixHeader.PosSag; % Sagittal coordinate of voxel [mm]
geometry.pos.PosTra         = TwixHeader.PosTra; % Transversal coordinate of voxel [mm]
geometry.pos.TablePosSag    = TwixHeader.TablePosSag; % Sagittal table position [mm]
geometry.pos.TablePosCor    = TwixHeader.TablePosCor; % Coronal table position [mm]
geometry.pos.TablePosTra    = TwixHeader.TablePosTra; % Transversal table position [mm]
geometry.rot.VoI_InPlaneRot = TwixHeader.VoI_InPlaneRot; % Voxel rotation in plane
geometry.rot.NormCor        = TwixHeader.NormCor; % Coronal component of normal vector of voxel
geometry.rot.NormSag        = TwixHeader.NormSag; % Sagittal component of normal vector of voxel
geometry.rot.NormTra        = TwixHeader.NormTra; % Transversal component of normal vector of voxel

%Calculate t and ppm arrays using the calculated parameters:
dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
spectralwidth=1/dwelltime;
f=(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)));
ppm=f/(Bo*42.577);
% Siemens data assumes the center frequency to be 4.7 ppm:
centerFreq = 4.7;
ppm=ppm + centerFreq;

t=0:dwelltime:(sz(1)-1)*dwelltime;

%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=fftshift(fft(fids,[],1));
out.sz=sz;
out.ppm=ppm;
out.t=t;
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date='';
out.dims.t = 1;
out.dims.coils = 2;
out.dims.averages = 3;
out.dims.subSpecs = 0;
out.dims.extras = 0;
out.Bo=Bo;
out.averages=0;
out.rawAverages=0;
out.subspecs=0;
out.rawSubspecs=0;
out.seq=seq;
out.te=0;
out.tr = twix_obj.hdr.MeasYaps.alTR{1};  %Franck Lamberton
out.pointsToLeftshift=leftshift;
out.centerFreq = centerFreq;
out.geometry = geometry;
out.PatientPosition = twix_obj.hdr.Config.PatientPosition;
out.flags = [];
out.Manufacturer = 'Siemens';if isfield(twix_obj.hdr.Dicom,'SoftwareVersions')
    out.software = [version ' ' twix_obj.hdr.Dicom.SoftwareVersions];
else
    out.software = version;
end
end