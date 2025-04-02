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
io_writeniimrs_mod(Spec_ste, OutLoc, {}, OspVer, 1);
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
io_writeniimrs_mod(Spec_Herc, OutLoc, {'DIM_EDIT'}, OspVer, 1);

%% Try and find ref data
InFile = strrep(InFile,'_act.','_ref.');

if exist(InFile,"file")
    Spec = io_loadspec_sdat(InFile,1);

    Spec_ste = op_takeaverages(Spec,1:32);
    Spec_ste.te = 35; 
    Spec_ste.seq = 'PRESS';
    OutLoc = fullfile(OutFolder,['sub-',SubID],['ses-',SesID],'mrs',['sub-',SubID,'_ses-',SesID,'_voi-',VOI,'_acq-press_te-35_mrsref.nii.gz']);
    io_writeniimrs_mod(Spec_ste, OutLoc, {}, OspVer, 0);

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
    io_writeniimrs_mod(Spec_Herc, OutLoc, {'DIM_EDIT'}, OspVer, 0);
end



end


function nii = io_writeniimrs_mod(in, outfile, UserNames, OspreyVersion, WS)
%function nii = io_writeniimrs(in, outfile, UserNames)

% Bring FID array into NIfTI-MRS shape
dims = in.dims;

% Check InNames supplied for additional dimensions
if ~exist('UserNames','var')
    if dims.subSpecs || dims.extras
        error('Need InNames for additional dimensions')
    else
        UserNames = {};
    end
end

%% Manage essential dimensions
if ~isfield(dims, 'x') || dims.x==0
    dims.x = 0;
    dim_1  = 1;
else
    dim_1 = in.sz(dims.x);
end
if ~isfield(dims, 'y') || dims.y==0
    dims.y = 0;
    dim_2  = 1;
else
    dim_2 = in.sz(dims.y);
end
if ~isfield(dims, 'z')|| dims.z==0
    dims.z = 0;
    dim_3  = 1;
else
    dim_3 = in.sz(dims.z);
end

dim_4 = in.sz(dims.t);

ArraySize = [dim_1,dim_2,dim_3,dim_4];
dimorder = [dims.x,dims.y,dims.z,in.dims.t];
dimname = {};

%% Manage non-essential dimensions, starting with defaults
if dims.coils
    ArraySize = [ArraySize, in.sz(dims.coils)];
    dimorder = [dimorder,dims.coils];
    dimname = [dimname,{'DIM_COIL'}];
end
if dims.averages
    ArraySize = [ArraySize, in.sz(in.dims.averages)];
    dimorder = [dimorder,dims.averages];
    dimname = [dimname,{'DIM_DYN'}];
end

%% Now manage user defined dimensions
Cnt = 1;
if in.dims.subSpecs
    ArraySize = [ArraySize, in.sz(dims.subSpecs)];
    dimorder = [dimorder,dims.subSpecs];
    dimname = [dimname,UserNames(Cnt)];
    Cnt = Cnt+1;
end
if in.dims.extras
    ArraySize = [ArraySize, in.sz(dims.extras)];
    dimorder = [dimorder,dims.extras];
    dimname = [dimname,UserNames(Cnt)];
    Cnt = Cnt+1;
end

%% Do some error handling

if length(dimorder)>7
    error('%i dimensions detected, but maximum allowed is 7',length(dimorder))
elseif ~(length(UserNames)==(Cnt-1))
    error('%i user-defined dimensions, and %i dimension labels',Cnt-1,length(UserNames))
elseif length(UserNames)>3
    error('%i user-defined names supplied. Maximum of 3 allowed',length(UserNames))
end

%% Format fid output and export to nii

% If dimensions are ordered correctly, then just add singletons. Otherwise,
% reshape the array to correct for this.
%
% Note: Matlab will remove trailing singletons
if issorted(dimorder)
    fids = reshape(in.fids,ArraySize);
else
    %Permute fid non-singleton dimensions according to default BIDs precedence
    NonSingletonOrdered = permute(in.fids,nonzeros(dimorder));
    
    %Define overall array size including singleton dimensions
    ArraySize_S = double(ArraySize==1);
    ArraySize_S(~(ArraySize==1)) = size(NonSingletonOrdered);
    
    % Reshape output to include correct singletons 
    fids = reshape(NonSingletonOrdered,ArraySize_S);
end

% Initialize the NIfTI struct using the dicm2nii toolbox
% (https://github.com/xiangruili/dicm2nii)
try
    nii = nii_tool('init', fids);
catch ME
    switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            error(['Cannot find the function ''nii_tool.m''.' ...
                ' Please ensure that you have downloaded the required', ...
                ' dcm2nii toolbox (https://github.com/xiangruili/dicm2nii)', ...
                ' and added it to your MATLAB path.']);
        otherwise
            rethrow(ME);
    end
end

%% Update NIfTI-MRS header
newDim = nii.hdr.dim;
newVoxOffset = nii.hdr.vox_offset;
% Add nii_mrs field if needed
if isfield(in,'nii_mrs')
    nii.hdr = in.nii_mrs.hdr;
else
    in.nii_mrs.hdr_ext = osp_generate_nii_hdr_ext(in,OspreyVersion);
    in.nii_mrs.hdr_ext.WaterSuppressed = WS;
    in.nii_mrs.hdr     = osp_generate_nii_hdr(in, nii.hdr,outfile);    
    nii.hdr = in.nii_mrs.hdr;
end
nii.img = conj(fids);

for JJ = 1:length(dimname)
    in.nii_mrs.hdr_ext.(sprintf('dim_%i',JJ+4)) = dimname{JJ};
end

if ~iscell(in.nii_mrs.hdr_ext.SpectrometerFrequency)
    in.nii_mrs.hdr_ext.SpectrometerFrequency = {in.nii_mrs.hdr_ext.SpectrometerFrequency};
end

nii.ext.ecode = 44;
nii.ext.edata_decoded = jsonencode(in.nii_mrs.hdr_ext);
len = int32(numel(nii.ext.edata_decoded));
myEdata = [nii.ext.edata_decoded'];
nii.ext.edata = myEdata;

% nii.hdr.pixdim     = in.nii_mrs.hdr.pixdim;
% nii.hdr.vox_offset     = in.nii_mrs.hdr.vox_offset;

% Save
nii_tool('save', nii, outfile);
end