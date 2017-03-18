function [spclib_relab] = readRELABdata( dir_path )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% USAGE
%   [spclib_relab] = readRELABdata( dir_path)
% Inputs
%   dir_path : path to the directory which contain relab database
%   sheetName: sheetname
%   wavelength: [w_strt, w_end]
%   prunLvl: the least number of memgbers in one class
% Outputs
%   spclib_relab: struct, having following fields
%       reflectance:
%       wavelength
%       std_measurement
%       fpath: path to the text files
%
%       ***** property found in the spectrum catalogue *****
%       spectrumID
%       specCode: show which sensor was used for each measurement
%       wavelength_strt
%       wavelength_end
%       resolution
%       incident: incident angles;
%       emission: emission angles;
%
%       ***** property found in the sample catalogue *****
%       sampleID: original sampleIDs in janice spectral library
%       sampleName
%       generalType1: mineral group
%       generalType2
%       type1
%       type2
%       subType
%       minSize : minimum particle size [um]
%       maxSize : maximum particle size [um]
%       particulate
%       texture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read specificaitons
if 0
    [sampleIDs,sampleNames,generalType1s,generalType2s,type1s,type2s,...
        subTypes,minSizes,maxSizes,particulates,textures,origins,...
        spectrumIDs,specCodes,wavelength_strts,wavelength_ends,resolutions,...
        incidents,emissions] = readRELABcatalogue(dir_path);
    
    
    save([dir_path 'sampleIDs'],'sampleIDs');
    save([dir_path 'sampleNames'],'sampleNames');
    save([dir_path 'generalType1s'],'generalType1s');
    save([dir_path 'generalType2s'],'generalType2s');
    save([dir_path 'type1s'], 'type1s');
    save([dir_path 'type2s'],'type2s');
    save([dir_path 'subTypes'],'subTypes');
    save([dir_path 'minSizes'],'minSizes');
    save([dir_path 'maxSizes'],'maxSizes');
    save([dir_path 'particulates'],'particulates');
    save([dir_path 'textures'],'textures');
    save([dir_path 'origins'],'origins');
    save([dir_path 'spectrumIDs'],'spectrumIDs');
    save([dir_path 'specCodes'],'specCodes');
    save([dir_path 'wavelength_strts'],'wavelength_strts');
    save([dir_path 'wavelength_ends'],'wavelength_ends');
    save([dir_path 'resolutions'],'resolutions');
    save([dir_path 'incidents'],'incidents');
    save([dir_path 'emissions'],'emissions');
end
% 
load([dir_path 'sampleIDs']);
load([dir_path 'sampleNames']);
load([dir_path 'generalType1s']);
load([dir_path 'generalType2s']);
load([dir_path 'type1s']);
load([dir_path 'type2s']);
load([dir_path 'subTypes']);
load([dir_path 'minSizes']);
load([dir_path 'maxSizes']);
load([dir_path 'particulates']);
load([dir_path 'textures']);
load([dir_path 'origins']);
load([dir_path 'spectrumIDs']);
load([dir_path 'specCodes']);
load([dir_path 'wavelength_strts']);
load([dir_path 'wavelength_ends']);
load([dir_path 'resolutions']);
load([dir_path 'incidents']);
load([dir_path 'emissions']);


% clean up some of the names of the categories
for i=1:length(generalType1s)
    i
    sampleIDs{i} = strtrim(sampleIDs{i});
    if ~isnumeric(sampleNames{i})
        sampleNames{i} = strtrim(sampleNames{i});
    end
    generalType1s{i} = strtrim(generalType1s{i});
    generalType2s{i} = strtrim(generalType2s{i});
    type1s{i} = strtrim(type1s{i});
    type2s{i} = strtrim(type2s{i});
    subTypes{i} = strtrim(subTypes{i});
    particulates{i} = strtrim(particulates{i});
    textures{i} = strtrim(textures{i});
    origins{i} = strtrim(origins{i});
    spectrumIDs{i} = strtrim(spectrumIDs{i});
    specCodes{i} = strtrim(specCodes{i});
end

% eliminate other than mineral
if 0
    idx_bool1 = cellfun(@(x) strcmp('mineral',lower(x)),generalType1s);
    idx_bool2 = cellfun(@(x) strcmp('mixture',lower(x)),type1s);
    idx_bool2 = (idx_bool2==0);
    idx_bool3 = cellfun(@(x) strcmp('yes',lower(x)),particulates);
    idx_bool4 = cellfun(@(x) strcmp('bd-vnir',lower(x)),specCodes);
    idx_bool5 = cellfun(@(x) ~isempty(strfind(lower(x),'particulate')),textures);
    idx_bool6 = cellfun(@(x) ~strcmp('mixture',lower(x)),type1s);
    idx_bool7 = cellfun(@(x) ~strcmp('olivine',lower(x)),subTypes);
    idx_bool8 = cellfun(@(x) ~strcmp('plagioclase',lower(x)),subTypes);
    idx_bool = and(idx_bool1,idx_bool2);
    idx_bool = and(idx_bool,idx_bool3);
    idx_bool = and(idx_bool,idx_bool4);
    idx_bool = and(idx_bool,idx_bool5);
    idx_bool = and(idx_bool,idx_bool6);
    idx_bool = and(idx_bool,idx_bool7);
    idx_bool = and(idx_bool,idx_bool8);

    sampleIDs = sampleIDs(idx_bool);
    sampleNames = sampleNames(idx_bool);
    generalType1s = generalType1s(idx_bool);
    generalType2s = generalType2s(idx_bool);
    type1s = type1s(idx_bool);
    type2s = type2s(idx_bool);
    subTypes = subTypes(idx_bool);
    minSizes = minSizes(idx_bool);
    maxSizes = maxSizes(idx_bool);
    particulates = particulates(idx_bool);
    textures = textures(idx_bool);
    origins = origins(idx_bool);

    spectrumIDs = spectrumIDs(idx_bool);
    specCodes = specCodes(idx_bool);
    wavelength_strts = wavelength_strts(idx_bool);
    wavelength_ends = wavelength_ends(idx_bool);
    resolutions = resolutions(idx_bool);
    incidents = incidents(idx_bool);
    emissions = emissions(idx_bool);
end

% extract only spectra acquired at 'San Carlo'
if 0
    ptr = lower('San Carlo');
    idx_bool = cellfun(@(x) ~isempty(strfind(lower(x),ptr)),origins);
    sampleIDs = sampleIDs(idx_bool);
    sampleNames = sampleNames(idx_bool);
    generalType1s = generalType1s(idx_bool);
    generalType2s = generalType2s(idx_bool);
    type1s = type1s(idx_bool);
    type2s = type2s(idx_bool);
    subTypes = subTypes(idx_bool);
    minSizes = minSizes(idx_bool);
    maxSizes = maxSizes(idx_bool);
    particulates = particulates(idx_bool);
    textures = textures(idx_bool);
    origins = origins(idx_bool);

    spectrumIDs = spectrumIDs(idx_bool);
    specCodes = specCodes(idx_bool);
    wavelength_strts = wavelength_strts(idx_bool);
    wavelength_ends = wavelength_ends(idx_bool);
    resolutions = resolutions(idx_bool);
    incidents = incidents(idx_bool);
    emissions = emissions(idx_bool);
end
    
%
% read the all spectral files
spclib = cell(1,length(spectrumIDs));
fpaths = cell(1,length(spectrumIDs));
err_idx = [];
err_idx2 = [];
warning off;
for i=1:length(spectrumIDs)
     try
        fname = spectrumIDs{i};
        % set path to the reflectance files
        PIcodes = strsplit(sampleIDs{i},'-');
        fpath = [dir_path,'/data/', lower(PIcodes{2}),'/',lower(PIcodes{1}),'/',fname,'.asc'];
        fpath_sub = ['data/', lower(PIcodes{2}),'/',lower(PIcodes{1}),'/',fname,'.asc'];
        % if isempty(strfind(fname,'bkr'))
            %fprintf([fname,'\n']);
            if exist(fpath,'file')
                spc = readRELABfile(fpath);
                spclib{i} = spc;
                fpaths{i} = fpath;
            else
                fprintf('warning: %s does not exist(%d)\n',fpath,i);
                err_idx = [err_idx i];
            end
        %else
        %  fprintf('warning: %s (%d)\n',fname,i);
        %        err_idx = [err_idx i];
        %end
     catch err
          err_idx2 = [err_idx2,i];
          fname = spectrumIDs{i};
          fprintf('%s error (%d)\n',fpath,i);
     end
end

valid_idx = setdiff(1:length(spectrumIDs),[err_idx err_idx2]);
% further prunning the invalid spectra
sampleIDs = sampleIDs(valid_idx);
sampleNames = sampleNames(valid_idx);
generalType1s = generalType1s(valid_idx);
generalType2s = generalType2s(valid_idx);
type1s = type1s(valid_idx);
type2s = type2s(valid_idx);
subTypes = subTypes(valid_idx);
minSizes = minSizes(valid_idx);
maxSizes = maxSizes(valid_idx);
particulates = particulates(valid_idx);
textures = textures(valid_idx);
origins = origins(valid_idx);

spectrumIDs = spectrumIDs(valid_idx);
specCodes = specCodes(valid_idx);
wavelength_strts = wavelength_strts(valid_idx);
wavelength_ends = wavelength_ends(valid_idx);
resolutions = resolutions(valid_idx);
incidents = incidents(valid_idx);
emissions = emissions(valid_idx);

spclib = spclib(valid_idx);
fpaths = fpaths(valid_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prun spectra based on start and end of the wavelength
if 0
    strt_cnd = (wavelength_strts <= wavelength(1));
    end_cnd = (wavelength(2) <= wavelength_ends);
    idx_bool = and(strt_cnd,end_cnd);

    sampleIDs = sampleIDs(idx_bool);
    sampleNames = sampleNames(idx_bool);
    generalType1s = generalType1s(idx_bool);
    generalType2s = generalType2s(idx_bool);
    type1s = type1s(idx_bool);
    type2s = type2s(idx_bool);
    subTypes = subTypes(idx_bool);
    minSizes = minSizes(idx_bool);
    maxSizes = maxSizes(idx_bool);
    particulates = particulates(idx_bool);
    textures = textures(idx_bool);
    origins = origins(idx_bool);

    spectrumIDs = spectrumIDs(idx_bool);
    specCodes = specCodes(idx_bool);
    wavelength_strts = wavelength_strts(idx_bool);
    wavelength_ends = wavelength_ends(idx_bool);
    resolutions = resolutions(idx_bool);
    incidents = incidents(idx_bool);
    emissions = emissions(idx_bool);

    spclib = spclib(idx_bool);
    fpaths = fpaths(idx_bool);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% further prunning spectra based on start and end of the wavelength
% because some of the spectra have limitted range in the specified range.
if 0
    valid_idx = [];
    for i=1:length(spclib)
        spc = spclib{i};
        fprintf([spectrumIDs{i} '\n']);
        if and(min(spc(:,1))<wavelength(1),wavelength(2)<max(spc(:,1)))
            valid_idx = [valid_idx i];
        end
    end

    sampleIDs = sampleIDs(valid_idx);
    sampleNames = sampleNames(valid_idx);
    generalType1s = generalType1s(valid_idx);
    generalType2s = generalType2s(valid_idx);
    type1s = type1s(valid_idx);
    type2s = type2s(valid_idx);
    subTypes = subTypes(valid_idx);
    minSizes = minSizes(valid_idx);
    maxSizes = maxSizes(valid_idx);
    particulates = particulates(valid_idx);
    textures = textures(valid_idx);

    spectrumIDs = spectrumIDs(valid_idx);
    specCodes = specCodes(valid_idxl);
    wavelength_strts = wavelength_strts(valid_idx);
    wavelength_ends = wavelength_ends(valid_idx);
    resolutions = resolutions(valid_idx);
    incidents = incidents(valid_idx);
    emissions = emissions(valid_idx);
    
    spclib = spclib(valid_idx);
    fpaths = fpaths(valid_idx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first, prun spectra which do not have enough bands.
% tmp = spc_data<10^(-10);
% tmp = sum(tmp,1);
% tmp = tmp <2;
% spc_data = spc_data(:,tmp);
% usgsOriNames = usgsOriNames(tmp);
% numRecords = numRecords(tmp);
% generalType = generalType(tmp);
% categories = categories(tmp);
% catNums = catNums(tmp);
% purityLvls = purityLvls(tmp);
% endmembers_boolean = endmembers_boolean(tmp);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prunning with regard to class size
if 0
    diff_clses = unique(subTypes);
    indicies = [];
    for i=1:length(diff_clses)
        cls = diff_clses{i};
        if ~strcmp(cls,'')
            cls_ar = find(cellfun(@(x) strcmp(x,cls),subTypes));
            if size(cls_ar(:),1)>=prunLvl
                indicies = [indicies cls_ar];
            end
        end
    end
    indicies = sort(indicies,'ascend');
    sampleIDs = sampleIDs(indicies);
    sampleNames = sampleNames(indicies);
    generalType1s = generalType1s(indicies);
    generalType2s = generalType2s(indicies);
    type1s = type1s(indicies);
    type2s = type2s(indicies);
    subTypes = subTypes(indicies);
    minSizes = minSizes(indicies);
    maxSizes = maxSizes(indicies);
    particulates = particulates(indicies);
    textures = textures(indicies);
    origins = origins(indicies);

    spectrumIDs = spectrumIDs(indicies);
    specCodes = specCodes(indicies);
    wavelength_strts = wavelength_strts(indicies);
    wavelength_ends = wavelength_ends(indicies);
    resolutions = resolutions(indicies);
    incidents = incidents(indicies);
    emissions = emissions(indicies);

    spclib = spclib(indicies);
    fpaths = fpaths(indicies);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spclib_relab = [];
for i = 1:length(sampleIDs)
    spclib_relab(i).reflectance = spclib{i}(:,2);
    spclib_relab(i).wavelength = spclib{i}(:,1);
    spclib_relab(i).std_measurement = spclib{i}(:,3);
    spclib_relab(i).fpath = fpaths{i};
    
    % spectrum catalogue properties
    spclib_relab(i).spectrumID = spectrumIDs{i};
    spclib_relab(i).specCode = specCodes{i};
    spclib_relab(i).wavelength_strt = wavelength_strts(i);
    spclib_relab(i).wavelength_end = wavelength_ends(i);
    spclib_relab(i).resolution = resolutions{i};
    spclib_relab(i).incident = incidents{i};
    spclib_relab(i).emission = emissions{i}; 
    
    
    % sample catalogue properties
    spclib_relab(i).sampleID = sampleIDs{i};
    spclib_relab(i).sampleName = sampleNames{i};
    spclib_relab(i).generalType1 = generalType1s{i};
    spclib_relab(i).generalType2 = generalType2s{i};
    spclib_relab(i).type1 = type1s{i};
    spclib_relab(i).type2 = type2s{i};
    spclib_relab(i).subType = subTypes{i};
    spclib_relab(i).minSize = minSizes{i};
    spclib_relab(i).maxSize = maxSizes{i};
    spclib_relab(i).particulate = particulates{i};
    spclib_relab(i).texture = textures{i};
    
end
end