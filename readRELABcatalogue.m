function [sampleIDs,sampleNames,generalType1s,generalType2s,type1s,type2s,...
    subTypes,minSizes,maxSizes,particulates,textures,origins,...
    spectrumIDs,specCodes,wavelength_strts,wavelength_ends,resolutions,...
    incidents,emissions] = readRELABcatalogue(relabLoc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
% [sampleIDs,numRecords,generalType,categories,wavelength_strt,wavelength_end,fnames]
%        = readRELABcatalogue(workbookFile, sheetName, range)
% Inputs
%   relabloc: path to the relab database
% Outputs
%   sampleIDs: original IDs in Janice library
%   sampleName:
%   generalType1: mineral group
%   generalType2: 
%   categories: class names
%   fnames: cell array of the file names
%   incidents: incident angles;
%   emissions: emission angles;
%   specCodes: show which sensor was used for each measurement
%   Type : showing type of the samples
%   minSize : minimum particle size [um]
%   maxSize : maximum particle size [um]
%   source : source
%   samplename : sample name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If no range is specified, read all data
%if nargin <= 2 || isempty(range)
    range = '';
%end

%% Import the data
%from location of Relab navigate to spectra and sample catalogue
catalogueLoc=strcat(relabLoc,'/catalogues');
sample_catalogue_loc=strcat(catalogueLoc,'/Sample_Catalogue.xlsx');
spectra_catalogue_loc=strcat(catalogueLoc,'/Spectra_Catalogue.xlsx');

%read both catalogues
[~,~,sample_catalogue] = xlsread(sample_catalogue_loc);
[~,~,spectra_catalogue] = xlsread(spectra_catalogue_loc);

%eliminate empty columns
sample_catalogue=sample_catalogue(:,1:20);
spectra_catalogue=spectra_catalogue(:,1:22);
%
% read information about samples
sample_catalogue(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),sample_catalogue)) = {''};
sampleIDs = sample_catalogue(2:end,1);
sampleNames = sample_catalogue(2:end,2);
generalType1s = sample_catalogue(2:end,7);
generalType2s = sample_catalogue(2:end,8);
type1s = sample_catalogue(2:end,9);
type2s = sample_catalogue(2:end,10);
subTypes = sample_catalogue(2:end,11);
minSizes = sample_catalogue(2:end,13);
maxSizes = sample_catalogue(2:end,14);
particulates = sample_catalogue(2:end,15);
textures = sample_catalogue(2:end,16);
origins = sample_catalogue(2:end,17);

% read spectral data 
spectra_catalogue(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),spectra_catalogue)) = {''};
spectrumIDs = spectra_catalogue(2:end,1);
sampleIDs_spectra = spectra_catalogue(2:end,2);
specCodes = spectra_catalogue(2:end,5);
wavelength_strts = spectra_catalogue(2:end,6);
wavelength_ends = spectra_catalogue(2:end,7);
resolutions = spectra_catalogue(2:end,8);
incidents = spectra_catalogue(2:end,9);
emissions = spectra_catalogue(2:end,10);

clear spectra_catalogue
clear sample_catalogue;

% matching the sample catalogue with spectra catalogue
sampleIDs_new = [];
sampleNames_new = [];
generalType1s_new = [];
generalType2s_new = [];
type1s_new = [];
type2s_new = [];
subTypes_new = [];
minSizes_new = [];
maxSizes_new = [];
particulates_new = [];
textures_new = [];
origins_new = [];

spectrumIDs_new = [];
specCodes_new = [];
wavelength_strts_new = [];
wavelength_ends_new = [];
resolutions_new = [];
incidents_new = [];
emissions_new = [];

for i = 1:length(sampleIDs)
    sampleID = sampleIDs{i};
    idx_bin = cellfun(@(x) strcmp(lower(x),lower(sampleID)),sampleIDs_spectra);
    idx = find(idx_bin);
    sampleIDs_new = [sampleIDs_new repmat({sampleID},[1, length(idx)]) ];
    sampleNames_new = [sampleNames_new repmat({sampleNames{i}}, [1, length(idx)])];
    generalType1s_new = [generalType1s_new repmat({generalType1s{i}}, [1, length(idx)])];
    generalType2s_new = [generalType2s_new repmat({generalType2s{i}}, [1, length(idx)])];
    type1s_new = [type1s_new repmat({type1s{i}}, [1, length(idx)])];
    type2s_new = [type2s_new repmat({type2s{i}}, [1, length(idx)])];
    subTypes_new = [subTypes_new repmat({subTypes{i}}, [1, length(idx)])];
    minSizes_new = [minSizes_new repmat({minSizes{i}}, [1, length(idx)])];
    maxSizes_new = [maxSizes_new repmat({maxSizes{i}}, [1, length(idx)])];
    particulates_new = [particulates_new repmat({particulates{i}}, [1, length(idx)]) ];
    textures_new = [textures_new repmat({textures{i}}, [1, length(idx)]) ];
    origins_new = [origins_new repmat({origins{i}}, [1, length(idx)]) ];
    
    spectrumIDs_new = [spectrumIDs_new {spectrumIDs{idx}}];
    specCodes_new = [specCodes_new {specCodes{idx}}];
    wavelength_strts_new = [wavelength_strts_new {wavelength_strts{idx}}];
    wavelength_ends_new = [wavelength_ends_new {wavelength_ends{idx}}];
    resolutions_new = [resolutions_new {resolutions{idx}}];
    incidents_new = [incidents_new {incidents{idx}}];
    emissions_new = [emissions_new {emissions{idx}}];
    
end
sampleIDs = sampleIDs_new;
sampleNames = sampleNames_new;
generalType1s = generalType1s_new;
generalType2s = generalType2s_new;
type1s = type1s_new;
type2s = type2s_new;
subTypes = subTypes_new;
minSizes = minSizes_new;
maxSizes = maxSizes_new;
particulates = particulates_new;
textures = textures_new;
origins = origins_new;

spectrumIDs = spectrumIDs_new;
specCodes = specCodes_new;
wavelength_strts = wavelength_strts_new;
wavelength_ends = wavelength_ends_new;
resolutions = resolutions_new;
incidents = incidents_new;
emissions = emissions_new;

clear *_new

% remove invalid data
wavelength_strts = cellfun(@(x) x(1),wavelength_strts,'ErrorHandler',@errorfun);
wavelength_ends = cellfun(@(x) x(1),wavelength_ends,'ErrorHandler',@errorfun);
 
valid_idx = find(and(~isnan(wavelength_strts),~isnan(wavelength_ends)));
 
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
 
 
end