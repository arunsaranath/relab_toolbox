%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
if ismac
    path2relab = '/Volumes/SED/data/relab/';
else
end

dir_path = [path2relab 'RelabDB2016Dec'];

sheetName = 'Sheet1';
prunLvl = 1;
wavelength = [1000 2600];
%[usgsOriNames,numRecords,broadCategories,categories,catNums,purityLvls,spc_data,wavelengths,endmembers_boolean]...
%                = usgs_all_read( fname,sheet,prunLvl);

[spclib_relab] = readRELABdata( dir_path );



% extract spectra which are only taken by BD-VNIR
ptr = 'BD-VNIR';
uspecCodes = unique(specCodes);
% remove invalid data base on angles
incidents_new = cellfun(@(x) isnumeric(x),incidents,'ErrorHandler',@errorfun);
emissions_new = cellfun(@(x) isnumeric(x),emissions,'ErrorHandler',@errorfun);


%%
valid_idx = find(and(incidents_new,emissions_new));
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
specCodes = specCodes(valid_idx);
wavelength_strts = wavelength_strts(valid_idx);
wavelength_ends = wavelength_ends(valid_idx);
resolutions = resolutions(valid_idx);
incidents = incidents(valid_idx);
emissions = emissions(valid_idx);
uspecCodes = unique(specCodes);

spclib = spclib(valid_idx);

incidents = cell2mat(incidents);
emissions = cell2mat(emissions);

%%
minSizes = cell2mat(minSizes);
maxSizes = cell2mat(maxSizes);
resolutions = cell2mat(resolutions);


% remove the minSize=0 and maxSize=0
valid_idx = maxSizes>eps;
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
specCodes = specCodes(valid_idx);
wavelength_strts = wavelength_strts(valid_idx);
wavelength_ends = wavelength_ends(valid_idx);
resolutions = resolutions(valid_idx);
incidents = incidents(valid_idx);
emissions = emissions(valid_idx);
uspecCodes = unique(specCodes);

spclib = spclib(valid_idx);



%%
mat_name = sprintf('Relab_mnrl_prtclt_bdvnir_prun%d.mat',prunLvl);

save(mat_name,'sampleIDs','sampleNames','generalType1s','generalType2s',...
    'type1s','type2s',...
    'subTypes','minSizes','maxSizes','particulates','textures',...
    'spectrumIDs','specCodes','wavelength_strts','wavelength_ends','resolutions',...
    'incidents','emissions','spclib');

%%

diff_cats = unique();
correctlabels = zeros(1,length(subTypes));
nClass = size(diff_cats);
diff_clses = 1:size(diff_cats);

%nBand = size(wavelengths(:),1);
N = size(correctlabels,2);
numspectra = N;
%bands=size(min_data,1);
%endmembers_idx = find(endmembers_boolean);

% reproduce correctlabels for setting continuing number
%correctlabels_new = zeros(size(correctlabels));
%for i=1:nClass
%    correctlabels_new(correctlabels==diff_clses(i)) = i;
%end
%correctlabels = correctlabels_new;

save(mat_name, 'spclib', 'correctlabels',...
     'nClass', 'N', 'sampleIDs','fnames', 'generalType',...
    'numRecords','categories','diff_cats','incidents','emissions',...
    'specCodes','Type','source','samplename','minSize','maxSize');


%%
prunLvl = 1;
mat_name = sprintf('Relab_mnrl_prtclt_bdvnir_prun%d.mat',prunLvl);
load(mat_name,'sampleIDs','sampleNames','generalType1s','generalType2s',...
    'type1s','type2s',...
    'subTypes','minSizes','maxSizes','particulates','textures',...
    'spectrumIDs','specCodes','wavelength_strts','wavelength_ends','resolutions',...
    'incidents','emissions','spclib');
%min_data = spc_data;
diff_clses = unique(subTypes);
diff_cats = unique(subTypes);
correctlabels = zeros(1,length(subTypes));
for i=1:size(diff_cats,2)
    cls = diff_clses{i};
    cls_ar = find(cellfun(@(x) strcmp(x,diff_cats{i}),subTypes));
    if ~strcmp(cls,'')
        cls_ar = find(cellfun(@(x) strcmp(x,cls),subTypes));
        if length(cls_ar)>=1
            % fprintf('%s: %d\n',diff_cats{i},length(cls_ar));
        
            correctlabels(cls_ar)=i;
            fprintf('%s: %d\n',diff_cats{i},length(cls_ar));
            fig = figure(1);
            set_figsize(fig,900,400);
            axh = subplot(1,1,1);
            hold(axh,'on');
            cols = distinguishable_colors(length(cls_ar));
            for j=1:size(cls_ar,2)
                idx = cls_ar(j);
                spc = spclib{idx};
                plot(axh,spc(:,1),spc(:,2),'Color',cols(j,:));
            end
            legend(spectrumIDs(cls_ar),'Location','EastOutside');
            title(diff_cats{i});
            axis tight;
            xlim([350,2600]);
            ylim([0,1]);
            xlabel('wavelength [nm]');
            ylabel('Reflectance');
            try
                saveas(fig,['./plots_cleaned/' diff_cats{i},'.fig']);
                export_fig(fig,['./plots/' diff_cats{i},'.png'],'-transparent');
            catch err
                saveas(fig,['./plots/' num2str(i),'_.fig']);
                export_fig(fig,['./plots/' num2str(i),'_.png'],'-transparent');
            end

            close(fig);
        end
    end
end

%%
usubTypes = unique(subTypes);
corType1 = [];
corType2 = [];
indicies = [];
l = [];
for i=1:length(diff_clses)
    cls = usubTypes{i};
    if ~strcmp(cls,'')
        cls_ar = find(cellfun(@(x) strcmp(x,cls),subTypes));
        if length(cls_ar)>=1
            fprintf('%s: %d\n',diff_cats{i},length(cls_ar));
            l = [l,length(cls_ar)];
            indicies = [indicies cls_ar];
            corType1 = [corType1 {type1s{cls_ar(1)}}];
            corType2 = [corType2 {type2s{cls_ar(1)}}];
        end
    end
    
end




figure();
hist(l,104);

ures = unique(res);
indices = [];
for i=1:length(ures)
    cls = ures(i);
    if ~strcmp(cls,'')
        cls_ar = find(res==cls);
        if length(cls_ar)>=1
            fprintf('%d: %d\n',ures(i),length(cls_ar));
        end
    end
end

%%
% construct the library
% class actinolite
categories = [];

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

spectrumIDs_new = [];
specCodes_new = [];
wavelength_strts_new = [];
wavelength_ends_new = [];
resolutions_new = [];
incidents_new = [];
emissions_new = [];

spclib_new = [];

for i=1:length(spclib)
    flag = 0;
    s = subTypes{i};
    % class actinolite
    if strcmp(s,'Actinolite') || strcmp(s,'Calcic Amphibole Actinolite')
        flag = 1;
        categories = [categories {'Actinolite'}];
    % class allophane
    elseif strcmp(s,'Allophane')
        flag = 1;
        categories = [categories {'Allophane'}];
    elseif strcmp(s,'Alunite')
        flag = 1;
        categories = [categories {'Alunite'}];
    elseif strcmp(s,'Amarantite')
        flag = 1;
        categories = [categories {'Amarantite'}];
    elseif strcmp(s,'Zeolite Analcime')
        flag = 1;
        categories = [categories {'Analcime'}];
    elseif strcmp(s,'Anglesite')
        flag = 1;
        categories = [categories {'Anglesite'}];
    elseif strcmp(s,'Anhydrite')
        flag = 1;
        categories = [categories {'Anhydrite'}];
    elseif strcmp(s,'Ankerite')
        flag = 1;
        categories = [categories {'Ankerite'}];
    elseif strcmp(s,'Annite (Biotite)')
        flag = 1;
        categories = [categories {'Annite'}];
    % class 4 Anorthite
    elseif strcmp(s,'Anorthite')
        flag = 1;
        categories = [categories {'Anorthite'}];
    % class 5 Aragonite
    elseif strcmp(s,'Aragonite')
        flag = 1;
        categories = [categories {'Aragonite'}];
    % class 6 Artinite
    elseif strcmp(s,'Artinite')
        flag = 1;
        categories = [categories {'Artinite'}];
    elseif strcmp(s,'Atacamite')
        flag = 1;
        categories = [categories {'Atacamite'}];
    elseif strcmp(s,'Attapulgite')
        flag = 1;
        categories = [categories {'Attapulgite'}];
    elseif strcmp(s,'Augite')
        flag = 1;
        categories = [categories {'Augite'}];
    elseif strcmp(s,'Azurite')
        flag = 1;
        categories = [categories {'Azurite'}];
    elseif strcmp(s,'Barite')
        flag = 1;
        categories = [categories {'Barite'}];
    elseif strcmp(s,'Beidellite')
        flag = 1;
        categories = [categories {'Beidellite'}];
    elseif strcmp(s,'Berthierine')
        flag = 1;
        categories = [categories {'Berthierine'}];
    elseif strcmp(s,'Beryl')
        flag = 1;
        categories = [categories {'Beryl'}];
    elseif strcmp(s,'Biotite') || strcmp(s,'Lepidomelane')
        flag = 1;
        categories = [categories {'Biotite'}];
    elseif strcmp(s,'Bornite')
        flag = 1;
        categories = [categories {'Bornite'}];
    elseif strcmp(s,'Botryogen')
        flag = 1;
        categories = [categories {'Botryogen'}];
    elseif strcmp(s,'Boussingaultite (NH4)2Mg(SO4)2-6(H2O)')
        flag = 1;
        categories = [categories {'Boussingaultite'}];
    % class 7 Brucite
    elseif strcmp(s,'Brucite')
        flag = 1;
        categories = [categories {'Brucite'}];
    elseif strcmp(s,'Brugnatellite')
        flag = 1;
        categories = [categories {'Brugnatellite'}];
    elseif strcmp(s,'Buddingtonite')
        flag = 1;
        categories = [categories {'Buddingtonite'}];
    % class 8 Calcite
    elseif strcmp(s,'Calcite')
        flag = 1;
        categories = [categories {'Calcite'}];
    elseif strcmp(s,'Celadonite')
        flag = 1;
        categories = [categories {'Celadonite'}];
    elseif strcmp(s,'Cerussite')
        flag = 1;
        categories = [categories {'Cerussite'}];
    elseif strcmp(s,'Chalcopyrite')
        flag = 1;
        categories = [categories {'Chalcopyrite'}];
    elseif strcmp(s,'Chamosite')
        flag = 1;
        categories = [categories {'Chamosite'}];
    % class 9 Chromite
    elseif strcmp(s,'Chromite')
        flag = 1;
        categories = [categories {'Chromite'}];
    elseif strcmp(s,'Alexandrite')
        flag = 1;
        categories = [categories {'Chrysoberyl'}];
    elseif strcmp(s,'Chrysocolla')
        flag = 1;
        categories = [categories {'Chrysocolla'}];
    elseif strcmp(s,'Clinochlore')
        flag = 1;
        categories = [categories {'Chlinochlore'}];
    elseif strcmp(s,'Coalingite')
        flag = 1;
        categories = [categories {'Coalingite'}];
    % class 10 Copiapite
    elseif strcmp(s,'Copiapite')
        flag = 1;
        categories = [categories {'Copiapite'}];
    elseif strcmp(s,'Coquimbite')
        flag = 1;
        categories = [categories {'Coquimbite'}];
    elseif strcmp(s,'Corundum')
        flag = 1;
        categories = [categories {'Corundum'}];
    elseif strcmp(s,'Cristobalite')
        flag = 1;
        categories = [categories {'Cristobalite'}];
    elseif strcmp(s,'Cronstedtite')
        flag = 1;
        categories = [categories {'Cronstedtite'}];
    elseif strcmp(s,'Dawsonite')
        flag = 1;
        categories = [categories {'Dawsonite'}];
    % class 11 Diopside
    elseif strcmp(s,'Diopside')
        flag = 1;
        categories = [categories {'Diopside'}];
    % class 12 Dolomite
    elseif strcmp(s,'Dolomite')
        flag = 1;
        categories = [categories {'Dolomite'}];
    % class 13 Enstatite
    elseif strcmp(s,'Enstatite') || strcmp(s,'Enstatite Orthopyroxene Pyroxene')...
        || strcmp(s,'Pyroxene Enstatite') || strcmp(s,'Pyroxene Orthopyroxene Enstatite')
        flag = 1;
        categories = [categories {'Enstatite'}];
    % class 14 Fayalite
    elseif strcmp(s,'Olivine Fayalite')
        flag = 1;
        categories = [categories {'Fayalite'}];
    elseif strcmp(s,'Ferrihydrite')
        flag = 1;
        categories = [categories {'Ferrihydrite'}];
    elseif strcmp(s,'Fibroferrite')
        flag = 1;
        categories = [categories {'Fibroferrite'}];
    elseif strcmp(s,'Franklinite')
        flag = 1;
        categories = [categories {'Franklinite'}];
    elseif strcmp(s,'Gahnite')
        flag = 1;
        categories = [categories {'Gahnite'}];
    elseif strcmp(s,'Gaspeite')
        flag = 1;
        categories = [categories {'Gaspeite'}];
    elseif strcmp(s,'Gaylussite')
        flag = 1;
        categories = [categories {'Gaylussite'}];
    elseif strcmp(s,'Geikielite')
        flag = 1;
        categories = [categories {'Geikielite'}];
    elseif strcmp(s,'Gibbsite')
        flag = 1;
        categories = [categories {'Gibbsite'}];
    elseif strcmp(s,'Zeolite Gismondine')
        flag = 1;
        categories = [categories {'Gismondine'}];
    elseif strcmp(s,'Glauconite')
        flag = 1;
        categories = [categories {'Glauconite'}];
    elseif strcmp(s,'Goethite')
        flag = 1;
        categories = [categories {'Goethite'}];
    elseif strcmp(s,'Graftonite')
        flag = 1;
        categories = [categories {'Graftonite'}];
    elseif strcmp(s,'Graphite')
        flag = 1;
        categories = [categories {'Graphite'}];
    elseif strcmp(s,'Greenalite')
        flag = 1;
        categories = [categories {'Greenalite'}];
    % class 15 Gypsum
    elseif strcmp(s,'Gypsum')
        flag = 1;
        categories = [categories {'Gypsum'}];
    elseif strcmp(s,'Gyrolite')
        flag = 1;
        categories = [categories {'Gyrolite'}];
    elseif strcmp(s,'Halite')
        flag = 1;
        categories = [categories {'Halite'}];
    elseif strcmp(s,'Halloysite')
        flag = 1;
        categories = [categories {'Halloysite'}];
    elseif strcmp(s,'Halotrichite')
        flag = 1;
        categories = [categories {'Halotrichite'}];
    elseif strcmp(s,'Hausmannite')
        flag = 1;
        categories = [categories {'Hausmannite'}];
    elseif strcmp(s,'Hexahydrite')
        flag = 1;
        categories = [categories {'Hexahydrite'}];
    elseif strcmp(s,'Hibonite')
        flag = 1;
        categories = [categories {'Hibonite'}];
    elseif strcmp(s,'Hissingerite')
        flag = 1;
        categories = [categories {'Hisingerite'}];
    % class 16 Hematite
    elseif strcmp(s,'Hematite')
        flag = 1;
        categories = [categories {'Hematite'}];
    % class 17 Hydromagnesite
    elseif strcmp(s,'Hydromagnesite')
        flag = 1;
        categories = [categories {'Hydromagnesite'}];
    % class 18 Hypersthene
    elseif strcmp(s,'Hypersthene') || strcmp(s,'Pyroxene Orthopyroxene Hypersthene')
        flag = 1;
        categories = [categories {'Hypersthene'}];
    % class 19 Jarosite
    elseif strcmp(s,'Jarosite')
        flag = 1;
        categories = [categories {'Jarosite'}];
    % class 20 Kaolinite
    elseif strcmp(s,'Clay Kaolinite') || strcmp(s,'Kaolinte')
        flag = 1;
        categories = [categories {'Kaolinite'}];
    % class 21 Magnesite
    elseif strcmp(s,'Magnesite')
        flag = 1;
        categories = [categories {'Magnesite'}];
    % class 22 Montmorillonite
    elseif strcmp(s,'Montmorillonite') || strcmp(s,'Beidellite Ca-saturated')...
        || strcmp(s,'Bentonite') || strcmp(s,'Montmorillonite Beidellite Ca-saturated')...
        || strcmp(s,'Smectite Na-Montmorillonite') || strcmp(s,'Smectite, Ca-Montmorillonite')
        flag = 1;
        categories = [categories {'Montmorillonite'}];
    % class 23 Nontronite
    elseif strcmp(s,'Nontronite')
        flag = 1;
        categories = [categories {'Nontronite'}]
    % class 24 Quartz
    elseif strcmp(s,'Quartz')
        flag = 1;
        categories = [categories {'Quartz'}];
    % class 25 Ripidolite
    elseif strcmp(s,'Ripidolite')
        flag = 1;
        categories = [categories {'Ripidolite'}];
    % class 26 Scapolite
    elseif strcmp(s,'Scapolite') || strcmp(s,'Scaporite')
        flag = 1;
        categories = [categories {'Scapolite'}];
    % class 27 Sepiolite
    elseif strcmp(s,'Sepiolite')
        flag = 1;
        categories = [categories {'Sepiolite'}];
    % class 28 Smithsonite
    elseif strcmp(s,'Smithsonite')
        flag = 1;
        categories = [categories {'Smithsonite'}];
    % class 29 Sphene
    elseif strcmp(s,'Sphene Titanite')
        flag = 1;
        categories = [categories {'Sphene'}];
    elseif strcmp(s,'Calcic Amphibole Tremolite')
        flag = 1;
        categories = [categories {'Tremolite'}];
    end
    
    if flag
        sampleIDs_new = [sampleIDs_new {sampleIDs{i}} ];
        sampleNames_new = [sampleNames_new {sampleNames{i}}];
        generalType1s_new = [generalType1s_new {generalType1s{i}}];
        generalType2s_new = [generalType2s_new {generalType2s{i}}];
        type1s_new = [type1s_new {type1s{i}}];
        type2s_new = [type2s_new {type2s{i}}];
        subTypes_new = [subTypes_new {subTypes{i}}];
        minSizes_new = [minSizes_new minSizes(i)];
        maxSizes_new = [maxSizes_new maxSizes(i)];
        particulates_new = [particulates_new {particulates{i}} ];
        textures_new = [textures_new {textures{i}}];

        spectrumIDs_new = [spectrumIDs_new {spectrumIDs{i}}];
        specCodes_new = [specCodes_new {specCodes{i}}];
        wavelength_strts_new = [wavelength_strts_new wavelength_strts(i)];
        wavelength_ends_new = [wavelength_ends_new wavelength_ends(i)];
        resolutions_new = [resolutions_new resolutions(i)];
        incidents_new = [incidents_new incidents(i)];
        emissions_new = [emissions_new emissions(i)];
        spclib_new = [spclib_new spclib(i)];
    end
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

spectrumIDs = spectrumIDs_new;
specCodes = specCodes_new;
wavelength_strts = wavelength_strts_new;
wavelength_ends = wavelength_ends_new;
resolutions = resolutions_new;
incidents = incidents_new;
emissions = emissions_new;

spclib = spclib_new;

prunLvl = 5;
mat_name = sprintf('Relab_mnrl_prtclt_bdvnir_prun%d_cleaned.mat',prunLvl);
save(mat_name,'sampleIDs','sampleNames','generalType1s','generalType2s',...
    'type1s','type2s',...
    'subTypes','minSizes','maxSizes','particulates','textures',...
    'spectrumIDs','specCodes','wavelength_strts','wavelength_ends','resolutions',...
    'incidents','emissions','spclib','categories');

%%
prunLvl = 5;
mat_name = sprintf('Relab_mnrl_prtclt_bdvnir_prun%d_cleaned.mat',prunLvl);
load(mat_name,'sampleIDs','sampleNames','generalType1s','generalType2s',...
    'type1s','type2s',...
    'subTypes','minSizes','maxSizes','particulates','textures',...
    'spectrumIDs','specCodes','wavelength_strts','wavelength_ends','resolutions',...
    'incidents','emissions','spclib');
%min_data = spc_data;
diff_clses = unique(categories);
for i=1:size(diff_clses,2)
    cls = diff_clses{i};
    cls_ar = find(cellfun(@(x) strcmp(x,diff_clses{i}),categories));
    if ~strcmp(cls,'')
        cls_ar = find(cellfun(@(x) strcmp(x,cls),categories));
        if length(cls_ar)>=1
            % fprintf('%s: %d\n',diff_cats{i},length(cls_ar));
        
            correctlabels(cls_ar)=i;
            fprintf('%s: %d\n',diff_clses{i},length(cls_ar));
            fig = figure(1);
            set_figsize(fig,900,400);
            axh = subplot(1,1,1);
            hold(axh,'on');
            cols = distinguishable_colors(length(cls_ar));
            for j=1:size(cls_ar,2)
                idx = cls_ar(j);
                spc = spclib{idx};
                plot(axh,spc(:,1),spc(:,2),'Color',cols(j,:));
            end
            legend(spectrumIDs(cls_ar),'Location','EastOutside');
            title(diff_clses{i});
            axis tight;
            xlim([350,2600]);
            ylim([0,1]);
            xlabel('wavelength [nm]');
            ylabel('Reflectance');
            try
                saveas(fig,['./plots_cleaned/' diff_clses{i},'.fig']);
                export_fig(fig,['./plots_cleaned/' diff_clses{i},'.png'],'-transparent');
            catch err
                saveas(fig,['./plots_cleaned/' num2str(i),'_.fig']);
                export_fig(fig,['./plots_cleaned/' num2str(i),'_.png'],'-transparent');
            end

            close(fig);
        end
    end
end



