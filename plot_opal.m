
if ismac
    path2relab = '/Volumes/SED/data/relab/';
else
end
dir_path = [path2relab 'RelabDB2016Dec'];
[spclib_relab] = relab_all_read( dir_path );

idList = {'bir1jb815','bir1jb816','bir1jb817','bir1jb818','bir1jb819'};

kaol1 = searchby('spectrumID',idList,spclib_relab,'COMP_FUNC','strcmpi');

sampIDList = {'JB-JLB-874','JB-JLB-A33','JM-TGS-066','NV-RBS-022','OP-MCG-001',...
    'OP-MCG-002','OP-MCG-003','OP-MCG-004','OP-MCG-005','OP-MCG-006',...
    'OP-MCG-007','OP-MCG-008','OP-MCG-009','OP-MCG-010','OP-MCG-011',...
    'OP-MCG-012','SI-DWM-008-A','SI-DWM-008-B'};

opal2 = searchby('sampleID',sampIDList,spclib_relab,'COMP_FUNC','strcmpi');

figure; hold on;
for i=1:length(opal2)
    wv = opal2(i).wavelength;
    spc = opal2(i).reflectance;
    plot(wv,spc);
end
xlim([1000 2600]);

figure; hold on;
for i=1:length(kaol1)
    wv = kaol1(i).wavelength;
    spc = kaol1(i).reflectance;
    plot(wv,spc);
end
xlim([1000 2600]);
xlabel('Wavelength [nm]');
set(gca,'XTick',[1000 1400 1900 2200 2500]);
ylabel('Reflectance');
grid on;
legend(idList,'Location','EastOutside');


% pdir = '../RelabDB2014Dec/data/jlb/jb/';

% pdir = '../RelabDB2014Dec/data/mcg/op/';
% idList = cellstr(num2str((1:12)', 'bkr1op%03d'));
% idList2 = cellstr(num2str((1:12)', 'c1op%02d'));
% 
% cols2 =    [0             0    1.0000;
%             0        0.5000         0;
%             1.0000        0         0;
%             0        0.7500    0.7500;
%             0.7500        0    0.7500;
%             0.7500   0.7500         0;
%             0.2500   0.2500    0.2500];
% 
% cols = distinguishable_colors(length(idList2));
%         
% figure;
% hold on;
% for i=1:length(idList2)
%     [ spc ] = janice_spc_read( [pdir idList2{i} '.asc'] );
%     plot(spc(:,1),spc(:,2)-0.1*i,'-','Color',cols(i,:));
% end
% 
% xlim([1000 2600]);

export_fig(gcf,'opal_mcgop_c1op.png','-r150','-transparent');



