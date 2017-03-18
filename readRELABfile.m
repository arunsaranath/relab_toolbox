function [ spc ] = readRELABfile( fpath )
%[spc] = readRELABfile(fpath)
%   fpath : filepath to the janice individual files
%   spc : L*3 matrix (1st col : wavelength, 2nd col : reflectance, 3rd col standard deviation)
%try
fid = fopen(fpath);
tmp = fgets(fid);
nRow = str2num(tmp);

spc = nan([nRow,3]);
data = textscan(fid,' %f %f %f ',nRow);
for i=1:size(data,2)
    spc(:,i) = data{i};
end
fclose(fid);
end

