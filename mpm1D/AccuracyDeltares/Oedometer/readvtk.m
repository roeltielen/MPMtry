function [S] = readvtk(this)
fileID = fopen(this);
D = textscan(fileID, '%f %f %f');
S=D{2};
%S1 = str2double(S);
%S1(end)
end