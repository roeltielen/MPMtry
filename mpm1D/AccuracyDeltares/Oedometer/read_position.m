function [S1] = readvtk(this)
fileID = fopen(this);
D = textscan(fileID, '%s');
S=D{1};
S1 = str2double(S);
%S1(end)
end
