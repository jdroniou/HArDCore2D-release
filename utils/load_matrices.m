function [mat]=load_matrices(N)
%
% Load computed matrices meshX_SystemMatrix.mtx
%

for i=1:N;
  matrix_file=strcat('mesh',num2str(i),'_SystemMatrix.mtx');
  mat{i} = mmread(matrix_file);
end;

