function [maxeig, mineig, cond, mat] = compute_eigs(N);

% Compute max and min eigenvalues of all matrices

for i=1:N;
  matrix_file = strcat('mesh',num2str(i),'_SystemMatrix.mtx');
  mat{i} = mmread(matrix_file);

  maxeig(i) = eigs(mat{i},1);
  mineig(i) = eigs(mat{i},1,'smallestabs');
  cond(i) = maxeig(i)/mineig(i);
end;


