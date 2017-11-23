% This is not a function file
1;

% Generates test data for numerical experiment
function [A, b, xtrue] = gen_data(N, k)
  [Q, R] = qr(randn(N));
  D = diag(10.^(1:(k-1)/(N-1):k));
  A = Q*D*Q';
  b = ones(N,1);
  xtrue = A\b;
end

% Biconjugate routine
