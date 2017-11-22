% This is not a function file
1;

% Generates test data for numerical experiment
function [A, b, xtrue] = gen_data()
  N = 400;
  k = 6;
  [Q, R] = qr(randn(N));
  D = diag(10.^(1:k/N:k));
  A = Q*D*Q';
  b = ones(N,1);
  xtrue = A/b;
end

% Steepest descent routine
function x = sd(A, b)
  n = length(b);
  x0 = rand(n, 1);
  x1 = A\b;
  nx1 = norm(x1);
  while true
    d = b - A*x0;
    x = x0 + ((d'*d)/(d'*A*d))*d;
    if norm(x-x1)/x1 < 1e-4
      break
    end
    x0 = x
  end
end

sd(rand(4), ones(4, 1))
