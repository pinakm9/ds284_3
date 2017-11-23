% This is not a function file
1;

% Generates test data for numerical experiment
function [A, b, xtrue] = gen_data(k)
  N = 400;
  [Q, R] = qr(randn(N));
  D = diag(10.^(1:(k-1)/(N-1):k));
  A = Q*D*Q';
  b = ones(N,1);
  xtrue = A\b;
end

% Steepest descent routine
function x = sd(A, b, x1)
  n = length(b);
  x0 = rand(n, 1);
  nx1 = norm(x1);
  while true
    d = b - A*x0;
    x = x0 + ((d'*d)/(d'*A*d))*d;
    if norm(x-x1)/nx1 < 1e-4
      break
    end
    x0 = x;
  end
end

% Conjugate gradient routine
function [x, re, rr] = cg(A, b, x1)
  n = length(b);
  x = zeros(n, 1);
  nx1 = norm(x1);
  r0 = b - A*x;
  p = r0; r = r0;
  while true
    a = (r0'*r0)/(p'*A*p);
    x =  x + a*p;
    r =  r - a*A*p;
    if norm(x-x1)/nx1 < 1e-4
      break
    end
    b =  (r'*r)/(r0'*r0);
    p = r + b*p;
    r0 = r;
  end
  re = 9;
  rr = 9; 
end

% 


[A, b, x1] = gen_data(5);
A*cg(A, b, x1)
