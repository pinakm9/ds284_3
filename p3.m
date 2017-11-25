% This is not a function file
1;

% Generates test data for numerical experiment
function [A, b, xtrue] = gen_data(k)
  N = 400;
  [Q, R] = qr(randn(N));
  D = diag(10.^(k/N:k/N:k));
  A = Q*D*Q';
  b = ones(N,1);
  xtrue = A\b;
end

% Steepest descent routine
function [x, re, rr] = sd(A, b, x1)
  n = length(b);
  nb = norm(b);
  x0 = rand(n, 1);
  nx1 = norm(x1);
  re = [norm(x0-x1)/nx1];
  rr = [norm(b-A*x0)/nb];
  while true
    d = b - A*x0;
    x = x0 + ((d'*d)/(d'*A*d))*d;
    e = norm(x-x1)/nx1;
    re = [re e];
    rr = [rr norm(b-A*x)/nb];
    if e < 1e-4
      break
    end
    x0 = x;
  end
end

% Conjugate gradient routine
function [x, re, rr] = cg(A, b, x1)
  n = length(b);
  nb = norm(b);
  x = zeros(n, 1);
  nx1 = norm(x1);
  r0 = b - A*x;
  p = r0; r = r0;
  re = [norm(x-x1)/nx1];
  rr = [norm(b-A*x)/nb];
  while true
    a = (r0'*r0)/(p'*A*p);
    x =  x + a*p;
    r =  r - a*A*p;
    e = norm(x-x1)/nx1;
    re = [re e];
    rr = [rr norm(b-A*x)/nb];
    if e < 1e-4
      break
    end
    beta =  (r'*r)/(r0'*r0);
    p = r + beta*p;
    r0 = r;
  end
end

% Plotter for errors, returns iterations required
function itr = plt(func, k, name)
  [A, b, x1] = gen_data(k);
  [~, re, rr] = func(A, b, x1);
  fig = figure('visible', 'off');
  itr = length(re);
  X = 1:1:itr;
  plot(X, re);
  hold on
  plot(X, rr);
  legend('relative error', 'relative residue');
  xlabel('iteration');
  saveas(fig, strcat(name, num2str(k)), 'png');
end

n = 4;
itrc = zeros(n, 1);
itrs = zeros(n, 1);
for k=1:n
  itrc(k) = plt(@cg, k, 'img3/con_gra_'); 
  itrs(k) = plt(@sd, k, 'img3/ste_des_'); 
end 
fig = figure('visible', 'off');
X = 1:1:n;
plot(X, itrc);
xlabel('k'), ylabel('iterations required');
saveas(fig, 'img3/con_gra_itr', 'png'); 
plot(X, itrs);
xlabel('k'), ylabel('iterations required');
saveas(fig, 'img3/ste_des_itr', 'png'); 
disp('Plots have been generated');                                      