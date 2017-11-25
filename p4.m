% This is not a function file
1;

% Generates test data for numerical experiment
function [A, b, xtrue, M, c] = gen_data(N, k)
  [Q1, R] = qr(randn(N));
  [Q2, R] = qr(randn(N));
  D = diag(10 .^(k/N:k/N:k));
  A = Q1*D*Q2;
  b = ones(N,1);
  z = zeros(N, N);
  xtrue = A\b;
  M = [z A; A' z];
  c = [b ; b];
end

% gmres routine
function x = gmr(A, b, max_iterations, threshold)
  n = length(A);
  m = max_iterations;
  x = zeros(n, 1);
  %use x as the initial vector
  r=b-A*x;

  b_norm = norm(b);
  error = norm(r)/b_norm;

  %initialize the 1D vectors
  sn = zeros(m,1);
  cs = zeros(m,1);
  e1 = zeros(n,1);
  e1(1) = 1;
  r_norm=norm(r);
  Q(:,1) = r/r_norm;
  beta = r_norm*e1;
  for k = 1:m                                   
    
    %run arnoldi
    [H(1:k+1,k) Q(:,k+1)] = arnoldi(A, Q, k);
    
    %eliminate the last element in H ith row and update the rotation matrix
    [H(1:k+1,k) cs(k) sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k);
    
    %update the residual vector
    beta(k+1) = -sn(k)*beta(k);
    beta(k)   = cs(k)*beta(k);
    error  = abs(beta(k+1)) / b_norm;
    
    if ( error <= threshold)
      break;
    end
  end

  %calculate the result
  y = H(1:k,1:k) \ beta(1:k);
  x = x + Q(:,1:k)*y; 
end

% Arnoldi Function                  

function [h q] = arnoldi(A, Q, k)
  q = A*Q(:,k);
  for i = 1:k,
    h(i)= q'*Q(:,i);
    q = q - h(i)*Q(:,i);
  end
  h(k+1) = norm(q);
  q = q / h(k+1);
end

%  Applying Givens Rotation to H col                  
function [h, cs_k sn_k] = apply_givens_rotation(h, cs, sn, k)
  %apply for ith column
  for i = 1:k-1,                              
    temp     =  cs(i)*h(i) + sn(i)*h(i+1);
    h(i+1) = -sn(i)*h(i) + cs(i)*h(i+1);
    h(i)   = temp;
  end
  
  %update the next sin cos values for rotation
  [cs_k sn_k] = givens_rotation( h(k), h(k+1));
  
  %eliminate H(i+1,i)
  h(k) = cs_k*h(k) + sn_k*h(k+1);
  h(k+1) = 0.0;
end

% Calculate the Given rotation matrix
function [cs sn] = givens_rotation(v1, v2)
  if (v1==0)
    cs = 0;
    sn = 1;
  else
    t=sqrt(v1^2+v2^2);
    cs = abs(v1) / t;
    sn = cs * v2 / v1;
  end
end

% Conjugate gradient routine
function [x, re, rr] = cg(A, b, x1, maxitr, tol)
  n = length(b);
  x = zeros(n, 1);
  nx1 = norm(x1);
  r0 = b - A*x;
  p = r0; r = r0;
  itr = 0;
  while true
    a = (r0'*r0)/(p'*A*p);
    x =  x + a*p;
    r =  r - a*A*p;
    re = norm(x-x1)/nx1;
    itr = itr + 1;
    if re < tol || itr > maxitr
      break
    end
    beta =  (r'*r)/(r0'*r0);
    p = r + beta*p;
    r0 = r;
  end
  rr = norm(b-a*x)/norm(b); 
end

% Computes errors for bicg and gmres
function e = err(A, b, x1)
  nx1 = norm(x1); nb = norm(b);
  [x, ~] = bicgstab(A, b, 1e-6, 20);
  e = [norm(x-x1)/nx1, norm(b-A*x)/nb];
  x = gmr(A, b, 20, 1e-6);
  e = [e [norm(x-x1)/nx1, norm(b-A*x)/nb]]; 
end

% Computes average error for bicg, gmres, cg
function [e, ecg] = avg_err(N, k, n)
  e = zeros(4, 1); ecg = [0 0];
  for j = 1:n
      [A, b, x1, M, c] = gen_data(N, k);
      e = e + err(A, b, x1);
      [~, re, rr] = cg(M, c, M\c, 20, 1e-6);
      ecg = ecg + [re, rr];
  end
  e = e/n; ecg = ecg/n;
end

function plt(k)
  bre = zeros(10, 1); % Relative error for bicg 
  gre = zeros(10, 1); % Relative error for gmres
  brr = zeros(10, 1); % relative residue for bicg
  grr = zeros(10, 1); % relative residue for gmres
  cre = zeros(10, 1); % Relative error for cg
  crr = zeros(10, 1); % relative residue for cg
  for i = 1:10
    [e, ecg] = avg_err(100*i, k, 1);
    bre(i) = e(1);
    brr(i) = e(2);
    gre(i) = e(3);
    grr(i) = e(4);
    cre(i) = ecg(1);
    crr(i) = ecg(2);
   end
   X = [100:100:1000];
   fig = figure('visible', 'off');
   plot(X, bre);
   xlabel('size of matrix');
   ylabel('relative error');
   saveas(fig, strcat('img4/Rel_err_bicg_', num2str(k)), 'png');
   plot(X, gre, '--');
   xlabel('size of matrix');
   ylabel('relative error');
   saveas(fig, strcat('img4/Rel_err_gmres_', num2str(k)), 'png');
   plot(X, brr);
   xlabel('size of matrix');
   ylabel('relative residue');
   saveas(fig, strcat('img4/Rel_res_bicg_', num2str(k)), 'png');
   plot(X, grr);
   xlabel('size of matrix');
   ylabel('relative residue');
   saveas(fig, strcat('img4/Rel_res_gmres_', num2str(k)), 'png');
   plot(X, cre, '--');
   xlabel('size of matrix');
   ylabel('relative error');
   saveas(fig, strcat('img4/Rel_err_cg_', num2str(k)), 'png');
   plot(X, crr);
   xlabel('size of matrix');
   ylabel('relative residue');
   saveas(fig, strcat('img4/Rel_res_cg_', num2str(k)), 'png');
end

for k = 3:6
  plt(k);
end

disp('Plots have been generated')