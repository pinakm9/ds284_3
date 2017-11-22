% This is not a function file
1;

% Random matrix generator
function A = mat(N)
  [Q, R] = qr(rand(N));
  B = diag(rand(1, N-1),1);
  A = Q*[diag(rand(1, N))+B+B']*Q';
end

% Power mathod for finding largest eigenvalue
function l = pow_l(A)
  [n, ~] = size(A);
  v = zeros(n,1);
  v(1) = 1;
  l0 = v'*A*v;
  while true
    w = A*v;
    v = w/norm(w);
    l = v'*A*v;
    if abs((l-l0)/l0) < 1e-3
      break
    end
    l0 = l;
  end
end

% Power mathod for finding smallest eigenvalue
function l = pow_s(A)
  [n, ~] = size(A);
  v = zeros(n,1);
  v(1) = 1;
  l0 = v'*A*v;
  while true
    w = A\v;
    v = w/norm(w);
    l = v'*A*v;
    if abs((l-l0)/l0) < 1e-3
      break
    end
    l0 = l;
  end
  eig(A)
end

% Rayleigh method for finding largest eigenvalue
function l = ray_l(A)
  [n, ~] = size(A);
  I = eye(n);
  v = zeros(n,1);
  v(1) = 1;
  l0 = 0;
  while true
    w = (A-l0*I)\v;
    v = w/norm(w);
    l = v'*A*v;
    if abs((l-l0)/l0) < 1e-3
      break
    end
    l0 = l;
  end
  eig(A)
end


ray_l(mat(5))