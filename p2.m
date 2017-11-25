% This is not a function file
1;

% Random matrix generator
function A = mat(N)
  [Q, R] = qr(rand(N));
  B = diag(rand(1, N-1),1);
  A = Q*(diag(rand(1, N))+B+B')*Q';
end

% Power mathod 
function [l, t] = pow(A)
  tic;
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
  t = toc;
end

% Rayleigh method
function [l, t] = ray(A)
  tic;
  [n, ~] = size(A);
  v = ones(n, 1);
  l0 = 0;
  I = eye(n);
  while true
    w = (A-l0*I)\v;
    v = w/norm(w);
    l = v'*A*v;
    if l0 ~= 0 && abs((l-l0)/l0) < 1e-3
      break
    end
    l0 = l;
  end
  t = toc;
end

% Returns average time taken by Power/Rayleigh routine for matrices of size N
function t = avg(func, N, n)
  T = ones(n, 1);
  b = ones(N, 1);
  for i = 1:n
    [~, T(i)] = func(mat(N));
  end
  t = mean(T);
 end
 
% Plots time taken by Power/Rayleigh routine
function plt(func, img_name)
  N = [10, 50, 100, 500];
  t = ones(4, 1);
  for i = 1:4
    t(i) = avg(func, N(i), 100);
  end
  fig = figure("visible", "off");
  plot(N, t);
  ylabel("time (s)"); xlabel("size of matrix"); 
  saveas(fig, img_name, "png");
  disp("Plot has been generated")
end

plt(@pow, "img2/Power");
plt(@ray, "img2/Rayleigh");