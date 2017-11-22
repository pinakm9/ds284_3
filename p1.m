% This is not a function file
1; 

% Random matrix generator
function A = mat(N)
  A = rand(N) + (2*N)*diag(randperm(N)); 
end

% Gauss Seidel routine
function t = gs(A, b)
  [n, ~] = size(A);
  x = rand(n, 1);
  tic;
  x1 = A\b;
  t0 = toc;
  tic;
  nx1 = norm(x1);
  while true 
    for i = 1:n
      x(i) = (b(i) - A(i, 1:n)*x + A(i,i)*x(i))/A(i,i);
    end
    if norm(x1-x)/nx1 < 1e-3
      break
    end
  end
  t = toc/t0;
end

% Jacobi routine
function t = ja(A, b)
  [n, ~] = size(A);
  x = ones(n, 1);
  tic;
  x1 = A\b;
  t0 = toc;
  tic;
  nx1 = norm(x1);
  x0 = ones(n, 1);
  while true 
    for i = 1:n
      x(i) = (b(i) - A(i, 1:n)*x0 + A(i,i)*x0(i))/A(i,i);
    end 
    if norm(x1-x)/nx1 < 1e-3
      break
    end
   x0 = x; 
  end
  t = toc/t0;
end

% Returns average normalized time taken by Gauss Seidel/Jacobi routine for matrices of size N
function t = avg(func, N, n)
  T = ones(n, 1);
  b = ones(N, 1);
  for i = 1:n
    T(i) = func(mat(N), b);
  end
  t = mean(T);
 end
 
% Plots normalized time taken by Gauss Seidel/Jacobi routine
function plt(func, img_name)
  N = [10, 50, 100, 500];
  t = ones(4, 1);
  for i = 1:4
    t(i) = avg(func, N(i), 100);
  end
  fig = figure("visible", "off");
  plot(N, t);
  xlabel('size of matrix'); ylabel('normalized time'); 
  saveas(fig, img_name, "png");
  disp("Plot has been generated")
end

plt(@gs, "Gauss_Seidel");
plt(@ja, "Jacobi");