function [xmin] = globalNewtonDescent(f, x0, eps, verbose)
%GLOBALNEWTONDESCENT Find minimal point of function using either Newton or steepest descent direction

%% Purpose:
% Find xmin to satisfy norm(gradf(xmin))<=eps
% Iteration: x_k = x_k + t_k * d_k
% d_k is the Newton direction if it is a descent direction check,
% otherwise choose d_k to be the steepest descent.
% t_k results from Wolfe-Powell


%% Input Definition:
% f: function handle of type [value, gradient, hessian] = f(x).
% x0: column vector in R^n (domain point), starting point.
% eps: positive value, tolerance for termination. Default value: 1.0e-3.
% verbose: bool, if set to true, verbose information is displayed

%% Output Definition:
% xmin: column vector in R^n (domain point), satisfies norm(gradf(xmin))<=eps

%% Required files:
% [x] = PrecCGSolver(A,b,delta)
% [t] = WolfePowellSearch(f, x, d, sigma, rho)

%% Test cases:
% [xmin]=globalNewtonDescent(@(x)nonlinearObjective(x), [-0.01;0.01], 1.0e-6, true);
% should return
% xmin close to [0.26;-0.21] (exact xmin depends on choice of delta in PrecCGSolver);

% [xmin]=globalNewtonDescent(@(x)nonlinearObjective(x), [-0.6;0.6], 1.0e-3, true);
% should return
% xmin close to [-0.26;0.21] (exact xmin depends on choice of delta in PrecCGSolver);

% [xmin]=globalNewtonDescent(@(x)nonlinearObjective(x), [0.6;-0.6], 1.0e-3, true);
% should return
% xmin close to [-0.26;0.21] (exact xmin depends on choice of delta in PrecCGSolver);

%% Input verification:

try
    [value, gradient, hessian] = f(x0);
catch
    error('evaluation of function handle failed!');
end

if (eps <= 0)
    error('range of eps is wrong!');
end

if nargin < 4
    verbose = false;
end

%% Implementation:
% Hints:
% 1. Whenever x changes, you need to update the objective value and
% gradient and hessian properly!
if verbose
    disp('Start globalNewtonDescent...');
    countIter = 0;
end

%Complete the code

x_k = x0;

while norm(gradient) > eps
  [value, gradient, hessian] = f(x_k);
  B = hessian;
  d = PrecCGSolver(B, -gradient, 1.0e-6);
  
  if gradient'*d >=0
    d = -gradient;
  endif
  
  t = WolfePowellSearch(f, x_k, d, 1.0e-3, 1.0e-2);
  
  x_k = x_k + t*d;
  [value, gradient, hessian] = f(x_k);
endwhile
xmin =x_k

if verbose
    [value, gradient] = f(xmin);
    disp(sprintf('globalNewtonDescent terminated after %i steps with norm of gradient =%d\n',countIter, norm(gradient)));
end

end