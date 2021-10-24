function [xmin] = projectedNewtonDescent(f, P, x0, eps, verbose)
%PROJECTEDNEWTONDESCENT Find minimal point of function using global Newton descent
%direction subject to projection

%% Purpose:
% Find xmin to satisfy norm(xmin - P(xmin - gradf(xmin)))<=eps
% Iteration: x_k = P(x_k + t_k * d_k)
% d_k is the Newton direction of the reduced hessian. If a descent direction check fails,
% choose d_k to be the steepest descent.
% t_k results from projected backtracking


%% Input Definition:
% f: function handle of type [value, gradient, hessian] = f(x).
% P: box projection handle of type [projectedX, activeIndexSet] = P(x).
% x0: column vector in R^n (domain point), starting point.
% eps: positive value, tolerance for termination. Default value: 1.0e-3.
% verbose: bool, if set to true, verbose information is displayed

%% Output Definition:
% xmin: column vector in R^n (domain point, box constraints), satisfies norm(xmin - P(xmin - gradf(xmin)))<=eps

%% Required files:
% [x] = PrecCGSolver(A,b,delta)
% [t] = projectedBacktrackingSearch(f, P, x, d, sigma)

%% Test cases:
% [xmin]=projectedNewtonDescent(@(x)simpleValleyObjective(x,[1;1]),@(x)projectionInBox(x,[1;1],[2;2],1.0e-6), [2;2], 1.0e-3, true);
% should return
% xmin close to [1;1] (exact xmin depends on choice of sigma in projectedbacktracking);
%
% [xmin]=projectedNewtonDescent(@(x)nonlinearObjective(x), @(x)projectionInBox(x,[1;1],[2;2],1.0e-6),[0.1;0.1], 1.0e-3, true);
% should return
% xmin close to [1;1]
%
% [xmin]=projectedNewtonDescent(@(x)nonlinearObjective(x), @(x)projectionInBox(x,[-2;-2],[2;2],1.0e-6),[1.5;2.0], 1.0e-3, true);
% should return
% xmin close to [-0.26;0.21] (if it is [0.2608;-0.2086] then maybe your
% reduction is done wrongly)
%
% [xmin]=projectedNewtonDescent(@(x)bananaValleyObjective(x), @(x)projectionInBox(x,[-10;-10],[10;10],1.0e-6),[0;1], 1.0e-6, true);
% should return
% xmin close to [1;1] in less than 25 iterations. If you have too much
% iterations, then maybe the hessian is used wrongly.

%% Input verification:

try
    [x,A]=P(x0);
    [value, gradient, hessian] = f(x);
catch
    error('evaluation of function handle or projection handle failed!');
end

if (eps <= 0)
    error('range of eps is wrong!');
end

if nargin < 5
    verbose = false;
end

%% Implementation:
% Hints:
% 1. Whenever x changes, you need to update dependend values properly!
% 2. Use [x,A] = P(x0) to get the current projected point and its active
% set
% 3. Hint for matrix reduction: You can access and overwrite specific
% rows/columns of a matrix by special calls. Execute the following example
% to get the idea:
%
% L = [1 2 3; 4 5 6; 7 8 9]
% R = [10 10 10; 10 10 10; 10 10 10 ]
% index = [1 3]
% L(:,index)=R(:,index)
% R(index,:)=L(index,:)
%
% 4. eye generates a unit matrix.
%
% 5. Keep track of the iterations with
% if verbose
%   countIter=countIter+1;
% end

  

if verbose
    disp('Start projectedBFGSDescent...');
    countIter = 0;
end

%static
delta=1.0e-6;
sigma=1.0e-4;
n=length(x0);

x_k=P(x0);

while norm(x_k - P(x_k - gradient) ) > eps
  [value, gradient, hessian] = f(x_k);
  B_k =  hessian;
  d_k = PrecCGSolver(B_k, -gradient, delta);
  
  if gradient' * d_k >=0
    d_k = -gradient;
  endif
  
  t_k = projectedBacktrackingSearch(f, P, x_k, d_k, sigma);
  x_k = P(x_k + t_k * d_k);
  
 
 endwhile
xmin = x_k

if verbose
    [value, gradient] = f(xmin);
    disp(sprintf('projectedNewtonDescent terminated after %i steps with norm of stationarity =%d\n',countIter, norm(x - P(x - gradient))));
end

end