function [t] = projectedBacktrackingSearch(f, P, x, d, sigma, verbose)
  %PROJECTEDBACKTRACKINGSEARCH Find stepsize t get sufficient decrease for multidimensional objective along line
  
  %% Purpose:
  % Find t to satisfy f(x+t*d)< f(x) - sigma/t*norm(x-P(x + t*d))^2
  
  %% Input Definition:
  % f: function handle of type [value, gradient] = f(x).
  % P: box projection handle of type [projectedX] = P(x).
  % x: column vector in R^n (domain point) 
  % d: column vector in R^n (search direction)
  % sigma: value in (0,1), marks quality of decrease. Default value:
  % 1.0e-4.
  % verbose: bool, if set to true, verbose information is displayed
  
  %% Output Definition:
  % t: t is set to the biggest 2^m, such that 2^m satisfies the projected sufficient decrease condition
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [t]=projectedBacktrackingSearch(@(x)simpleValleyObjective(x,[1;1]), @(x)projectionInBox(x,[-2;1],[2;2],1.0e-6),[1;1], [-1.99;0], 0.5, true);
  % should return
  % t=0.0625;
  
  %% Input verification:
  
  try
      xk=P(x);
    [value, gradient] = f(xk);
  catch
    error('evaluation of function handle or projection handle failed!'); 
  end
  
  if (gradient'*d>= 0)
    error('descent direction check failed!');    
  end 
  if (sigma <= 0 || sigma >= 1)
    error('range of sigma is wrong!');    
  end  
  if nargin < 6
    verbose = false;
  end
  
  %% Implementation:
  if verbose
    disp('Start projectedBacktrackingSearch...');
  end
  beta = 0.5;
 
  t = 1;
  phi=@(t)(f(P(xk+t*d)));
  f_x = phi(0);
     
  while (phi(t) > f_x - sigma/t * norm(xk - P(xk + t* d))^2 ) 
   t = t * beta;
  end 
   
  if verbose
    disp(sprintf('projectedBacktrackingSearch terminated with t=%d',t));
  end
  
end

