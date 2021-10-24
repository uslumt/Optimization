function [t] = WolfePowellSearch(f, x, d, sigma, rho, verbose)
  %WOLFEPOWELLSEARCH Find stepsize t get sufficient decrease and steepness for multidimensional objective along line
  
  %% Purpose:
  % Find t to satisfy f(x+t*d)<=f(x) + t*sigma*gradf(x)'*d and
  % gradf(x+t*d)'*d >= rho*gradf(x)'*d
  
  %% Input Definition:
  % f: function handle of type [value, gradient] = f(x).
  % x: column vector in R^n (domain point) 
  % d: column vector in R^n (search direction)
  % sigma: value in (0,1/2), marks quality of decrease. Default value: 1.0e-3
  % rho: value in (sigma,1), marks quality of steepness. Default value:1.0e-2
  % verbose: bool, if set to true, verbose information is displayed
  
  %% Output Definition:
  % t: t is set, such that t satisfies both Wolfe-Powell conditions
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [t] = WolfePowellSearch(@(x)simpleValleyObjective(x,[0;1]), [-1.01;1], [1;1], 1.0e-3, 1.0e-2, true);
  % should return
  % t=1;
  %
  % [t] = WolfePowellSearch(@(x)simpleValleyObjective(x,[0;1]), [-1.2;1], [0.1;1], 1.0e-3, 1.0e-2, true);
  % should return
  % t=16;
  %
  % [t] = WolfePowellSearch(@(x)simpleValleyObjective(x,[0;1]), [-0.2;1], [1;1], 1.0e-3, 1.0e-2, true);
  % should return
  % t=0.25;
  %
  % [t] = WolfePowellSearch(@(x)nonlinearObjective(x), [0.53;-0.29], [-3.88;1.43], 1.0e-3, 1.0e-2, true);
  % should return
  % t=0.0938;
  
  %% Input verification:
  
  try
    [value, gradient] = f(x);
  catch
    error('evaluation of function handle failed!'); 
  end
  
  if (gradient'*d>= 0)
    error('descent direction check failed!');    
  end 
  if (sigma <= 0 || sigma >= 0.5)
    error('range of sigma is wrong!');    
  end  
  if (rho <= sigma || rho >= 1)
    error('range of rho is wrong!');    
  end   
  if nargin < 6
    verbose = false;
  end
  
  %% Implementation:
  % Hints: 
  % 1. Whenever t changes, you need to update the objective value and
  % gradient properly!
  % 2. Use the return keyword (see documentation)
  if verbose
    disp('Start WolfePowellSearch...');
  end
      
%Complete the code
 
  function w1 = W1(t)
    w1 = f(x + t*d) <= f(x) + t*sigma*gradient'*d
  end
  
  function w2=W2(t)
    [value1, gradient1] = f(x+t*d)
    w2 = gradient1'*d  >= rho*gradient'*d
  end
  
  t = 1
  
  if W1(t) == false
    t = t/2
    while W1(t) == false
      t = t/2
    endwhile
    t_min = t
    t_max = 2*t
     
  elseif W2(t) == true
    t_opt = t
  return ;
  
  else 
  t = 2*t
  
  while W1(t) == true
    t = 2*t
   endwhile
   
   t_min = t/2
   t_max = t
   
 endif
 t = t_min
 
   while W2(t) == false
     t = (t_min + t_max)/2
     
     if W1(t) == true
       t_min = t
     else
       t_max = t
     endif
   endwhile
   
  t  = t_min
   
  if verbose
    disp(sprintf('WolfePowellSearch terminated with t=%d',t));
  end
end

