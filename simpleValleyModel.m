function [value, gradient_x, gradient_p] = simpleValleyModel(x,p)
  %SIMPLEVALLEYMODEL Nonlinear function R^2->R
  
  %% Purpose: 
  % test function depending on parameter vector p and mapping x -> cosh(x(1)) + p(1)*(x(2)-1)^2+ p(2)
  
  %% Input Definition:
  % x: column vector in R^2 (domain space)
  % p: column vector in R^2 (parameter space), p(1) is scaling in x(2) direction, p(2) is elevation
  
  %% Output Definition:
  % value: real number, evaluation of simpleValleyModel at x for parameters p 
  % gradient_x: real column vector in R^2, evaluation of the gradient with respect to x at x for parameters p
  % gradient_p: real column vector in R^2, evaluation of the gradient with respect to p at p for domain points x
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient_x,myGradient_p]=simpleValleyModel([0;1],[1;2]);
  % should return
  % myValue=3; myGradient_x=[0;0]; myGradient_p=[0;1];
  
  %% Input verification:
  
  if ~isequal(size(x), [2,1])
    error('Size of x is wrong.');    
  end
  if ~isequal(size(p), [2,1])
    error('Size of p is wrong.');    
  end
  
  %% Implementation:
  
  value = cosh(x(1)) + p(1)*(x(2)-1)^2 + p(2);
  
  if nargout > 1
    
    f_dx1 = sinh(x(1));
    f_dx2 = 2*p(1)*(x(2)-1);
    
    gradient_x = [f_dx1;f_dx2];
    
  end
  
  if nargout > 2
      
     gradient_p = [(x(2)-1)^2;1];
      
  end
end
