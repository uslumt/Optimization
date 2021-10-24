function [value, gradient_x, gradient_p] = exponentialModel(x,p)
  %EXPONENTIALMODEL Nonlinear function R->R
  
  %% Purpose: 
  % test function depending on parameter vector p and mapping x -> p(1)+p(2)*exp(p(3)*x)
  
  %% Input Definition:
  % x: scalar in R (domain space)
  % p: column vector in R^3 (parameter space), p(1) is constant value, p(2)
  % is scaling, p(3) is growth rate
  
  %% Output Definition:
  % value: real number, evaluation of exponentialModel at x for parameters p 
  % gradient_x: real column vector in R^2, evaluation of the gradient with respect to x at x for parameters p
  % gradient_p: real column vector in R^2, evaluation of the gradient with respect to p at p for domain points x
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient_x,myGradient_p]=exponentialModel(0,[1;-1;1]);
  % should return
  % myValue=0; myGradient_x=-1; myGradient_p=[1;1;0];
  
  %% Input verification:
  
  if ~isequal(size(x), [1,1])
    error('Size of x is wrong.');    
  end
  if ~isequal(size(p), [3,1])
    error('Size of p is wrong.');    
  end
  
  %% Implementation:
  
  value = p(1)+p(2)*exp(p(3)*x);
  
  if nargout > 1    
   
    gradient_x = p(2)*exp(p(3)*x)*p(3);
    
  end
  
  if nargout > 2
      
     gradient_p = [1;exp(p(3)*x);p(2)*exp(p(3)*x)*x];
      
  end
end
