function [value, gradient, hessian] = simpleValleyObjective(x,p)
  %SIMPLEVALLEYOBJECTIVE Nonlinear function R^2->R
  
  %% Purpose: 
  % test function depending on parameter vector p and mapping x -> cosh(x(1)) + p(1)*(x(2)-1)^2+ p(2)
  
  %% Input Definition:
  % x: column vector in R^2 (domain space)
  % p: column vector in R^2 (parameter space), p(1) is scaling in x(2) direction, p(2) is elevation
  
  %% Output Definition:
  % value: real number, evaluation of simpleValleyObjective at x for parameters p 
  % gradient: real column vector in R^2, evaluation of the gradient with respect to x at x for parameters p
  % hessian: real 2x2 matrix, evaluation of the hessian with respect to x at x for parameters p
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient,myHessian]=simpleValleyObjective([0;1],[1;2]);
  % should return
  % myValue=3; myGradient=[0;0]; myHessian=[1,0;0,2];
  
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
    
    gradient = [f_dx1;f_dx2];
    
  end
  
  if nargout > 2
    f_dx11 = cosh(x(1));
    f_dx12 = 0;
    f_dx22 = 2*p(1);
    
    hessian = [f_dx11 f_dx12; f_dx12 f_dx22];
  end
end
