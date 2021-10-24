function [value, gradient, hessian] = bananaValleyObjective(x)
  %BANANAVALLEYOBJECTIVE Nonlinear function R^2->R
  
  %% Purpose: 
  % test function mapping x -> f=100*(x(2)-x(1)^2)^2+(1-x(1))^2+2;
  
  %% Input Definition:
  % x: column vector in R^2 (domain space)
  
  %% Output Definition:
  % value: real number, evaluation of bananaValleyObjective at x 
  % gradient: real column vector in R^2, evaluation of the gradient with respect to x at x
  % hessian: real 2x2 matrix, evaluation of the hessian with respect to x at x
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient,myHessian]=bananaValleyObjective([1;1]);
  % should return
  % myValue=2; myGradient=[0;0]; myHessian=[802,-400;-400,200];
  
  %% Input verification:
  
  if ~isequal(size(x), [2,1])
    error('Size of x is wrong.');    
  end
  
  %% Implementation:
  
  value = 100*(x(2)-x(1)^2)^2+(1-x(1))^2+2;
  
  if nargout > 1
    
    f_dx1 = -400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
    f_dx2 = 200*(x(2)-x(1)^2);
    
    gradient = [f_dx1;f_dx2];
    
  end
  
  if nargout > 2
    f_dx11 = -400*x(2)+1200*x(1)^2+2;
    f_dx12 = -400*x(1);
    f_dx22 = 200;
    
    hessian = [f_dx11 f_dx12; f_dx12 f_dx22];
  end
end
