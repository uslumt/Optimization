function [value, gradient, hessian] = benchmarkObjective(x,p)
  %BENCHMARKOBJECTIVE Nonlinear function R^3->R
  
  %% Purpose: 
  % test function depending on parameters p and mapping x ->p(1)*(x(2)+1)*x(1)^2+exp(p(2)*x(3)+1)*x(2)^2+p(3)*sqrt(x(1)+1)*x(3)^2
  
  %% Input Definition:
  % x: column vector in R^3 (domain space), x(1) > -1 must hold 
  % p: column vector in R^3 (parameter space)
  
  %% Output Definition:
  % value: real number, evaluation of benchmarkObjective at x for parameters p 
  % gradient: real column vector in R^3, evaluation of the gradient with respect to x of benchmarkObjective at x for parameters p
  % hessian: real 3x3 matrix, evaluation of the hessian with respect to x of benchmarkObjective at x for parameters p
    
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient,myHessian]=benchmarkObjective([0;0;-1/2],[3;2;16]);
  % should return
  % myValue=4; myGradient=[2;0;-16]; myHessian=[5,0,-8;0,2,0;-8,0,32];
  
  %% Input verification:
  if ~isequal(size(x), [3,1])
    error('Size of x is wrong.');    
  end
  
  if (x(1)<=-1)
    error('x(1)<=-1, f is not continuously differentiable.');    
  end
  
  if ~isequal(size(p), [3,1])
    error('Size of p is wrong.');    
  end
  
  %% Implementation:
  
  value = p(1)*(x(2)+1)*x(1)^2+exp(p(2)*x(3)+1)*x(2)^2+p(3)*sqrt(x(1)+1)*x(3)^2;
  
  if nargout > 1
    f_dx1 = 2*p(1)*(x(2)+1)*x(1)+p(3)/2*x(3)^2/sqrt(x(1)+1);
    f_dx2 = p(1)*x(1)^2+2*exp(p(2)*x(3)+1)*x(2);
    f_dx3 = p(2)*exp(p(2)*x(3)+1)*x(2)^2+2*p(3)*sqrt(x(1)+1)*x(3);
    
    gradient = [f_dx1;f_dx2;f_dx3];
  end
  
  if nargout > 2
    f_dx11 = 2*p(1)*(x(2)+1)-p(3)/4*x(3)^2/sqrt((x(1)+1)^3);
    f_dx12 = 2*p(1)*x(1);
    f_dx13 = p(3)*x(3)/sqrt(x(1)+1);
    
    f_dx22 = 2*exp(p(2)*x(3)+1);
    f_dx23 = 2*p(2)*exp(p(2)*x(3)+1)*x(2);
    
    f_dx33 = p(2)^2*exp(p(2)*x(3)+1)*x(2)^2+2*p(3)*sqrt(x(1)+1);
    
    hessian = [f_dx11 f_dx12 f_dx13; f_dx12 f_dx22 f_dx23; f_dx13 f_dx23 f_dx33];
  end
  
end