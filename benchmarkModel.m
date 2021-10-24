function [value, gradient_x, gradient_p] = benchmarkModel(x,p)
  %BENCHMARKMODEL Nonlinear function R^3->R
  
  %% Purpose: 
  % test function depending on parameters p and mapping x ->p(1)*(x(2)+1)*x(1)^2+exp(p(2)*x(3)+1)*x(2)^2+p(3)*sqrt(x(1)+1)*x(3)^2
  
  %% Input Definition:
  % x: column vector in R^3 (domain space), x(1) > -1 must hold 
  % p: column vector in R^3 (parameter space)
  
  %% Output Definition:
  % value: real number, evaluation of benchmarkModel at x for parameters p 
  % gradient_x: real column vector in R^3, evaluation of the gradient with respect to x at x for parameters p
  % gradient_p: real column vector in R^3, evaluation of the gradient with respect to p at p for domain points x
    
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient_x,myGradient_p]=benchmarkModel([0;0;-1/2],[3;2;16]);
  % should return
  % myValue=4; myGradient_x=[2;0;-16]; myGradient_p=[0;0;0.25];
  
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
    
    gradient_x = [f_dx1;f_dx2;f_dx3];
  end
  
  if nargout > 2
      R_dp1=(x(2)+1)*x(1)^2;
      R_dp2=x(3)*exp(p(2)*x(3)+1)*x(2)^2;
      R_dp3=sqrt(x(1)+1)*x(3)^2;
      gradient_p=[R_dp1; R_dp2; R_dp3];
  end
  
end