function [value, gradient, hessian] = nonlinearObjective(x)
  %NONLINEAROBJECTIVE Nonlinear function R^2->R
  
  %% Purpose: 
  % 2-dimensional nonlinear function mapping x -> -0.03/((x(1) +0.25)^2 +(x(2) -0.2)^2 +0.03) -0.1/((x(1) -0.25)^2 +(x(2)+ 0.2)^2 +0.04) + 0.1/(x(1)^2 +x(2)^2 +0.05) +1 +x(1)^2 +x(2)^2 + 1;
  % Has a local maximizing point at approx [-0.0158;0.0126], and a local minimizing point at approx [-0.265;0.212] and a global minimizing point at approx [0.261;-0.209]
  
  %% Input Definition:
  % x: column vector in R^2 (domain space)
  
  %% Output Definition:
  % value: real number, evaluation of nonlinearObjective at x
  % gradient: real column vector in R^2, evaluation of the gradient with respect to x at x
  % hessian: real 2x2 matrix, evaluation of the hessian with respect to x at x
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient,myHessian]=nonlinearObjective([-0.015793;0.012647]);
  % should return
  % myValue=3.0925; myGradient close to zero; myHessian=[-85.299, 16.795; 16.795 -77.739];
  
  %% Input verification:
  if ~isequal(size(x), [2,1])
    error('x has wrong dimension.');    
  end
  
  %% Implementation:
  x1=x(1);
  x2=x(2);
  
  value = -0.03/((x1 +0.25)^2 +(x2 -0.2)^2 +0.03) -0.1/((x1 -0.25)^2 +(x2 + 0.2)^2 +0.04) + 0.1/(x1^2 +x2^2 +0.05) +1 +x1^2 +x2^2 + 1;
  
  if nargout > 1
    dx1 = 2*(x1 +0.25)*0.03/((x1 +0.25)^2 +(x2 -0.2)^2 +0.03)^2 +2*(x1 -0.25)*0.1/((x1 -0.25)^2 +(x2 + 0.2)^2 +0.04)^2 -2*x1*0.1/(x1^2 +x2^2 +0.05)^2 +2*x1;
    dx2 = 2*(x2 -0.2)*0.03/((x1 +0.25)^2 +(x2 -0.2)^2 +0.03)^2 +2*(x2 +0.2)*0.1/((x1 -0.25)^2 +(x2 + 0.2)^2 +0.04)^2 -2*x2*0.1/(x1^2 +x2^2 +0.05)^2 +2*x2;
    
    gradient = [dx1;dx2];
    
  end
  
  if nargout > 2
    
    dx11 = -.6e-1/((x1+.25)^2+(x2-.2)^2+.3e-1)^3*(2*x1+.50)^2+.6e-1/((x1+.25)^2+(x2-.2)^2+.3e-1)^2-.2/((x1-.25)^2+(x2+.2)^2+.4e-1)^3*(2*x1-.50)^2+.2/((x1-.25)^2+(x2+.2)^2+.4e-1)^2+.8/(x1^2+x2^2+.5e-1)^3*x1^2-.2/(x1^2+x2^2+.5e-1)^2+2;
    dx12 = -.6e-1/((x1+.25)^2+(x2-.2)^2+.3e-1)^3*(2*x1+.50)*(2*x2-.4)-.2/((x1-.25)^2+(x2+.2)^2+.4e-1)^3*(2*x1-.50)*(2*x2+.4)+.8/(x1^2+x2^2+.5e-1)^3*x1*x2;
    
    dx22 = -.6e-1/((x1+.25)^2+(x2-.2)^2+.3e-1)^3*(2*x2-.4)^2+.6e-1/((x1+.25)^2+(x2-.2)^2+.3e-1)^2-.2/((x1-.25)^2+(x2+.2)^2+.4e-1)^3*(2*x2+.4)^2+.2/((x1-.25)^2+(x2+.2)^2+.4e-1)^2+.8/(x1^2+x2^2+.5e-1)^3*x2^2-.2/(x1^2+x2^2+.5e-1)^2+2;
    hessian = [dx11, dx12; dx12, dx22];
    %warning("You called evaluation of Hessian of nonlinearObjective. Try to avoid this.")
    
  end
  
end

