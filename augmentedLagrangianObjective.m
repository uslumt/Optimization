function [A_Value, A_Gradient,A_Hessian] = augmentedLagrangianObjective(f, h, x, alpha, gamma)
  %AUGMENTEDLAGRANGIANOBJECTIVE builds the augmented Lagrangian out of the objective and one constraint function
  
  %% Purpose: 
  % provides value, gradient and hessian of the augmented Lagrangian mapping x -> A = f(x) + alpha*h(x)+ 0.5*gamma*h(x)^2
  
  %% Input Definition:
  % f: function handle of type [fvalue, fgradient, fHessian] = f(x), objective function.
  % h: function handle of type [hvalue, hgradient, hHessian] = h(x), equality constraint.
  % x: column vector in R^n (domain point).
  % alpha: real value, current guess for Lagrangian multiplier for h.
  % gamma: positive value, penalty parameter.
  
  %% Output Definition:
  % A_value: real number, evaluation of augmentedLagrangianObjective at x 
  % A_gradient: real column vector in R^n, evaluation of the gradient with respect to x at x
  % A_Hessian: real matrix in R^nxn, evaluation of the Hessian with respect to x at x
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient,myHessian]=augmentedLagrangianObjective(@(x)quadraticConstraint(x,[2,0;0,2],[0;0],1), @(x)quadraticConstraint(x,[2,0;0,2],[0;0],-1), [2;2], -1, 10);
  % should return
  % myValue=247 and myGradient=[280;280] and myHessian = [300, 160; 160, 300]
  
  %% Input verification:
  try
    [fvalue, fgradient, fhessian] = f(x);
    [hvalue, hgradient, hhessian] = h(x);
  catch
    error('evaluation of function handle or constraint handle failed!'); 
  end
  
  if (gamma <= 0)
    error('range of gamma is wrong!');    
  end
  
  %% Implementation:  
  % Hints:
  % 1. nabla(h(x)^2) is 2*h(x)*nabla(h(x))
  % 2. nabla^2(h(x)^2) is nabla(2*h(x)*nabla(h(x))) is 2*h(x)*nabla^2(h(x))+2*nabla(h(x))*nabla(h(x))'
   
   A_Value = f(x) + alpha*h(x)+ 0.5*gamma*h(x)^2
   A_Gradient = fgradient + alpha * hgradient + gamma * h(x) * hgradient
   A_Hessian = fhessian + alpha*hhessian + 0.5*gamma*( 2*h(x)*hhessian + 2*(hgradient*hgradient'))
  
end
