function [value, gradient, hessian] = quadraticConstraint(x,A,b,c)
  %QUADRATICOBJECTIVE Nonlinear quadratic function R^n->R
  
  %% Purpose: 
  % n-dimensional quadratic function mapping x -> 0.5*x'*A*x + b'*x +c
  
  %% Input Definition:
  % x: column vector in R^n (domain space)
  % A: real valued matrix nxn
  % b: column vector in R^n
  % c: real number
  
  %% Output Definition:
  % value: real number, evaluation of quadraticConstraint at x for parameters A,b,c 
  % gradient: real column vector in R^n, evaluation of the gradient with respect to x at x for parameters A,b,c
  % hessian: real nxn matrix, evaluation of the hessian with respect to x at x for parameters A,b,c
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [myValue,myGradient,myHessian]=quadraticConstraint([1;1],[1,0;0,1],[1;1],1);
  % should return
  % myValue=4; myGradient=[2;2]; myHessian=[1,0;0,1];
  
  %% Input verification:
  [n,m]=size(x);
  
  if ~isequal(size(A), [n,n])
    error('A has wrong dimension.');    
  end
  if ~isequal(size(b), [n,1])
    error('b has wrong dimension.');    
  end
  if ~isequal(size(c), [1,1])
    error('c has wrong dimension.');    
  end
  
  %% Implementation:
  
  value = 0.5*x'*A*x + b'*x +c;
  
  if nargout > 1
    gradient = A*x + b;
    
  end
  
  if nargout > 2
    hessian = A;
    
  end
  
end
