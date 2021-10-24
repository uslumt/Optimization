function [y] = LLTSolver(L,r,verbose)
  %LLTSOLVER Solves y = (L*L')^(-1)*r.
  
  %% Purpose:
  % LLTSolver solves  (L*L')*y=r for y using forward and backward substitution
  
  %% Input Definition:
  % L: real valued lower triangle matrix nxn with nonzero diagonal elements
  % r: column vector in R^n
  % verbose: bool, if set to true, verbose information is displayed
  
  %% Output Definition:
  % y: column vector in R^n (solution in domain space)
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [y]=LLTSolver([2,0,0;0.5,sqrt(15/4),0;0,0,2],[5;5;4],true);
  % should return
  % y=[1;1;1];
  
  % [y]=LLTSolver([22,0,0,0,0;17,13,0,0,0;13,-2,17,0,0;8,-4,-7,18,0;4,-5,-4,-5,19],[1320;773;1192;132;1405],true);
  % should return
  % y=[1;0;2;0;3];
  
  %% Input verification:
  n=length(r);
  
  if ~isequal(size(L), [n,n])
    error('A has wrong dimension.');
  end
  
  if nargin < 3
    verbose = false;
  end
  
  %% Implementation:
  if verbose
    disp('Start LLTSolver...');
  end
  
  %static
  %none
  
  %dynamic
  
  for i=1:n
    s (i,1)= r(i,1);
    for j=1:i-1
      s(i,1)=s(i,1)-L(i,j)*s(j);
    end
    if (L(i,i)==0)
      disp('Zero diagonal element detected...');
    end
    s(i,1)=s(i,1)/L(i,i);
  end
  
  y=s;
  for i=n:-1:1
    for j=n:-1:i+1
      y(i,1)=y(i,1)-L(j,i)*y(j);
    end
    y(i,1)=y(i,1)/L(i,i);
  end
  
  
  if verbose
    residual = norm((L*L')*y-r)
    disp(sprintf('LLTSolver terminated with norm of residual =%d\n', norm(residual)));
  end
  
end
