function [projectedX, activeIndexSet] = projectionInBox(x,a,b,eps)
  %PROJECTIONINBOX Projects x into box constraints defined by a and b and computes eps-active index set
  
  %% Purpose: 
  % if x(i) is bigger or smaller than the bounds, x(i) is set to the closest boundary
  
  %% Input Definition:
  % x: column vector in R^n (domain space)
  % a: column vector in R^n, lower bounds for x, must be smaller than b by at least eps in each component.
  % b: column vector in R^n, upper bounds for x
  % eps: nonnegative value, tolerance for accepting being active. Default value: 1.0e-6.
  
  %% Output Definition:
  % projectedX: column vector in R^n, satisfies box constraints 
  % activeIndexSet: row vector with indices, collected indices mark x(i) components with projectedX(i)-a(i) <= eps or projectedX(i) - b(i) >= -eps 
  
  %% Required files:
  % <none>
  
  %% Test cases:
  % [projectedX, activeIndexSet] = projectionInBox([1;1;1;1;1],[0;2;0.9;0;0],[2;3;3;0.5;1.1],0.2);
  % should return
  % projectedX=[1;2;1;0.5;1]; activeIndexSet = [2 3 4 5];
  %
  % [projectedX, activeIndexSet] = projectionInBox([1;1;1;-0.5],[0;0;0.9;-1],[2;3;3;-0.25],0.2);
  % should return
  % projectedX=[1;1;1;-0.5]; activeIndexSet = [3];
  
  %% Input verification:
  [n,m] = size(x);
  
  if ~isequal(size(a), [n,m])
    error('Size of a is wrong.');    
  end
  
  if ~isequal(size(b), [n,m])
    error('Size of b is wrong.');    
  end
  
  if min(b-a)<eps
    error('a and b forming box is degenerate.');    
  end
  
  %% Implementation:
  n=length(x);  
  projectedX = x;
  activeIndexSet = [];
  
  for i=1:n 
    if((x(i) <= a(i)+eps) || (x(i) >= b(i)-eps) )
    activeIndexSet=[activeIndexSet i];
      if (x(i) < a(i))
        projectedX(i) = a(i);
      end
      
      if (x(i) > b(i))
        projectedX(i) = b(i);
      end
    end
  end

end
