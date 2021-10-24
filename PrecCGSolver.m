function [x] = PrecCGSolver(A,b,delta,verbose)
  %PrecCGSOLVER Solves linear system A*y = b, only works reliably, if A is symmetric positive or negative definite.
  
  %% Purpose: 
  % CG Solver finds y such that norm(A*y - b) <= delta using incompleteCholesky as preconditioner
  
  %% Input Definition:
  % A: real valued matrix nxn
  % b: column vector in R^n
  % delta: positive value, tolerance for termination. Default value: 1.0e-6.
  % verbose: bool, if set to true, verbose information is displayed
  
  %% Output Definition:
  % x: column vector in R^n (solution in domain space)
  
  %% Required files:
  % [L] = incompleteCholesky(A,1.0e-3,delta)
  % [y] = LLTSolver(L,r)
  
  %% Test cases:
  % [x]=PrecCGSolver([4,1,0;1,7,0;0,0,3],[5;8;3],1.0e-6,true);
  % should return
  % x=[1;1;1];
  
  % [x]=PrecCGSolver([484,374,286,176,88;374,458,195,84,3;286,195,462,-7,-6;176,84,-7,453,-10;88,3,-6,-10,443],[1320;773;1192;132;1405],1.0e-6,true);
  % should return approx
  % x=[1;0;2;0;3];
  
  % [x]=PrecCGSolver([1,2,3,4;2,4,-100,-100;3,-100,7,0;4,-100,0,3],[0;5;8;3],1.0e-6,true);
  % should require MORE than n=4 steps (because the matrix is not SPD)!
  
  %% Input verification:
  [n,m]=size(b);
  
  if ~isequal(size(A), [n,n])
    error('Matrix has wrong dimension.');    
  end
  
  if ~issymmetric(A)
    error('Matrix is not symmetric.');    
  end
  
  if (delta <= 0)
    error('range of delta is wrong!');    
  end 
  
  if nargin < 4
    verbose = false;
  end
  
  %% Implementation:
  if verbose
    disp('Start PrecCGSolver...');
    countIter = 0;
  end
  
  %static

  L = incompleteCholesky(A,1.0e-3,delta);

    x = LLTSolver(L, b);
    r = A* x  - b;
    d= -LLTSolver(L , r );

  while norm(r) > delta 
    Ad= A*d;
    ro = d'*Ad;
    t= r' *LLTSolver(L ,r )/ ro;
    x= x+t*d;
    rold=r;
    r= rold + t*Ad;
    beta = r'*LLTSolver(L , r)/( rold'*LLTSolver(L, rold) );
    d = - LLTSolver(L , r) + beta*d;
    couter = countIter++
  endwhile

 
  if verbose
    disp(sprintf('precCGSolver terminated after %i steps with norm of residual =%d\n',countIter, norm(A*x-b)));
  end
  
end
