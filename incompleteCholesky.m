function [L] = incompleteCholesky(A,lambda,delta,verbose)
  %INCOMPLETECHOLESKY approximates A=L*L^T.
  
  %% Purpose: 
  % incompleteCholesky finds lower triangle matrix L such that A - L*L^T is small, but eigenvalues are positive and sparsity is preserved
  
  %% Input Definition:
  % A: real valued symmetric matrix nxn
  % lambda: nonnegative scalar, lower bound for eigenvalues of L*L^T. Default value: 1.0e-3.
  % delta: scalar, if positive it is tolerance for recognizing nonsparse entry. If negative, do complete cholesky. Default value: 1.0e-6.
  % verbose: bool, if set to true, verbose information is displayed
  
  %% Output Definition:
  % L: real valued lower triangle matrix nxn
  
  %% Required files:
  % <none>
  
  %% Test cases:
  %[L]=incompleteCholesky([5,4,3,2,1;4,5,2,1,0;3,2,5,0,0;2,1,0,5,0;1,0,0,0,5],0,-1,true);
  % executes complete Cholesky decomposition with 
  % norm of residual approx. 8.89e-16
  % the warning is okay.
  
  % [L]=incompleteCholesky([4,1,0;1,4,0;0,0,4],1.0e-3,1.0e-6,true);
  % should return approximately
  % L=[2 0 0;0.5 1.94 0;0 0 2];
  
  % [L]=incompleteCholesky([4,1,0;1,4,0;0,0,4],4,1.0e-6,true);
  % should return approximately
  % L=[2 0 0;0.5 2 0;0 0 2];
  % (there was a mistake in the previous version of this file)
  
  % [L]=incompleteCholesky([4,1,0;1,4,0;0,0,4],1.0e-3,1,true);
  % should return approximately
  % L=[2 0 0;0 2 0;0 0 2];
  
  %% Input verification:
  [n,m]=size(A);

  if ~isequal(n,m)
    error('A has wrong dimension.');    
  end
  
  if ~issymmetric(A)
    error('Matrix is not symmetric.');    
  end
  
  if (lambda < 0)
    error('range of lambda is wrong!');    
  end 
  
  if (delta < 0)
    warning('negative delta detected, sparsity is not preserved.');    
  end 
  
  if nargin < 4
    verbose = false;
  end
  
  %% Implementation:
  if verbose
    disp('Start incompleteCholesky...');
    A_old=A;
  end
  
  %static
  %none
  
  %dynamic
  
	for k=1:n
		A(k,k) = sqrt(max(A(k,k), lambda));
		for i=(k+1):n
		    if (A(i,k) > delta)
		        A(i,k) = A(i,k)/A(k,k);
        else  A(i,k) = 0;
      endif
		endfor
		for j=(k+1):n
		    for i=j:n
		        if (abs(A(i,j)) > delta)
		            A(i,j) = A(i,j)-A(i,k)*A(j,k);  
		        endif
		    endfor
		endfor
	endfor

  for i=1:n
        for j=i+1:n
            A(i,j) = 0;
        endfor
   endfor    
   L=A;
  
  if verbose
     
    residual=norm(A_old-L*L');
    disp(sprintf('incompleteCholesky terminated with norm of residual =%d\n',residual));
  end
  
end
