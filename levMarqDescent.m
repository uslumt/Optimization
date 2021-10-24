function [pmin] = levMarqDescent(R, p0, eps, alpha0, beta, verbose)
  %LEVMARQDESCENT Minimize nonlinear least squares with Levenberg-Marquardt-Algorithm
  
  %% Purpose:
  % Find pmin to satisfy norm(jacobian_R'*R(pmin))<=eps
  
  %% Input Definition:
  % R: error vector handle of type [errorVector, jacobian] = R(p).
  % p0: column vector in R^n (parameter point), starting point.
  % eps: positive value, tolerance for termination. Default value: 1.0e-3.
  % alpha0: positive value, starting value for damping. Default value: 1.0e-3.
  % beta: positive value bigger than 1, scaling factor for alpha. Default value: 100.
  % verbose: bool, if set to true, verbose information is displayed.
  
  %% Output Definition:
  % pmin: column vector in R^n (parameter point), satisfies norm(jacobian_R'*R(pmin))<=eps
  
  %% Required files:
  % [x] = PrecCGSolver(A,b,delta)
  
  %% Test cases:
  % [pmin]=levMarqDescent(@(p)leastSquaresObjective(@simpleValleyModel,p,[0, 0; 1, 2],[2, 3]), [180;0], 1.0e-4, 1.0e-3, 100, true);
  % should return
  % pmin close to [1;1];
  
  % [pmin]=levMarqDescent(@(p)leastSquaresObjective(@exponentialModel,p,[-2, -1, 0, 1, 2],[0.86,0.63,0,-1.72,-6.39]), [0;0;0], 1.0e-4, 1.0e-3, 10, true);
  % should return
  % pmin close to [1;-1;1];
  
  %% Input verification:
  
  try
    [errorVector, jacobian] = R(p0);
  catch
    error('evaluation of function handle failed!'); 
  end
  
  if (eps <= 0)
    error('range of eps is wrong!');    
  end
  if (alpha0 <= 0)
    error('range of alpha0 is wrong!');    
  end
  if (beta <= 1)
    error('range of beta is wrong!');    
  end  
  
  if nargin < 6
    verbose = false;
  end
  
  %% Implementation:
  % Hints: 
  % 1. Remember the connection f = 1/2 R'*R and grad_f = J'*R
  % 2. Use eye(n) to get the unit matrix.
  
  if verbose
    disp('Start levMarqDescent...');
    countIter = 0;
  end
 
  p_k = p0;
  alpha_k = alpha0;

 while norm (jacobian' * R(p_k)) >eps
   
   n = length(p0)
   E = eye(n);
   
   A = jacobian' * jacobian + alpha_k * E;
   b = -jacobian' * R(p_k);
   d_k = PrecCGSolver(A ,b,1.0e-6,true);
    
   if 0.5 * R(p_k + d_k)' * R(p_k + d_k) < 0.5 *  R(p_k)' * R(p_k) 
     p_k = p_k + d_k;
     [errorVector, jacobian] = R(p_k);
     alpha_k = alpha0;
   else
     alpha_k = beta * alpha_k;
   endif
  
   countIter = countIter+1;
endwhile

pmin = p_k
  
  if verbose
    disp(sprintf('levMarqDescent terminated after %i steps with norm of gradient =%d\n',countIter, norm(jacobian'*errorVector)));
  end
end