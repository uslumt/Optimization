function [xmin, lambda] = augmentedLagrangianDescent(f, P, h, x0, alpha0, eps, delta, verbose)
  %AUGMENTEDLAGRANGIANDESCENT Find minimal point of function subject to box constraints and equality constraints by minimizing the augmented Lagrangian
  
  %% Purpose:
  % Find xmin to satisfy the stationarity condition for the augmented Lagrangian: norm(xmin - P(xmin - gradA(xmin)))<=eps
  % And xmin also satisfies the feasibility condition norm(h(xmin))<=delta. 
  
  %% Input Definition:
  % f: function handle of type [fvalue, fgradient] = f(x), objective function.
  % P: box projection handle of type [projectedX, activeIndexSet] = P(x).
  % h: function handle of type [hvalue, hgradient] = h(x), equality constraint.
  % x0: column vector in R^n (domain point), starting point.
  % alpha0: real value, starting guess for Lagrangian multiplier for h. Default value: 0.
  % eps: positive value, tolerance for termination. Default value: 1.0e-3.
  % delta: positive value in (0,eps), tolerance for feasibility. Default value: 1.0e-6.
  % verbose: bool, if set to true, verbose information is displayed.
  
  %% Output Definition:
  % xmin: column vector in R^n (domain point), satisfies norm(xmin - P(xmin - gradA(xmin)))<=eps and norm(h(xmin))<=delta.
  % lambda: real value, approximates Lagrangian multiplier.
  
  %% Required files:
  % [xmin_k] = projectedNewtonDescent(f, P, x0, eps)
  % [A_value, A_gradient, A_hessian] = augmentedLagrangianObjective(f, h, x, alpha, gamma)
  
  %% Test cases:
  % [xmin, lambda] = augmentedLagrangianDescent(@(x)quadraticConstraint(x,[4,0;0,2],[0;0],1), @(x)projectionInBox(x,[0;0],[2;2],1.0e-6), @(x)quadraticConstraint(x,[2,0;0,2],[0;0],-1), [1;1], 0, 1.0e-3, 1.0e-6, true);
  % should return
  % xmin close to [0;1] and lambda close to -1;
  
  %% Input verification:
  
  try
    [fvalue, fgradient, fhessian] = f(P(x0));
    [hvalue, hgradient, hhessian] = h(P(x0));
  catch
    error('evaluation of function handle or projection handle or constraint handle failed!'); 
  end
  
  if (eps <= 0)
    error('range of eps is wrong!');    
  end
  
  if (delta <= 0 || delta >= eps)
    error('range of delta is wrong!');    
  end  
  
  if nargin < 8
    verbose = false;
  end
  
  %% Implementation:
  % Hints:
  % 1. do not confuse eps_k (tolerance for projectedNewtonDescent) and eps (tolerance for outer loop)
  % 2. do not forget to update all depended variables (like Agradient) when
  % x changes and use augmentedLagrangianObjective after the changes of
  % alphak, gammak to update the augmented Lagrangian Ak.
   
  
  
  if verbose
    disp('Start augmentedLagrangianDescent...');
    countIter = 0;
  end
     
  %dynamic
  xk = P(x0);
  hk = h(xk);
  alphak = alpha0;
  gammak = 10;
  Ak = @(x)augmentedLagrangianObjective(f, h, x, alphak, gammak);
  [Avalue, Agradient] = Ak(xk);
  epsk = 1 / gammak;
  deltak = epsk^0.1;
 
 while norm(xk - P(xk - Agradient)) > eps || norm (hk) > delta 
    xk = projectedNewtonDescent(Ak, P, xk, epsk, verbose);
    hk = h(xk);
    [Avalue, Agradient] = Ak(xk);
    if (norm (hk) <= deltak)
      alphak = alphak + (gammak * hk);
      epsk = max(epsk / gammak ,eps);
      deltak = max(deltak / gammak, delta);
    else
      gammak = max(10, sqrt(gammak)) * gammak;
      epsk = 1 / gammak;
      deltak = 1 / gammak;
    end
    Ak = @(x)augmentedLagrangianObjective(f, h, x, alphak, gammak);
    [Avalue, Agradient] = Ak(xk);
    countIter = countIter+1;
 end
    
  xmin=xk;
  lambda=alphak;
 
  if verbose
    disp(sprintf('augmentedLagrangianDescent terminated after %i steps with stationarity =%d and feasibility =%d\n',countIter, norm(xk - P(xk - Agradient)), norm(hk)));
  end
 
end