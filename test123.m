%[xmin]=globalNewtonDescent(@(x)nonlinearObjective(x), [-0.01;0.01], 1.0e-6, true)
% should return
% xmin close to [0.26;-0.21] (exact xmin depends on choice of delta in PrecCGSolver);

[xmin]=globalNewtonDescent(@(x)nonlinearObjective(x), [-0.6;0.6], 1.0e-3, true);
% should return
% xmin close to [-0.26;0.21] (exact xmin depends on choice of delta in PrecCGSolver);

[xmin]=globalNewtonDescent(@(x)nonlinearObjective(x), [0.6;-0.6], 1.0e-3, true);
% should return
% xmin close to [-0.26;0.21] (exact xmin depends on choice of delta in PrecCGSolver);