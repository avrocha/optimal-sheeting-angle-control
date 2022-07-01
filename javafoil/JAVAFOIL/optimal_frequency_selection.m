%% GB-ESC 4D
clear

% Maximum frequency
max_f = 0.2; % T = 5s

% Criterion (f < 0 -> Maximization problem)
% x = [f_dither(1), f_dither(2), f_dither(3), f_dither(4), y]
f      = -ones(4+9, 1);
intcon = 5:(4+9);

% Lower/Upper bounds
lb = zeros(5, 1);
ub = [0.7*max_f; 0.8*max_f; 0.9*max_f; max_f; ones(9,1)]; % choose order according to direction of AWA

% Inequality constraints
tol = 0.1 * max_f; % numerical tolerance to avoid implementation-induced coupling
M   = 10 * max_f; % large constant for integer-programming

% x1 < x2 < x3 < x4
A1 = [1 -1 0 0 0;
      0 1 -1 0 0;
      0 0 1 -1 0];

b1 = [-tol;
      -tol;
      -tol];

% xi + xj ~= xk, forall i,j,k

A2  = [2 -1 0 0 -M;
       -2 1 0 0 M;
       2 0 -1 0 -M;
       -2 0 1 0 M;
       2 0 0 -1 -M;
       -2 0 0 1 M;
       0 2 -1 0 -M;
       0 -2 1 0 M;
       0 2 0 -1 -M;
       0 -2 0 1 M;
       0 0 2 -1 -M;
       0 0 -2 1 M;
       1 1 -1 0 -M;
       -1 -1 1 0 M;
       1 0 1 -1 -M;
       -1 0 -1 1 M;
       0 1 1 -1 -M;
       0 -1 -1 1 M];

b2 = [-tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M];

A = [A1; A2];
b = [b1; b2];

x = intlinprog(f, intcon, A, b, [], [], lb, ub)

%% NB-ESC 4D [WIP]
clear

% Maximum frequency
max_f = 0.2; % T = 5s

% Criterion (f < 0 -> Maximization problem)
% x = [f_dither(1), f_dither(2), f_dither(3), f_dither(4), y]
f      = -ones(5, 1);
intcon = 5;

% Lower/Upper bounds
lb = zeros(5, 1);
ub = [0.7*max_f; 0.8*max_f; 0.9*max_f; max_f; 1]; % choose order according to direction of AWA

% Inequality constraints
tol = 0.1 * max_f; % numerical tolerance to avoid implementation-induced coupling
M   = 100 * max_f; % large constant for integer-programming

% x1 < x2 < x3 < x4
A1 = [1 -1 0 0 0;
      0 1 -1 0 0;
      0 0 1 -1 0];

b1 = [-tol;
      -tol;
      -tol];

% xi + xj ~= xk, forall i,j,k
A2  = [2 -1 0 0 -M;
       -2 1 0 0 M;
       2 0 -1 0 -M;
       -2 0 1 0 M;
       2 0 0 -1 -M;
       -2 0 0 1 M;
       0 2 -1 0 -M;
       0 -2 1 0 M;
       0 2 0 -1 -M;
       0 -2 0 1 M;
       0 0 2 -1 -M;
       0 0 -2 1 M;
       1 1 -1 0 -M;
       -1 -1 1 0 M;
       1 0 1 -1 -M;
       -1 0 -1 1 M;
       0 1 1 -1 -M;
       0 -1 -1 1 M];

b2 = [-tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M];

% 2xi ~= xj + xk, forall i,j,k distinct
A3 = [-1 2 -1 0 -M;
      1 -2 1 0 M;
      -1 0 2 -1 -M;
      1 0 -2 1 M;
      0 -1 2 -1 -M;
      0 1 -2 1 M];

b3 = [-tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M];

% 2xi + xj ~= xk, forall i,j,k distinct
A4 = [2 1 -1 0 -M;
      -2 -1 1 0 M;
      2 0 1 -1 -M;
      -2 0 -1 1 M;
      1 2 -1 0 -M;
      -1 -2 1 0 M;
      1 2 0 -1 -M;
      -1 -2 0 1 M;
      0 2 1 -1 -M;
      0 -2 -1 1 M;
      1 0 2 -1 -M;
      -1 0 -2 1 M;
      0 1 2 -1 -M
      0 -1 -2 1 M];

b4 = [-tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M];

% xi ~= xj + xk +- xl, forall i,j,k distinct
A5 = [-1 -1 -1 1 -M;
      1 1 1 -1 M;
      1 -1 1 -1 -M;
      -1 1 -1 1 M;
      1 -1 -1 1 -M;
      -1 1 1 -1 M];

b5 = [-tol;
      -tol + M;
      -tol;
      -tol + M;
      -tol;
      -tol + M];

% As in paper, constraint 2 is relaxed since they have the same phase
A = [A1; A3; A4; A5];
b = [b1; b3; b4; b5];

x = intlinprog(f, intcon, A, b, [], [], lb, ub)