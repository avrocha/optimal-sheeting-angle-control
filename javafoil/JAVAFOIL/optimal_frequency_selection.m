%% GB-ESC 4D
clear

% State dimensions
n  = 4; % Number of dither signals
nc = 9; % Number of control variables

% Maximum frequency
max_f = 0.2; % T = 5s

% Criterion | x = [f_dither(1), f_dither(2), f_dither(3), f_dither(4), y]
f      = -ones(n+nc, 1);
intcon = (n+1):(n+nc);

% Lower/Upper bounds
lb = zeros(n+nc, 1);
ub = [0.7*max_f; 0.8*max_f; 0.9*max_f; max_f; ones(nc,1)]; % choose order according to direction of AWA

% Inequality constraints
tol = 0.1 * max_f; % numerical tolerance to avoid implementation-induced coupling
M   = 10 * max_f; % large constant for integer-programming

% x1 < x2 < x3 < x4
A1 = [1 -1 0 0;
      0 1 -1 0;
      0 0 1 -1];
  
A1 = [A1, zeros(3, nc)];

b1 = -tol*ones(3, 1);

% xi + xj ~= xk, forall i,j,k
A2 = [2 -1 0 0;
      2 0 -1 0;
      2 0 0 -1;
      0 2 -1 0;
      0 2 0 -1
      0 0 2 -1
      1 1 -1 0
      1 0 1 -1
      0 1 1 -1];

A2 = [A2, -M*eye(nc);
      -A2, M*eye(nc)];

b2 = [-tol * ones(nc,1);
      (-tol + M) * ones(nc,1)];

A = [A1; A2];
b = [b1; b2];

x = intlinprog(f, intcon, A, b, [], [], lb, ub)

%% NB-ESC 4D 
clear

% State dimensions
n  = 4; % Number of dither signals
nc = 22; % Number of control variables

% Maximum frequency
max_f = 0.2; % T = 5s

% Criterion | x = [f_dither(1), f_dither(2), f_dither(3), f_dither(4), y]
f      = -ones(n+nc, 1);
intcon = (n+1):(n+nc);

% Lower/Upper bounds
lb = zeros(n+nc, 1);
ub = [0.7*max_f; 0.8*max_f; 0.9*max_f; max_f; ones(nc,1)]; % choose order according to direction of AWA

% Inequality constraints
tol = 0.05 * max_f; % numerical tolerance to avoid implementation-induced coupling
M   = 10 * max_f; % large constant for integer-programming

% x1 < x2 < x3 < x4
A1 = [1 -1 0 0;
      0 1 -1 0;
      0 0 1 -1];

n1 = size(A1, 1);
  
A1 = [A1, zeros(n1, nc)];

b1 = -tol*ones(n1, 1);

% xi + xj ~= xk, forall i,j,k
A2 = [2 -1 0 0;
      2 0 -1 0;
      2 0 0 -1;
      0 2 -1 0;
      0 2 0 -1
      0 0 2 -1
      1 1 -1 0
      1 0 1 -1
      0 1 1 -1];
  
n2 = size(A2, 1);

A2 = [A2, -M*eye(n2), zeros(n2, nc-n2);
      -A2, M*eye(n2), zeros(n2, nc-n2)];

b2 = [-tol * ones(n2,1);
      (-tol + M) * ones(n2, 1)];


% 2xi ~= xj + xk, forall i,j,k distinct
A3 = [-1 2 -1 0;
      -1 0 2 -1;
      0 -1 2 -1];
  
n3 = size(A3, 1);

A3 = [A3, zeros(n3, n2), -M*eye(n3), zeros(n3, nc-n2-n3);
      -A3, zeros(n3, n2), M*eye(n3), zeros(n3, nc-n2-n3)];

b3 = [-tol * ones(n3, 1);
      (-tol + M) * ones(n3, 1)];

% 2xi + xj ~= xk, forall i,j,k distinct
A4 = [2 1 -1 0;
      2 0 1 -1;
      1 2 -1 0;
      1 2 0 -1;
      0 2 1 -1;
      1 0 2 -1;
      0 1 2 -1];

n4 = size(A4, 1);

A4 = [A4, zeros(n4, n2+n3), -M*eye(n4), zeros(n4, nc-n2-n3-n4);
      -A4, zeros(n4, n2+n3), M*eye(n4), zeros(n4, nc-n2-n3-n4)];

b4 = [-tol * ones(n4, 1);
      (-tol + M) * ones(n4, 1)];
  
% % xi ~= xj + xk +- xl, forall i,j,k distinct
A5 = [-1 -1 -1 1;
      1 -1 1 -1;
      1 -1 -1 1];

n5 = size(A5, 1);

A5 = [A5, zeros(n5, n2+n3+n4), -M*eye(n5);
      -A5, zeros(n5, n2+n3+n4), M*eye(n5)];

b5 = [-tol * ones(n5, 1);
      (-tol + M) * ones(n5, 1)];

% As in paper, constraint 2 is relaxed since they have the same phase
A = [A1; A2; A3; A4; A5];
b = [b1; b2; b3; b4; b5];

x = intlinprog(f, intcon, A, b, [], [], lb, ub)