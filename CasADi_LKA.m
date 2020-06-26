% Implement tasks that need to be performed only once,
% such as pre-computed constants.

import casadi.*

% define parameters
M = 1575;       % kg, mass
Iz = 2800;      % kg-m2, moment of inertia about vertical axis
lf = 1.2;       % m, distance from center of gravity to front axle
lr = 1.6;       % m, distance from center of gravity to rear axle
cf = 19000;     % N/rad, cornering stiffnesses of front tires
cr = 33000;     % N/rad, cornering stiffnesses of rear tires

% sample time
Ts = 0.1;

% problem size
nx = 4;
ny = 1;
nu = 1;
nc = 1;

% objective function parameters
Q = 1*eye(nx+ny);
R = 10000*eye(nu);

% prediction horizon
N = 10;

% constraints on states
x_min = [-inf ; -inf ; -inf ; -inf ; -inf];
x_max = [ inf ;  inf ;  inf ;  inf ;  inf];
x0 = zeros(nx+ny,1);

% constraints on inputs
u_min = [-0.7];
u_max = [0.7];
u0 = [0];

% preview curvative
CNmeas = zeros(nc*N,1);

% start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "lift" initial conditions
Xk = MX.sym('X0', nx+ny);
w = {w{:}, Xk};
lbw = [lbw; x0];
ubw = [ubw; x0];
w0 = [w0; x0];

% "lift" curvature preview
CurvN = MX.sym('CurvN', nc*N);
w = {w{:}, CurvN};
lbw = [lbw; CNmeas];
ubw = [ubw; CNmeas];
w0 = [w0; CNmeas];

% "lift" velocity Vx
Vx = MX.sym('Vx', 1);
w = {w{:}, Vx};
lbw = [lbw; 15];
ubw = [ubw; 15];
w0 = [w0; 15];

% calculate state space matrices for given bounds on Vx
Vx_min = 10;
Vx_max = 40;
Ac_fcn = @(Vx)[0, 1, 0, 0 ; 0, -(2*cf+2*cr)/(M*Vx), (2*cf+2*cr)/M, (-2*lf*cf+2*cr*lr)/(M*Vx) ; 0, 0, 0, 1 ; 0, (-2*lf*cf+2*lr*cr)/(Iz*Vx), (2*lf*cf-2*lr*cr)/Iz, -(2*lf^2*cf+2*lr^2*cr)/(Iz*Vx)];
B1c_fcn = @(Vx)[0 ; 2*cf/M ; 0 ; 2*lf*cf/Iz];
B2c_fcn = @(Vx)[0 ; (-2*lf*cf+2*lr*cr)/(M*Vx)-Vx ; 0 ; -(2*lf^2*cf+2*lr^2*cr)/(Iz*Vx)];
Ac_min = Ac_fcn(Vx_min);
B1c_min = B1c_fcn(Vx_min);
B2c_min = B2c_fcn(Vx_min);
Ac_max = Ac_fcn(Vx_max);
B1c_max = B1c_fcn(Vx_max);
B2c_max = B2c_fcn(Vx_max);

% discretization of continuous time ss matrices
sysc = ss(Ac_min,[B1c_min,B2c_min],[1,0,0,0],0);
sysd = c2d(sysc,Ts);
A_min = sysd.A;
B1_min = sysd.B(:,1);
B2_min = sysd.B(:,2);
sysc = ss(Ac_max,[B1c_max,B2c_max],[1,0,0,0],0);
sysd = c2d(sysc,Ts);
A_max = sysd.A;
B1_max = sysd.B(:,1);
B2_max = sysd.B(:,2);

% augment matrices
C = [1, 0, 0, 0];
Aaug_min = [A_min, zeros(nx,ny) ; C, eye(ny)];
B1aug_min = [B1_min ; zeros(ny,nu)];
B2aug_min = [B2_min ; zeros(ny,nc)];
Aaug_max = [A_max, zeros(nx,ny) ; C, eye(ny)];
B1aug_max = [B1_max ; zeros(ny,nu)];
B2aug_max = [B2_max ; zeros(ny,nc)];
Caug = [C, zeros(ny,ny)];

% create casadi representation of A, B1, B2
lambda = (Vx-Vx_min)/(Vx_max-Vx_min);
Aaug = Aaug_min*(1-lambda) + Aaug_max*lambda;
B1aug = B1aug_min*(1-lambda) + B1aug_max*lambda;
B2aug = B2aug_min*(1-lambda) + B2aug_max*lambda;

% formulate quadratic program
for k = 0:N-1
    % new NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], nu);
    w = {w{:}, Uk};
    lbw = [lbw; u_min];
    ubw = [ubw; u_max];
    w0 = [w0;  u0];
    
    % add contribution to objective
    J = J + Xk'*Q*Xk + Uk'*R*Uk;
    
    % dynamics
    Xk_end = Aaug*Xk + B1aug*Uk + B2aug*CurvN(1+nc*k:nc+nc*k);
    
    % new NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nx+ny);
    w = {w{:}, Xk};
    lbw = [lbw; x_min];
    ubw = [ubw; x_max];
    w0 = [w0; x0];
    
    % add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; zeros(nx+ny,1)];
    ubg = [ubg; zeros(nx+ny,1)];
end
% add terminal cost
J = J + Xk'*Q*Xk;

% offset vectors
offset_x0 = 0;
offset_curv = nx + ny;
offset_vx = nx + ny + nc*N;
offset_u0 = nx + ny + nc*N + 1;

% create an QP solver
opts.print_iter = 0;
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = qpsol('solver', 'qrqp', prob, opts);

% solve the QP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
w0 = w_opt;

obj.casadi_solver = solver;
obj.q0 = 0;
obj.x0 = w0;
obj.lbx = lbw;
obj.ubx = ubw;
obj.lbg = lbg;
obj.ubg = ubg;
obj.offset_x0 = offset_x0;
obj.offset_curv = offset_curv;
obj.offset_vx = offset_vx;
obj.offset_u0 = offset_u0;
obj.nx = nx+ny;
obj.N = N;
obj.nu = nu;

% Call function for extracting arguments, calculating augmented state,
% updating the bounds, and solving the MPC ==> u = stepImpl(obj,x,r,Vx)

% extract arguments
w0 = obj.x0;
lbw = obj.lbx;
ubw = obj.ubx;
solver = obj.casadi_solver;

% Define MX variables for solver/s-function
s0 = MX.sym('s0',2);
x = MX.sym('x',4);
r = MX.sym('r',1);
Vx = MX.sym('Vx',1);
w0 = MX(w0);
lbw = MX(lbw);
ubw = MX(ubw);

% calculate augmented state and curvature
xaug = [x ; obj.q0 + x(1)];
curv = Vx./r;

% update bounds for measured state
w0(obj.offset_x0+1:obj.offset_x0+obj.nx) = xaug;
lbw(obj.offset_x0+1:obj.offset_x0+obj.nx) = xaug;
ubw(obj.offset_x0+1:obj.offset_x0+obj.nx) = xaug;

% update bounds on curvature preview
w0(obj.offset_curv+1:obj.offset_curv+obj.N) = curv;
lbw(obj.offset_curv+1:obj.offset_curv+obj.N) = curv;
ubw(obj.offset_curv+1:obj.offset_curv+obj.N) = curv;
 fi
% update bounds on velocity
w0(obj.offset_vx+1:obj.offset_vx+1) = Vx;
lbw(obj.offset_vx+1:obj.offset_vx+1) = Vx;
ubw(obj.offset_vx+1:obj.offset_vx+1) = Vx;


% solve mpc problem
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', obj.lbg, 'ubg', obj.ubg);
w_opt = full(sol.x);
u = w_opt(1+obj.offset_u0:obj.nu+obj.offset_u0);
obj.w0 = w_opt; % for warm starting


%% CasADi function with three inputs and one output
% Mapping from initial state to control action
% *No mapping for continuous to discrete

% function_name = 'f';
% f = Function(function_name,{s0},{sol.x(3)});
% 
% file_name = 'f.casadi';
% f.save(file_name);

lib_path = GlobalOptions.getCasadiPath();
inc_path = GlobalOptions.getCasadiIncludePath();
mex('-v',['-I' inc_path],['-L' lib_path],'-lcasadi', 'casadi_fun.c')