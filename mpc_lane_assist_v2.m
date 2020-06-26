
clear
import casadi.*

% define parameters
M = 1575;       % kg, mass
Iz = 2875;      % kg-m2, moment of inertia about vertical axis
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
R = 10*eye(nu);

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
Vx_max = 30;
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

% create function for plant model
x = MX.sym('x', nx+ny);
u = MX.sym('u', nu);
curv = MX.sym('curv', nc);
fplant = Function('fplant', {x, u, Vx, curv}, {Aaug*x + B1aug*u + B2aug*curv});

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

% closed-loop simulation
Nsim = 50;
X = zeros(nx+ny,Nsim+1);
U = zeros(nu,Nsim);
Vx_list = 15*ones(1,Nsim);
Curv = Vx_list./10000; Curv(Nsim/2:end) = Vx_list(Nsim/2:end)./10000;
Tmpc = zeros(1,Nsim);
X(:,1) = [0.2 ; 0 ; 0 ; 0 ; 0];
for k = 1:Nsim
    % update bounds for measured state
    w0(offset_x0+1:offset_x0+nx+ny) = X(:,k);
    lbw(offset_x0+1:offset_x0+nx+ny) = X(:,k);
    ubw(offset_x0+1:offset_x0+nx+ny) = X(:,k);
    
    % update bounds on curvature preview
    w0(offset_curv+1:offset_curv+nc*N) = Curv(k)*ones(N,1);
    lbw(offset_curv+1:offset_curv+nc*N) = Curv(k)*ones(N,1);
    ubw(offset_curv+1:offset_curv+nc*N) = Curv(k)*ones(N,1);

    % update bounds on velocity
    w0(offset_vx+1:offset_vx+1) = Vx_list(k);
    lbw(offset_vx+1:offset_vx+1) = Vx_list(k);
    ubw(offset_vx+1:offset_vx+1) = Vx_list(k);

    % solve NMPC problem
    tic
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    Tmpc(k) = toc;
    w_opt = full(sol.x);
    w0 = w_opt; % for warm starting
    
    % extract optimal inputs
    U(:,k) = w_opt(1+offset_u0:nu+offset_u0);
    
    % give information to plant
    X(:,k+1) = full(fplant(X(:,k), U(:,k), Vx_list(k), Curv(:,k)));
end

% time profile
time = (0:Nsim)*Ts;

% plot states
figure; hold on;
plot(time,X(1,:),'-b','linewidth',2)
plot([time(1) time(end)], [x_min(1), x_min(1)],'--k','linewidth',2)
plot([time(1) time(end)], [x_max(1), x_max(1)],'--k','linewidth',2)
plot([time(Nsim/2+1) time(Nsim/2+1)], [-0.2, 1],'--r','linewidth',2)
xlabel('time (seconds)')
ylabel('lateral deviation (m)')
set(gcf,'color','w');
set(gca,'FontSize',20)

% plot inputs
figure; hold on;
stairs(time(1:end-1),U(1,:),'-b','linewidth',2)
plot([time(1) time(end)], [u_min(1), u_min(1)],'--k','linewidth',2)
plot([time(1) time(end)], [u_max(1), u_max(1)],'--k','linewidth',2)
plot([time(Nsim/2+1) time(Nsim/2+1)], [-1, 1],'--r','linewidth',2)
xlabel('time (seconds)')
ylabel('steering angle (rad)')
set(gcf,'color','w');
set(gca,'FontSize',20)
