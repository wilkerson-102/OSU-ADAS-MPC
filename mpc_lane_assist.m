
clear
import casadi.*

% define parameters
M = 1575;       % kg, mass
Iz = 2875;      % kg-m2, moment of inertia about vertical axis
lf = 1.2;       % m, distance from center of gravity to front axle
lr = 1.6;       % m, distance from center of gravity to rear axle
cf = 19000;     % N/rad, cornering stiffnesses of front tires
cr = 33000;     % N/rad, cornering stiffnesses of rear tires
Vx = 15;        % m/s, longitudinal speed

% continuous time system state space matrices
%state vector = [e1, e1dot, e2, e2dot]' where e1 is lateral position wrt road and e2 is yaw angle wrt road
%input vector = [delta]' where delta is front wheel steering angle
%feedforward vector = [phidot_des] where phidot_des is the desired yaw rate determined from curvature of raod = Vx/R
Ac = [0, 1, 0, 0 ; 0, -(2*cf+2*cr)/(M*Vx), (2*cf+2*cr)/M, (-2*lf*cf+2*cr*lr)/(M*Vx) ; 0, 0, 0, 1 ; 0, (-2*lf*cf+2*lr*cr)/(Iz*Vx), (2*lf*cf-2*lr*cr)/Iz, -(2*lf^2*cf+2*lr^2*cr)/(Iz*Vx)]; % state to state mapping
B1c = [0 ; 2*cf/M ; 0 ; 2*lf*cf/Iz]; % input to state mapping
B2c = [0 ; (-2*lf*cf+2*lr*cr)/(M*Vx)-Vx ; 0 ; -(2*lf^2*cf+2*lr^2*cr)/(Iz*Vx)]; % feedforward to state mapping
C = [1, 0, 0, 0];

% system sizes
nx = size(Ac,1);
nu = size(B1c,2);
nc = size(B2c,2);
ny = size(C,1);

% get discrete time state space matrices
Ts = 0.1;
sysc = ss(Ac,[B1c,B2c],C,0);
sysd = c2d(sysc,Ts);
A = sysd.A;
B1 = sysd.B(:,1);
B2 = sysd.B(:,2);
C = sysd.C; 

% augment dynamics
Aaug = [A, zeros(nx,ny) ; C, eye(ny)];
B1aug = [B1 ; zeros(ny,nu)];
B2aug = [B2 ; zeros(ny,nc)];
Caug = [C, zeros(ny,ny)];

% objective function parameters
Q = zeros(nx+ny); Q(1,1) = 10;
R = eye(nu);

% prediction horizon
N = 20;

% constraints on states
x_min = [-inf ; -4 ; -inf ; -inf ; -inf];
x_max = [inf ; 4 ; inf ; inf ; inf];
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

% formulate quadratic program
for k = 0:N-1
    % new NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
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

% store offset for first input
offset_u0 = nx+ny+nc*N;

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = qpsol('solver', 'qpoases', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);
w0 = w_opt;

% closed-loop simulation
Nsim = 50;
X = zeros(nx+ny,Nsim+1);
U = zeros(nu,Nsim);
Curv = 0*ones(nc,Nsim); Curv(Nsim/2+1:end) = Vx/20;
Tmpc = zeros(1,Nsim);
X(:,1) = [1; 0; 0; 0 ; 0];
for k = 1:Nsim
    % update bounds for measured state
    w0(1:nx+ny) = X(:,k);
    lbw(1:nx+ny) = X(:,k);
    ubw(1:nx+ny) = X(:,k);
    
    % update bounds on curvature preview
    w0(nx+ny+1:nx+ny+N*nc) = Curv(k)*ones(N,1);
    lbw(nx+ny+1:nx+ny+N*nc) = Curv(k)*ones(N,1);
    ubw(nx+ny+1:nx+ny+N*nc) = Curv(k)*ones(N,1);
    
    % solve NMPC problem
    tic
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);
    Tmpc(k) = toc;
    w_opt = full(sol.x);
    w0 = w_opt; % for warm starting
    
    % extract optimal inputs
    U(:,k) = w_opt(1+offset_u0:nu+offset_u0);
    
    % give information to plant
    X(:,k+1) = Aaug*X(:,k) + B1aug*U(:,k) + B2aug*Curv(:,k);    
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
