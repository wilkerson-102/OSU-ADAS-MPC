
% Implement tasks that need to be performed only once,
% such as pre-computed constants.
import casadi.*

% define parameters
r = 1;
tauh = 1;
vf_high = 40;
vf_low = 10;
vfmean = (vf_high+vf_low)/2;
TL = 0.45;
KL = 1.0;
d0 = 10;

% state space matrices
Ac_high = [0, 1, -tauh-r*(2*vf_high-vfmean) ; 0, 0, -1 ; 0, 0, -1/TL];
Ac_low = [0, 1, -tauh-r*(2*vf_low-vfmean) ; 0, 0, -1 ; 0, 0, -1/TL];
B1c = [0 ; 0 ; KL/TL];
B2c = [0 ; 1 ; 0];

% system sizes
nx = size(Ac_high,1);
nu = size(B1c,2);
nc = size(B2c,2);

% get discrete time state space matrices
Ts = 0.1;
sysc_high = ss(Ac_high,[B1c,B2c],eye(nx),0);
sysd_high = c2d(sysc_high,Ts);
sysc_low = ss(Ac_low,[B1c,B2c],eye(nx),0);
sysd_low = c2d(sysc_low,Ts);
A_high = sysd_high.A;
A_low = sysd_low.A;
B1 = sysd_high.B(:,1);
B2 = sysd_high.B(:,2);

% objective function parameters
wdd = 0.02;
wdv = 0.025;
waf = 1;
wafdes = 5;

% constraint parameters
TTC = -2.5;
ds0 = 10;
u_min = -2;
u_max = 2;

% prediction horizon
N = 20;

% initial condition
x0 = [-10 ; 0 ; 0];

% initial guess
x_init = x0;
u_init = 0;

% build MPC problem using opti
opti = casadi.Opti();
X = cell(N+1,1);
U = cell(N,1);
Eps = cell(N,1);
X{1} = opti.parameter(nx);
opti.set_value(X{1},x0);
Vf = opti.parameter();
opti.set_value(Vf,vf_high);
AfRef = opti.parameter();
opti.set_value(AfRef,0);
Ap = opti.parameter();
opti.set_value(Ap,0);
J = 0;
Lambda = (Vf-vf_low)/(vf_high-vf_low);
% Lambda = 1;
for k = 1:N
    % inputs
    U{k} = opti.variable(nu);
    opti.subject_to( u_min <= U{k} <= u_max )
    opti.set_initial(U{k}, u_init)
    
    % define new state variables
    X{k+1} = opti.variable(nx);
    opti.set_initial(X{k+1}, x_init)
    
    % define penalty function variables
    Eps{k} = opti.variable();
    opti.subject_to( 0 <= Eps{k} )
    opti.set_initial(Eps{k}, 0)
    
    % enforce dynamics
    A = Lambda*A_high + (1-Lambda)*A_low;
    opti.subject_to( X{k+1} == A*X{k} + B1*U{k} + B2*Ap )
    
    % safety collision avoidance constraint
    dsafe = TTC*X{k}(2) + ds0;
    ddes = r*(Vf-vfmean)+tauh*Vf+d0;
    opti.subject_to( X{k}(1) + ddes + Eps{k} >= dsafe )
    
    % objective function
    J = J + 1000*Eps{k}^2 + wdd*(X{k}(1))^2 + wdv*(X{k}(2))^2 + waf*(X{k}(3)-AfRef)^2 + wafdes*(U{k})^2;
end
opti.minimize( J )

p_opts = struct('expand',false);
s_opts = struct('max_iter',100);
s_opts.tol = 1e-4;
s_opts.print_level = 0;
opts.qpsol = 'qrqp';
opti.solver('ipopt',p_opts,s_opts);
sol = opti.solve();
opti.set_initial(sol.value_variables());

ACC = opti.to_function('ACC',{X{1}, Vf, AfRef, Ap},{U{1}});

file_name = 'ACC.casadi';
ACC.save(file_name);

lib_path = GlobalOptions.getCasadiPath();
inc_path = GlobalOptions.getCasadiIncludePath();
mex('-v',['-I' inc_path],['-L' lib_path],'-lcasadi', 'casadi_fun.c')
