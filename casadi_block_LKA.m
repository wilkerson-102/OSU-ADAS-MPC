classdef casadi_block_LKA < matlab.System & matlab.system.mixin.Propagates
    % untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    properties
        % Public, tunable properties.
    end

    properties (DiscreteState)
    end

    properties (Access = private)
        % Pre-computed constants.
        casadi_solver
        q0
        x0
        lbx
        ubx
        lbg
        ubg
        Vx
    end

    methods (Access = protected)
        function num = getNumInputsImpl(~)
            num = 2;
        end
        function num = getNumOutputsImpl(~)
            num = 1;
        end
        function dt1 = getOutputDataTypeImpl(~)
        	dt1 = 'double';
        end
        function dt1 = getInputDataTypeImpl(~)
        	dt1 = 'double';
        end
        function sz1 = getOutputSizeImpl(~)
        	sz1 = [1,1];
        end
        function sz1 = getInputSizeImpl(~)
        	sz1 = [1,1];
        end
        function cp1 = isInputComplexImpl(~)
        	cp1 = false;
        end
        function cp1 = isOutputComplexImpl(~)
        	cp1 = false;
        end
        function fz1 = isInputFixedSizeImpl(~)
        	fz1 = true;
        end
        function fz1 = isOutputFixedSizeImpl(~)
        	fz1 = true;
        end
        function setupImpl(obj,~,~)
            % Implement tasks that need to be performed only once, 
            % such as pre-computed constants.
            
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
            R = 1*eye(nu);
            
            % prediction horizon
            N = 10;
            
            % constraints on states
            x_min = [-inf ; -inf ; -inf ; -inf ; -inf];
            x_max = [inf ; inf ; inf ; inf ; inf];
            x0 = zeros(nx+ny,1);
            
            % constraints on inputs
            u_min = [-0.7];
            u_max = [0.7];
            u0 = [0];
            
            % preview curvative
            CNmeas = zeros(nc,1);
            
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
            CurvN = MX.sym('CurvN', nc);
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
                Xk_end = Aaug*Xk + B1aug*Uk + B2aug*CurvN;
                
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
            
%             % store offset for first input
%             offset_u0 = nx+ny+nc;
            
            % Create an NLP solver
            prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
            solver = qpsol('solver', 'qpoases', prob);
            
            obj.casadi_solver = solver;
            obj.q0 = 0;
            obj.x0 = w0;
            obj.lbx = lbw;
            obj.ubx = ubw;
            obj.lbg = lbg;
            obj.ubg = ubg;
            obj.Vx = Vx;
        end
        
        function u = stepImpl(obj,x,r)
            % radius update
            r = 1./r(1);
            
            % extract arguments
            w0 = obj.x0;
            lbw = obj.lbx;
            ubw = obj.ubx;
            solver = obj.casadi_solver;
            
            % calculate augmented state and curvature
            xaug = [x ; obj.q0 + x(1)];            
            curvN = obj.Vx./r;

            % integrate error in lateral deviation
            obj.q0 = obj.q0 + x(1);            
            
            % update bounds for measured state
            w0(1:5) = xaug;
            lbw(1:5) = xaug;
            ubw(1:5) = xaug;
            
            % update bounds on curvature preview
            w0(5+1:5+1) = curvN;
            lbw(5+1:5+1) = curvN;
            ubw(5+1:5+1) = curvN;
            
            % solve NMPC problem
            sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', obj.lbg, 'ubg', obj.ubg);
            w_opt = full(sol.x);
            u = w_opt(6+1:6+1);
            obj.x0 = w_opt; % for warm starting            
        end
        
        function resetImpl(obj)
            % Initialize discrete-state properties.
        end
    end
end
