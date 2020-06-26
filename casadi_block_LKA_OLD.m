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
        opti
        X
        U
        CurvN
        Vx
        q0
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
            
            % build MPC problem using opti
            opti = casadi.Opti();
            X = cell(N+1,1);
            U = cell(N,1);
            X{1} = opti.parameter(nx+ny);
            opti.set_value(X{1},x0);
            CurvN = opti.parameter(nc*N);
            opti.set_value(CurvN,CNmeas);
            J = 0;
            for k = 1:N
                % inputs
                U{k} = opti.variable(nu);
                opti.subject_to( u_min <= U{k} <= u_max )
                opti.set_initial(U{k}, u0)
                
                % define new state variables
                X{k+1} = opti.variable(nx+ny);
                opti.set_initial(X{k+1}, x0)
                
                % enforce dynamics
                opti.subject_to( X{k+1} == Aaug*X{k} + B1aug*U{k} + B2aug*CurvN(nc*(k-1)+1:nc*k) );
                
                % add contribution to objective
                J = J + X{k}'*Q*X{k} + U{k}'*R*U{k};
            end
            opti.minimize( J )
            
            % define solver
%             opts.qpsol = 'qpoases';
%             opts.tol_pr = 1e-6;
%             opts.tol_du = 1e-6;
%             opts.min_step_size = 1e-10;
%             opti.solver('sqpmethod',opts);
            p_opts = struct('expand',false);
            s_opts = struct('max_iter',1000);
            s_opts.tol = 1e-3;
            s_opts.print_level = 0;
            opts.qpsol = 'qrqp';
            opti.solver('ipopt',p_opts,s_opts);
            sol = opti.solve();
            opti.set_initial(sol.value_variables())
            
            obj.opti = opti;
            obj.X = X;
            obj.U = U;
            obj.CurvN = CurvN;
            obj.Vx = Vx;
            obj.q0 = 0;            
        end
        
        function u = stepImpl(obj,x,r)
            r = 1./r;
            
            % calculate augmented state and curvature
            xaug = [x ; obj.q0 + x(1)];            
            curvN = obj.Vx./r;
            
            % update bounds for measured state
            obj.opti.set_value(obj.X{1},xaug)

            % update bounds on curvature preview
            obj.opti.set_value(obj.CurvN,curvN)        
            
            % integrate error in lateral deviation
            obj.q0 = obj.q0 + x(1);

            % solve mpc problem
            sol = obj.opti.solve();
            
            % extract optimal inputs
            u = sol.value(obj.U{1});
            
            % warm start
            obj.opti.set_initial(sol.value_variables())
        end
        
        function resetImpl(obj)
            % Initialize discrete-state properties.
        end
    end
end
