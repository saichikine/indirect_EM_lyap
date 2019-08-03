function chidot = state_costate_full_dynamics(t,chi,u_func,params)
    
    % dynamics for 14-dimensional state+costate
    % state is R,V,m (mass)
    % c is exhaust velocity (c = Isp*g0)
    
    % Check for length
    if length(chi) ~= 14
        error("State vector chi is wrong length, should be 14 dimensional.")
    end
    
    %% Unpack
    
    mu = params.mu;
    Tmax = params.Tmax; % nondim
    c = params.c; % nondim
    L_EM = params.L_EM;
    T_EM = params.T_EM;
    
    chi = reshape(chi,[],1); % make column vector
    chidot = NaN(14,1);
    
    R = chi(1:3);
    x = R(1);
    y = R(2);
    z = R(3);
    V = chi(4:6);
    vx = V(1);
    vy = V(2);
    vz = V(3);
    m = chi(7);
    
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x+mu-1)^2 + y^2 + z^2);
    
    lambdaR = chi(8:10);
    lambdaV = chi(11:13);
    if norm(lambdaV)==0
        lambdaVhat = zeros(3,1);
    else
        lambdaVhat = lambdaV/norm(lambdaV);
    end
    lambdam = chi(14);
    
    % g(R), part of CR3BP dynamics that depend on position
    g = [x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
        y - (1-mu)*y/r1^3 - mu*y/r2^3;
        -(1-mu)*z/r1^3 - mu*z/r2^3];
    
    % h(R), part of CR3BP dynamics that depend on velocity
    h = [2*vy; -2*vx; 0];
    
    % G, =dgdR
    G = zeros(3,3);
    G(1,1) = 1 - (1-mu)/r1^3 + 3*(1-mu)*(x+mu)^2/r1^5 - mu/r2^3 + 3*mu*(x+mu-1)^2/r2^5;
    G(2,2) = 1 - (1-mu)/r1^3 + 3*(1-mu)*y^2/r1^5 - mu/r2^3 + 3*mu*y^2/r2^5;
    G(3,3) = -(1-mu)/r1^3 + 3*(1-mu)*z^2/r1^5 - mu/r2^3 + 3*mu*z^2/r2^5;
    G(1,2) = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5;
    G(2,1) = G(1,2);
    G(1,3) = 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x+mu-1)*z/r2^5;
    G(3,1) = G(1,3);
    G(2,3) = 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5;
    G(3,2) = G(2,3);
    
    % H, = dh/dV
    H = [0, 2, 0; -2, 0, 0; 0, 0, 0];
    
    % Set what u is
    u = u_func(chi,params);
    %fprintf("u = %d\n",u);
    
    control_accel = -lambdaVhat*u*Tmax/m;
    %control_accel = -lambdaVhat*u*params.amax*m/params.m0;
    %control_accel = -lambdaVhat*u*params.amax*params.m0/m;
    %fprintf("Magnitude of control accel is: %d\n",norm(control_accel));
    
    %u = 0;
    
    % Build dynamics function
    chidot(1:3) = V; % Rdot
%     chidot(4:6) = g + h - (lambdaVhat*u*Tmax/m)/params.accel_norm;
%     chidot(4:6) = g + h - (lambdaVhat*u*Tmax/m)/1000/L_EM*(T_EM/(2*pi))^2; % Vdot (nondimensionalize control acceleration)
    chidot(4:6) = g + h + control_accel;
    chidot(7) = -u*Tmax/(c); % mdot
    
    chidot(8:10) = -G'*lambdaV;
    chidot(11:13) = -lambdaR - H'*lambdaV;
    chidot(14) = -norm(lambdaV)*u*Tmax/m^2;
    %chidot(14) = -norm(lambdaV)*u*params.amax*m/params.m0/m;
    %chidot(14) = -norm(lambdaV)*u*params.amax*params.m0/m^2;
    %chidot(14) = -norm(lambdaV)*u*Tmax/(m/params.m0)^2;
    
    if any(isnan(chidot))
        error("Invalid dynamics.")
    end
end