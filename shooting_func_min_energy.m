function WVec = shooting_func_min_energy(T, free_vars, X0, m0, Xf, tf, params)

    % Shooting function for minimum energy case
    % Computes boundary point at tf 
    % W(lambda_i, tf) = [r(tf) - rf; v(tf) - vf; lambdam(tf)]
    % Goal is W(lambda_i) = 0
    % T is so that we can use numjac
    
    %% Unpack and setup
    mu = params.mu;
    Tmax = params.Tmax;
    c = params.c;
    L_EM = params.L_EM;
    T_EM = params.T_EM;
    
    Rf = Xf(1:3);
    Vf = Xf(4:6);
    
    lambda_i = free_vars;
    
    t0 = 0;
    ode_opts = odeset('RelTol',3e-14,'AbsTol',1e-20);
    %% Integrate forward with given lambda_i and tf
    
    chi0 = [X0; m0; lambda_i];
    
    [t_shoot, chi_shoot] = ode113(@(t,chi) state_costate_full_dynamics(t,chi,@opt_control_min_energy,params), [t0, tf], chi0, ode_opts);
    
    %% Compute result
    
    WVec = NaN(7,1);
    
    WVec(1:3) = chi_shoot(end,1:3)' - Rf;
    WVec(4:6) = chi_shoot(end,4:6)' - Vf;
    WVec(7) = chi_shoot(end,14);
    
end