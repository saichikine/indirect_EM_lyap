function lyapunov_approx_orbit = lyapunov_approx_ICs(Ax_km, L_km, mu, L_point_string, L_point, gamma)
    % Use Szebehely's approximation to produce ICs for Lyapunov orbit
    
    x0_L_center = Ax_km/L_km;
    
    A_mat = CR3BPLinA(L_point, mu);
    
    Omega_xx = A_mat(4,1);
    Omega_yy = A_mat(5,2);
    
    beta_1 = 2 - (Omega_xx + Omega_yy)/2;
    beta_2 = sqrt(-Omega_xx*Omega_yy);
    
    s = sqrt(beta_1 + sqrt(beta_1^2 + beta_2^2));
    
    beta_3 = (s^2 + Omega_xx)/(2*s);
    
    ydot0_L_center = -beta_3*x0_L_center*s;
    
    orbit_period = 2*pi/s;
    
    % Transform back into usual CR3BP coordinate system
    if strcmpi(L_point_string,"L1")
        x0 = (1-mu) + gamma*(x0_L_center - 1);
    elseif strcmpi(L_point_string,"L2")
        x0 = (1-mu) + gamma*(x0_L_center + 1);
    else
        error('This is bad.')
    end
    ydot0 = ydot0_L_center*gamma;

    X0 = [x0; 0; 0; 0; ydot0; 0];
    
    lyapunov_approx_orbit = struct('ICs',X0,'orbit_period',orbit_period);
end