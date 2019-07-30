function results = linear_planar_L1L2_dynamics(t, mu, L_point_num, x0_L_center)

    %% Set up constants for sln
    L_points = lagrangePoints(mu);
    gammas = CR3BP_L1_L2_gammas(mu);
    
    if L_point_num == 1
        gamma = gammas.one;
    elseif L_point_num == 2
    	gamma = gammas.two;
    else
        error("Critical mission failure.")
    end

    L_point = L_points(:,L_point_num);
    

    A_mat = CR3BPLinA(L_point, mu);

    Omega_xx = A_mat(4,1);
    Omega_yy = A_mat(5,2);

    beta_1 = 2 - (Omega_xx + Omega_yy)/2;
    beta_2 = sqrt(-Omega_xx*Omega_yy);

    s = sqrt(beta_1 + sqrt(beta_1^2 + beta_2^2));
    beta_3 = (s^2 + Omega_xx)/(2*s);

    orbit_period = 2*pi/s;
    
    %% Dynamics
    if L_point_num==1
        x_approx = (1-mu) + gamma*(-x0_L_center*cos(s*t) - 1);
        y_approx = beta_3*x0_L_center*sin(s*t)*gamma;
        z_approx = 0;
        xdot_approx = gamma*(x0_L_center*s*sin(s*t));
        ydot_approx = beta_3*x0_L_center*s*cos(s*t)*gamma;
        zdot_approx = 0;
    elseif L_point_num==2
        x_approx = (1-mu) + gamma*(x0_L_center*cos(s*t) + 1);
        y_approx = -beta_3*x0_L_center*sin(s*t)*gamma;
        z_approx = 0;
        xdot_approx = gamma*(-x0_L_center*s*sin(s*t));
        ydot_approx = -beta_3*x0_L_center*s*cos(s*t)*gamma;
        zdot_approx = 0;
    else
        error("this is bad")
    end
    
    %% Output
    results = [x_approx; y_approx; z_approx; xdot_approx; ydot_approx; zdot_approx];
    
end