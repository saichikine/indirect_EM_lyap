function orbit_period = linear_planar_L1L2_orbit_period(mu, L_point_num)

    %% Set up constants for sln
    L_points = lagrangePoints(mu);
    L_point = L_points(:,L_point_num);
    
    A_mat = CR3BPLinA(L_point, mu);

    Omega_xx = A_mat(4,1);
    Omega_yy = A_mat(5,2);

    beta_1 = 2 - (Omega_xx + Omega_yy)/2;
    beta_2 = sqrt(-Omega_xx*Omega_yy);

    s = sqrt(beta_1 + sqrt(beta_1^2 + beta_2^2));

    %% Results
    orbit_period = 2*pi/s;
    
end