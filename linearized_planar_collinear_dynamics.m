function X = linearized_planar_collinear_dynamics(t, mu, gamma, L_point_num, x0_L_center, s, beta_3)
    
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

    X = [x_approx; y_approx; z_approx; xdot_approx; ydot_approx; zdot_approx];
    
end