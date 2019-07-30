function H = hamiltonian_min_time(chi,params)

    % Computes minimum-time optimal Hamiltonian
    
    %% Unpack
    R = chi(1:3);
    V = chi(4:6);
    m = chi(7);
    lambdaR = chi(8:10);
    lambdaV = chi(11:13);
    if norm(lambdaV) == 0
        lambdaVhat = zeros(3,1);
    else
        lambdaVhat = lambdaV/norm(lambdaV);
    end
    lambdam = chi(14);
    
    mu = params.mu;
    Tmax = params.Tmax;
    c = params.c;
    L_EM = params.L_EM;
    T_EM = params.T_EM;
    m_norm = params.m_norm;
    
    Tmax_normalized = Tmax/(m_norm/1000/L_EM*(T_EM/(2*pi))^2);
    m_normalized = m/m_norm;
    c_normalized = c/L_EM*T_EM/(2*pi);
    
    %% Compute result
    u = opt_control_min_time(chi,params);
%     H = lambdaR'*V + lambdaV'*(gfunc(R) + hfunc(V) - u*Tmax_normalized/m_normalized*lambdaVhat) - lambdam*u*Tmax_normalized/c_normalized + 1;
    H = lambdaR'*V + lambdaV'*(gfunc(R,mu) + hfunc(V,mu) - u*Tmax/m*lambdaVhat) - lambdam*u*Tmax/c + 1;
    
end
    