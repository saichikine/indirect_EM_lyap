function L_points = CR3BP_L_Points(mu)
    % Locates equilibrium points of CR3BP in nondimensional coordinates
    % See Grebow 2006 thesis (pg 13,14)
    
    % L1
    L1_x = @(x) x - (1-mu)/(x+mu) + mu/(x-1+mu);
    xL1 = fzero(L1_x, 1-mu-eps);
    L1 = [xL1; 0; 0];
    
    % L2
    L2_x = @(x) x - (1-mu)/(x+mu) - mu/(x-1+mu);
    xL2 = fzero(L2_x, 1-mu+eps);
    L2 = [xL2; 0; 0];
    
    % L3
    L3_x = @(x) x + (1-mu)/(x+mu) + mu/(x-1+mu);
    xL3 = fzero(L3_x, -1);
    L3 = [xL3; 0; 0];
    
    % L4
    L4 = [1/2 - mu; sqrt(3)/2; 0];
    
    % L5
    L5 = [1/2 - mu; -sqrt(3)/2; 0];
    
    L_points = [L1 L2 L3 L4 L5];
    
end