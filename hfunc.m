function h = hfunc(V, mu)
    
    vx = V(1);
    vy = V(2);
    vz = V(3);
    
    h = [2*vy; -2*vx; 0];
    
end