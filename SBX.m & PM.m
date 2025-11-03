% SBX.m
function [c1, c2] = SBX(p1, p2, eta_c, VarMin, VarMax)
    nVar = numel(p1);
    c1 = zeros(1, nVar);
    c2 = zeros(1, nVar);
    
    for i = 1:nVar
        u = rand();
        if u <= 0.5
            beta = (2*u)^(1/(eta_c+1));
        else
            beta = (1/(2*(1-u)))^(1/(eta_c+1));
        end
        
        c1(i) = 0.5 * ((1+beta)*p1(i) + (1-beta)*p2(i));
        c2(i) = 0.5 * ((1-beta)*p1(i) + (1+beta)*p2(i));
    end
    
    c1 = max(c1, VarMin); c1 = min(c1, VarMax);
    c2 = max(c2, VarMin); c2 = min(c2, VarMax);
end

% PM.m
function child = PM(parent, pMutation, eta_m, VarMin, VarMax)
    nVar = numel(parent);
    child = parent;
    
    for i = 1:nVar
        if rand() < pMutation
            u = rand();
            if u < 0.5
                delta = (2*u)^(1/(eta_m+1)) - 1;
            else
                delta = 1 - (2*(1-u))^(1/(eta_m+1));
            end
            child(i) = parent(i) + delta * (VarMax(i) - VarMin(i));
        end
    end
    
    child = max(child, VarMin);
    child = min(child, VarMax);
end
