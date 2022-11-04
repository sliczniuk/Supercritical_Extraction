function y = modelSFE_out(x, theta)

% theta = [nstages, C0solid, V, epsi, dp, L, rho_s, km, mi];
%          1        2        3  4     5   6  7      8   9
    nstages = theta(1);
    C0solid = theta(2);
    V       = theta(3);
    
    Csolid = x(1 + nstages : 2 * nstages,:);    % Extract concentration is solid
                                              % Each stage
    
    stageV = V / nstages;                     % Volume of one stage
    
    Mextracted = stageV * (C0solid - Csolid); % Mass of extract extracted
                                              % Each stage
   
    y = sum(Mextracted) / C0solid / V * 100;  % Percentage of extract extracted
                                              % Total, from all stages
    
end    