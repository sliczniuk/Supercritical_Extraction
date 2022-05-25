function Di = Di_of_T(T, theta,k)

% I will get back to you, promised!
% Dataset comes in as theta

    R   = theta{12};     % R = 8.31446;
    %EA = theta{15};
    betah = theta{16};

    % k0 = [EA_Di, betah_Di, EA_km, betah_km];
    EA = k(1);
    %betah = k(2);

    Di = betah .* exp( -EA ./ R ./ T);                   % Estimated Di for some T
  
end