function Di = Di_of_T(T, theta)

% I will get back to you, promised!
% Dataset comes in as theta

    R   = theta{12};     % R = 8.31446;
    EA = theta{15};
    betah = theta{16};

    Di = betah .* exp( -EA ./ R ./ T);                   % Estimated Di for some T
  
end