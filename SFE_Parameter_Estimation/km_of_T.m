function km = km_of_T(T, theta)

% I will get back to you, promised!
% Dataset comes in as theta

    R  = theta{12};     % R = 8.31446;
    EA = theta{22};
    betah = theta{23};

    km = betah .* exp( -EA ./ R ./ T);                   % Estimated Di for some T
  
end