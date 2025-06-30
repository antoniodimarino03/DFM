function [X_forecast] = gdfm17_forecast(CL, v, h_max, sigma, mu)
%The function performs forecast as in Forni, Hallin, Lippi, Zaffaroni 2017
%forecasting equation:
% x̂ᵢ,t+h = Σⱼ₌₁q Σₖ₌₀∞ ĉᵢⱼ,h+k ûⱼ,t₋k

% INPUT:
%   CL          : n x q x nlagsimp matrix of impulse responses
%   v           : T x q matrix of common shocks/factors 
%   h_max       : maximum forecast horizon
% OUTPUT:
%   X_forecast  : n x h_max matrix of forecasts for each variable and horizon

[n, q, nlagsimp] = size(CL);
[T, ~] = size(v);

% Initialize forecast matrix
X_forecast = zeros(n, h_max);

% Direct forecasting using impulse responses
for h = 1:h_max
    for i = 1:n
        forecast_sum = 0;
        
        % Sum over factors j=1 to q
        for j = 1:q
            factor_contribution = 0;
            
            % Sum over lags k=0 to min(T-1, nlagsimp-h)
            max_lag = min(T-1, nlagsimp-h);
            
            for k = 0:max_lag
                if (h + k) <= nlagsimp && (T - k) > 0
                    % c_ij,h+k * u_j,T-k
                    factor_contribution = factor_contribution + ...
                        CL(i, j, h + k) * v(T - k, j);
                end
            end
            
            forecast_sum = forecast_sum + factor_contribution;
        end
        
        X_forecast(i, h) = forecast_sum;
    end
end

X_forecast = X_forecast';
X_forecast = X_forecast.*sigma+mu;
end