function forecasts = forecast_SW(X, Fhat, h, s, p)
% =========================================================================
% DESCRIPTION
% This function performs multi-step direct forecasting using factor-augmented
% models as in Stock and Watson (2002). For each variable and each horizon,
% it estimates the forecasting equation:
% y_{t+h} = α_h(L^s)F_t + β_h(L)y_t + ε_{t+h}
%
% where:
% α_h(L^s)F_t = α_{0h}F_t + α_{1h}F_{t-1} + ... + α_{sh}F_{t-s}
% β_h(L)y_t = β_{1h}y_{t-1} + β_{2h}y_{t-2} + ... + β_{ph}y_{t-p}
%
% -------------------------------------------------------------------------
% INPUTS
%   X       = T x N matrix of original data (T observations, N variables)
%   Fhat    = T x q matrix of estimated factors (from factors_SW function)
%   h       = maximum forecast horizon (integer)
%   s       = number of lags of factors to include in forecasting equation
%   p       = number of lags of dependent variable to include
%
% OUTPUTS
%   forecasts = h x N matrix of forecasts where forecasts(i,j) is the
%              forecast of variable j at horizon i
%
% -------------------------------------------------------------------------
% NOTES
% - Uses direct method (not iterative)
% - Forecast origin is the last observation in the sample
% - Requires sufficient observations for lags: T > max(s, p) + h
% =========================================================================

[T, N] = size(X);
[T_factors, q] = size(Fhat);

if T ~= T_factors
    error('X and Fhat must have the same number of time periods.');
end

if T <= max(s, p) + h
    error('Insufficient observations. Need T > max(s, p) + h = %d', max(s, p) + h);
end

if h < 1 || s < 0 || p < 0
    error('Invalid lag specifications. h >= 1, s >= 0, p >= 0 required.');
end

% Lagged factors and lagged dependent variables

max_lag = max(s, p);
effective_start = max_lag + 1;
effective_end = T - h;
effective_T = effective_end - effective_start + 1;

if effective_T <= 0
    error('No observations available for estimation after accounting for lags and forecast horizon.');
end

fprintf('Using %d factors for forecasting.\n', q);
fprintf('Effective sample for estimation: %d observations (from %d to %d)\n', ...
        effective_T, effective_start, effective_end);

forecasts = NaN(h, N);

fprintf('Performing direct forecasting...\n');

for var_idx = 1:N
    if mod(var_idx, 10) == 0
        fprintf('Processing variable %d of %d...\n', var_idx, N);
    end
    
    y = X(:, var_idx);  % Current variable
    
    for horizon = 1:h
 
        % C (y_{t+h})

        y_future = y((effective_start + horizon):(effective_end + horizon));
        n_regs = (s + 1) * q + p;  % Number of regressors
        X_reg = NaN(effective_T, n_regs);
        reg_idx = 1;
        
        % Add factor lags: F_t, F_{t-1}, ..., F_{t-s}
        for lag = 0:s
            factor_lag = Fhat((effective_start - lag):(effective_end - lag), :);
            X_reg(:, reg_idx:(reg_idx + q - 1)) = factor_lag;
            reg_idx = reg_idx + q;
        end
        
        % Add dependent variable lags: y_{t-1}, y_{t-2}, ..., y_{t-p}
        for lag = 1:p
            y_lag = y((effective_start - lag):(effective_end - lag));
            X_reg(:, reg_idx) = y_lag;
            reg_idx = reg_idx + 1;
        end
        
        % Handle case where there are no regressors
        if n_regs == 0
            forecasts(horizon, var_idx) = mean(y_future);
        else
            % Add constant term
            X_reg_const = [ones(effective_T, 1), X_reg];
            
            % OLS estimation: beta = (X'X)^(-1)X'y
            try
                beta = (X_reg_const' * X_reg_const) \ (X_reg_const' * y_future);

                % Create regressor vector for forecast (using data up to time T)
                X_forecast = NaN(1, n_regs);
                reg_idx = 1;
                
                % Factor values: F_T, F_{T-1}, ..., F_{T-s}
                for lag = 0:s
                    if T - lag >= 1
                        X_forecast(reg_idx:(reg_idx + q - 1)) = Fhat(T - lag, :);
                    else
                        X_forecast(reg_idx:(reg_idx + q - 1)) = 0;  % Pad with zeros if not enough data
                    end
                    reg_idx = reg_idx + q;
                end
                
                % Dependent variable lags: y_{T-1}, y_{T-2}, ..., y_{T-p}
                for lag = 1:p
                    if T - lag >= 1
                        X_forecast(reg_idx) = y(T - lag);
                    else
                        X_forecast(reg_idx) = 0;  % Pad with zeros if not enough data
                    end
                    reg_idx = reg_idx + 1;
                end
                
                % Generate forecast
                X_forecast_const = [1, X_forecast];
                forecasts(horizon, var_idx) = X_forecast_const * beta;
                
            catch ME
                warning('OLS estimation failed for variable %d, horizon %d: %s', ...
                        var_idx, horizon, ME.message);
                forecasts(horizon, var_idx) = y(T);  % Use last observation as forecast
            end
        end
    end
end
end