function yt = prepare_missing(rawdata,tcode)
% =========================================================================
% DESCRIPTION:
% This function transforms raw data based on each series' transformation
% code. Optimized for datasets without missing values.
%
% -------------------------------------------------------------------------
% INPUT:
% rawdata = raw data (no missing values)
% tcode = transformation codes for each series
%
% OUTPUT:
% yt = transformed data
%
% -------------------------------------------------------------------------
% SUBFUNCTION:
% transxf: transforms a single series as specified by a
% given transfromation code
%
% =========================================================================
% APPLY TRANSFORMATION:
% Get dimensions
[T, N] = size(rawdata);
% Initialize output variable with proper size
yt = zeros(T, N);
% Perform transformation using subfunction transxf
for i = 1:N
    yt(:,i) = transxf(rawdata(:,i), tcode(i));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
function y = transxf(x, tcode)
% =========================================================================
% DESCRIPTION:
% This function transforms a single series (in a column vector) as specified
% by a given transformation code. Optimized for data without missing values.
%
% -------------------------------------------------------------------------
% INPUT:
% x = series (in a column vector) to be transformed (no missing values)
% tcode = transformation code (1-7)
%
% OUTPUT:
% y = transformed series (as a column vector)
%
% =========================================================================
% SETUP:
% Number of observations
n = length(x);
% Value close to zero for log transformations
small = 1e-6;

% =========================================================================
% TRANSFORMATION:
% Determine case 1-7 by transformation code
switch tcode
    case 1 % Level (i.e. no transformation): x(t)
        y = x;
        
    case 2 % First difference: x(t)-x(t-1)
        y = [NaN; diff(x)];
        
    case 3 % Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
        y = [NaN; NaN; diff(x, 2)];
        
    case 4 % Natural log: ln(x)
        if min(x) <= small
            error('Cannot take log of non-positive values in series');
        end
        y = log(x);
        
    case 5 % First difference of natural log: ln(x)-ln(x-1)
        if min(x) <= small
            error('Cannot take log of non-positive values in series');
        end
        log_x = log(x);
        y = [NaN; diff(log_x)];
        
    case 6 % Second difference of natural log
        if min(x) <= small
            error('Cannot take log of non-positive values in series');
        end
        log_x = log(x);
        y = [NaN; NaN; diff(log_x, 2)];
        
    case 7 % First difference of percent change
        pct_change = [NaN; diff(x) ./ x(1:end-1)];
        y = [NaN; NaN; diff(pct_change(2:end))];
        
    otherwise
        error('Invalid transformation code. Must be 1-7.');
end
end