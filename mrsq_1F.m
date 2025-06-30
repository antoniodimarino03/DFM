function [mR2_omega, mR2_F_omega, freqs, t10_s_avg, t10_mR2_avg] = mrsq_1F(chi, D_X, series, m, freq_band)
% =========================================================================
% This function computes frequency-dependent variance decompositions for a
% SINGLE dynamic factor model estimated via the FHLR (2000) method. It is
% -------------------------------------------------------------------------
% INPUTS
% chi          = T x N common component matrix, an output from 'gdfm_twosided'.
% D_X          = N x (2h+1) matrix of ALL dynamic eigenvalues of the original
%                data, an output from 'gdfm_twosided'.
% series       = N x 1 cell array of series names.
% m            = (Optional) Covariogram truncation. Should match the 'm' used
%                in 'gdfm_twosided'. Default: calculated from T.
% freq_band    = (Optional) A 1x2 vector specifying the frequency band [w_low, w_high]
%                in radians (0 to pi) over which to average for the 'top 10'
%                series ranking. Default: business cycle band [pi/16, pi/3].
% OUTPUTS
% mR2_omega    = Marginal R-squared for each series (N x (h+1)). Represents
%                the proportion of the factor's power at a given frequency
%                attributable to each series.
% mR2_F_omega  = Marginal R-squared for the factor at each frequency (1 x (h+1)).
%                This is the fraction of total system power at each frequency
%                explained by the factor.
% freqs        = Vector of frequencies (1 x (h+1)) from 0 to pi at which
%                the analysis is performed.
% t10_s_avg    = Top 10 series that contribute most to the factor's power,
%                averaged over the specified frequency band.
% t10_mR2_avg  = Average marginal R-squared for the top 10 series over the band.
%
%% --- 1. Preliminary Settings and Input Validation ---
[T, N] = size(chi);
[N_eig, H] = size(D_X);

if N ~= N_eig
    error('Dimension mismatch: The number of series in chi (columns) must match the number of eigenvalues in D_X (rows).');
end

% The FHLR code computes for 2h+1 frequencies. We derive h.
if mod(H, 2) == 0
    error('The number of frequency points (columns of D_X) must be odd (2h+1).');
end
h = (H - 1) / 2;

% Set default for covariogram truncation 'm'
if nargin < 4 || isempty(m)
    m = floor(sqrt(T));
    fprintf('Using default covariogram truncation m = %d.\n', m);
end

if nargin < 5 || isempty(freq_band)
    % Default to business cycle frequencies for quarterly data (periods of 6 to 32 quarters)
    % freq = 2*pi/period, so period = 2*pi/freq.
    freq_band = [2*pi/32, 2*pi/6]; 
    fprintf('Using default business cycle frequency band for Top 10 ranking.\n');
end

%% --- 2. Generate Frequency Vector ---
freqs = (0:h) * (2*pi / H);
freq_indices = (h+1):H;
lambda_omega = D_X(:, freq_indices);

%% --- 3. Spectral Analysis of the Common Component 'chi' ---
% We need the power spectrum of the common component for each series. This is
% the diagonal of the spectral density matrix of 'chi'. We can get this by
% calling the 'spectral' helper function on the 'chi' matrix.
fprintf('Performing spectral analysis on the common component...\n');
[~, ~, Sigma_chi] = spectral(chi, N, m, h);

% Extract the diagonal elements (power of each series) for positive frequencies
chi_power_omega = zeros(N, h + 1);
for i = 1:(h+1)
    % The diagonal contains the power spectrum for each series.
    % Take real part to discard negligible imaginary parts due to numerical error.
    chi_power_omega(:, i) = real(diag(Sigma_chi(:, :, freq_indices(i))));
end

%% --- 4. Calculate Frequency-Domain R-squared Metrics ---

% --- Marginal R-squared for the Factor (mR2_F_omega) ---
% This is the fraction of total system variance at each frequency explained
% by the first dynamic factor.
total_power_omega = sum(lambda_omega, 1);
factor_power_omega = lambda_omega(1, :); % Power of the first factor is the first eigenvalue

mR2_F_omega = factor_power_omega ./ total_power_omega;
mR2_F_omega(total_power_omega == 0) = 0; % Handle division by zero

% --- Marginal R-squared for each Series (mR2_omega) ---
% This is each series' contribution to the *factor's* total power at each frequency.
total_chi_power_omega = sum(chi_power_omega, 1);

mR2_omega = bsxfun(@rdivide, chi_power_omega, total_chi_power_omega);
mR2_omega(isnan(mR2_omega)) = 0; % Handle cases where total factor power is zero

%% --- 5. Identify Top 10 Series based on Average Contribution ---
% We average each series' contribution over the specified frequency band
band_indices = find(freqs >= freq_band(1) & freqs <= freq_band(2));
if isempty(band_indices)
    warning('The specified frequency band contains no frequency points. Top 10 will not be computed.');
    t10_s_avg = cell(10, 1);
    t10_mR2_avg = NaN(10, 1);
    return;
end

mR2_avg_in_band = mean(mR2_omega(:, band_indices), 2);
[sorted_mR2_avg, ind] = sort(mR2_avg_in_band, 'descend');
num_top = min(10, N);
t10_s_avg = cell(num_top, 1);
t10_mR2_avg = NaN(num_top, 1);
for i = 1:num_top
    t10_s_avg{i} = series{ind(i)};
    t10_mR2_avg(i) = sorted_mR2_avg(i);
end

fprintf('Analysis complete.\n');

end

%===================================================================================
function [P_chi, D_chi, Sigma_chi] = spectral(X, q, m, h)
% Computes spectral decomposition for the data matrix in input.
[T,n] = size(X);
if nargin < 2, disp('ERROR MESSAGE: Too few input arguments'); return; end
if nargin == 2, m = floor(sqrt(T)); h = m; end
if nargin == 3, h = m; end

M = 2*m+1;
B = 1 - abs(-m:m)/m;
Gamma_k = zeros(n,n,M);
for k = 1:m+1
    Gamma_k(:,:,m+k) = B(m+k)*(X(k:T,:))'*(X(1:T+1-k,:))/(T-k+1);
    Gamma_k(:,:,m-k+2) = Gamma_k(:,:,m+k)';
end

H = 2*h+1;
Factor = exp(-1i*(-m:m)'*(-2*pi*h/H:2*pi/H:2*pi*h/H));
Sigma_X = zeros(n,n,H);
for j = 1:n
    Sigma_X(j,:,:) = squeeze(Gamma_k(j,:,:))*Factor;
end

P_chi = zeros(n,q,H);
D_chi = zeros(q,H);
Sigma_chi = zeros(n,n,H);

opt.disp = 0;
opt.isreal = false; 
opt.issym = true;   


for j = 1:H
    [P, D] = eig(squeeze(Sigma_X(:,:,j)), 'vector');
    [sorted_D, IX] = sort(real(D), 'descend');
    P_chi(:,:,j) = P(:, IX(1:q));
    D_chi(:,j) = sorted_D(1:q);
    Sigma_chi(:,:,j) = P_chi(:,:,j) * diag(D_chi(:,j)) * P_chi(:,:,j)';
end
end