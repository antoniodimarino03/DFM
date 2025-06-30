function [ehat,Fhat,lamhat,ve2,x2,icstar] = factors_SW(x,kmax,jj,DEMEAN)
% =========================================================================
% DESCRIPTION
% This program estimates a set of factors for a given dataset using
% principal component analysis. The number of factors estimated is
% determined by an information criterion specified by the user. 
% This version assumes NO missing values in the dataset.
%
% -------------------------------------------------------------------------
% INPUTS
%           x       = dataset (one series per column) - NO MISSING VALUES
%           kmax    = an integer indicating the maximum number of factors
%                     to be estimated; if set to 99, the number of factors
%                     selected is forced to equal 8
%           jj      = an integer indicating the information criterion used 
%                     for selecting the number of factors; it can take on 
%                     the following values:
%                           1 (information criterion PC_p1)
%                           2 (information criterion PC_p2)
%                           3 (information criterion PC_p3)      
%           DEMEAN  = an integer indicating the type of transformation
%                     performed on each series in x before the factors are
%                     estimated; it can take on the following values:
%                           0 (no transformation)
%                           1 (demean only)
%                           2 (demean and standardize)
%                           3 (recursively demean and then standardize) 
%
% OUTPUTS
%           ehat    = difference between x and values of x predicted by
%                     the factors
%           Fhat    = set of factors
%           lamhat  = factor loadings
%           ve2     = eigenvalues of x2'*x2 (where x2 is the transformed dataset)
%           x2      = transformed dataset (same as input since no missing values)
%
% -------------------------------------------------------------------------
% SUBFUNCTIONS
%
% baing() - selects number of factors
%
% pc2() - runs principal component analysis
%
% minindc() - finds the index of the minimum value for each column of a
%       given matrix
%
% transform_data() - performs data transformation
%
% =========================================================================
% PART 1: CHECKS

% Check that kmax is an integer between 1 and the number of columns of x, or 99
if ~((kmax<=size(x,2) && kmax>=1 && floor(kmax)==kmax) || kmax==99)
    error('Input kmax is specified incorrectly.');
end

% Check that jj is one of 1, 2, 3
if jj~=1 && jj~=2 && jj~=3
    error('Input jj is specified incorrectly.');
end

% Check that DEMEAN is one of 0, 1, 2, 3
if DEMEAN ~= 0 && DEMEAN ~= 1 && DEMEAN ~= 2 && DEMEAN ~= 3
    error('Input DEMEAN is specified incorrectly.');
end

% Check for missing values and warn user
if any(isnan(x(:)))
    warning('Dataset contains missing values. Consider using factors_em() instead.');
end

% =========================================================================
% PART 2: DATA TRANSFORMATION AND FACTOR ESTIMATION

% Transform the data according to DEMEAN specification
[x2,mut,sdt] = transform_data(x,DEMEAN);

% Determine the number of factors to estimate
if kmax ~= 99
    [icstar,~,~,~] = baing(x2,kmax,jj);
else
    icstar = 8;
end

% Display selected number of factors
fprintf('Number of factors selected: %d\n', icstar);

% Run principal components analysis
[chat,Fhat,lamhat,ve2] = pc2(x2,icstar);

% Calculate residuals (difference between original and predicted values)
ehat = x - chat.*sdt - mut;

% Note: x2 is returned as the transformed dataset (not the original x since no missing values)
% If you want the original untransformed data, use x instead of x2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ic1,chat,Fhat,eigval]=baing(X,kmax,jj)
% =========================================================================
% DESCRIPTION
% This function determines the number of factors to be selected for a given
% dataset using one of three information criteria specified by the user.
% The user also specifies the maximum number of factors to be selected.
%
% -------------------------------------------------------------------------
% INPUTS
%           X       = dataset (one series per column)
%           kmax    = an integer indicating the maximum number of factors
%                     to be estimated
%           jj      = an integer indicating the information criterion used 
%                     for selecting the number of factors; it can take on 
%                     the following values:
%                           1 (information criterion PC_p1)
%                           2 (information criterion PC_p2)
%                           3 (information criterion PC_p3)    
%
% OUTPUTS
%           ic1     = number of factors selected
%           chat    = values of X predicted by the factors
%           Fhat    = factors
%           eigval  = eigenvalues of X'*X (or X*X' if N>T)
%
% -------------------------------------------------------------------------
% SUBFUNCTIONS USED
%
% minindc() - finds the index of the minimum value for each column of a
%       given matrix
%
% =========================================================================

% Number of observations per series (i.e. number of rows)
T=size(X,1);

% Number of series (i.e. number of columns)
N=size(X,2);

% Total number of observations
NT=N*T;

% Number of rows + columns
NT1=N+T;

% =========================================================================
% PART 2: OVERFITTING PENALTY

% Allocate memory for overfitting penalty
CT=zeros(1,kmax);

% Array containing possible number of factors that can be selected (1 to kmax)
ii=1:1:kmax;

% The smaller of N and T
GCT=min([N;T]);

% Calculate penalty based on criterion determined by jj. 
switch jj
    case 1  % Criterion PC_p1
        CT(1,:)=log(NT/NT1)*ii*NT1/NT;
    case 2  % Criterion PC_p2
        CT(1,:)=(NT1/NT)*log(min([N;T]))*ii;
    case 3  % Criterion PC_p3
        CT(1,:)=ii*log(GCT)/GCT;
end

% =========================================================================
% PART 3: SELECT NUMBER OF FACTORS

% Run principal component analysis
if T<N 
    % Singular value decomposition
    [ev,eigval,~]=svd(X*X'); 
    % Components
    Fhat0=sqrt(T)*ev;
    % Loadings
    Lambda0=X'*Fhat0/T;
else
    % Singular value decomposition
    [ev,eigval,~]=svd(X'*X);
    % Loadings
    Lambda0=sqrt(N)*ev;
    % Components
    Fhat0=X*Lambda0/N;
end

% Preallocate memory
Sigma=zeros(1,kmax+1); % sum of squared residuals divided by NT
IC1=zeros(size(CT,1),kmax+1); % information criterion value

% Loop through all possibilities for the number of factors 
for i=kmax:-1:1
    % Identify factors as first i components
    Fhat=Fhat0(:,1:i);
    % Identify factor loadings as first i loadings
    lambda=Lambda0(:,1:i);
    % Predict X using i factors
    chat=Fhat*lambda';
    % Residuals from predicting X using the factors
    ehat=X-chat;
    % Sum of squared residuals divided by NT
    Sigma(i)=mean(sum(ehat.*ehat/T));
    % Value of the information criterion when using i factors
    IC1(:,i)=log(Sigma(i))+CT(:,i);
end

% Sum of squared residuals when using no factors
Sigma(kmax+1)=mean(sum(X.*X/T));
% Value of the information criterion when using no factors
IC1(:,kmax+1)=log(Sigma(kmax+1));

% Number of factors that minimizes the information criterion
ic1=minindc(IC1')';

% Set ic1=0 if ic1>kmax (i.e. no factors selected)
ic1=ic1 .*(ic1 <= kmax);

% Save other output
Fhat=Fhat0(:,1:kmax); % factors
Lambda=Lambda0(:,1:kmax); % factor loadings
chat=Fhat*Lambda'; % predict X using kmax factors
eigval=diag(eigval); % eigenvalues


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [chat,fhat,lambda,ss]=pc2(X,nfac)
% =========================================================================
% DESCRIPTION
% This function runs principal component analysis.
%
% -------------------------------------------------------------------------
% INPUTS
%           X      = dataset (one series per column)
%           nfac   = number of factors to be selected
%
% OUTPUTS
%           chat   = values of X predicted by the factors
%           fhat   = factors scaled by (1/sqrt(N)) where N is the number of series
%           lambda = factor loadings scaled by number of series
%           ss     = eigenvalues of X'*X 
%
% =========================================================================

% Number of series in X (i.e. number of columns)
N=size(X,2);

% Singular value decomposition: X'*X = U*S*V'
[U,S,~]=svd(X'*X);

% Factor loadings scaled by sqrt(N)
lambda=U(:,1:nfac)*sqrt(N);

% Factors scaled by 1/sqrt(N)
fhat=X*lambda/N;

% Estimate initial dataset X using the factors
chat=fhat*lambda';

% Identify eigenvalues of X'*X
ss=diag(S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pos=minindc(x)
% =========================================================================
% DESCRIPTION
% This function finds the index of the minimum value for each column of a
% given matrix.
%
% -------------------------------------------------------------------------
% INPUT
%           x   = matrix 
%
% OUTPUT
%           pos = column vector with pos(i) containing the row number
%                 corresponding to the minimum value of x(:,i)
%
% =========================================================================

% Number of rows and columns of x
nrows=size(x,1);
ncols=size(x,2);

% Preallocate memory for output array
pos=zeros(ncols,1);

% Create column vector 1:nrows
seq=(1:nrows)';

% Find the index of the minimum value of each column in x
for i=1:ncols
    % Minimum value of column i
    min_i=min(x(:,i));
    % Column vector containing row number of minimum value
    colmin_i= seq.*((x(:,i)-min_i)==0);
    % Produce an error if the minimum value occurs more than once
    if sum(colmin_i>0)>1
        error('Minimum value occurs more than once.');
    end
    % Obtain the index of the minimum value
    pos(i)=sum(colmin_i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x22,mut,sdt]=transform_data(x2,DEMEAN)
% =========================================================================
% DESCRIPTION
% This function transforms a given set of series based upon the input
% variable DEMEAN.
%
% -------------------------------------------------------------------------
% INPUTS
%           x2      = set of series to be transformed (one series per column)
%           DEMEAN  = transformation type:
%                           0 (no transformation)
%                           1 (demean only)
%                           2 (demean and standardize)
%                           3 (recursively demean and then standardize) 
%
% OUTPUTS
%           x22     = transformed dataset
%           mut     = matrix containing the values subtracted from x2
%           sdt     = matrix containing the values that x2 was divided by
%
% =========================================================================

% Number of observations and series
T=size(x2,1);
N=size(x2,2);

% Perform transformation based on DEMEAN
switch DEMEAN
    case 0  % No transformation
        mut=repmat(zeros(1,N),T,1);
        sdt=repmat(ones(1,N),T,1);
        x22=x2;
        
    case 1  % Demean only
        mut=repmat(mean(x2),T,1);
        sdt=repmat(ones(1,N),T,1);
        x22=x2-mut;
        
    case 2  % Demean and standardize 
        mut=repmat(mean(x2),T,1);
        sdt=repmat(std(x2),T,1);
        x22=(x2-mut)./sdt;
        
    case 3  % Recursively demean and then standardize
        mut=NaN(size(x2));
        for t=1:T
            mut(t,:)=mean(x2(1:t,:),1);
        end
        sdt=repmat(std(x2),T,1);
        x22=(x2-mut)./sdt; 
end