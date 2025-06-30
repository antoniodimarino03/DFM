% =========================================================================
% PART 1: LOAD AND LABEL DATA

% Load data from CSV file
dum=importdata("import the data");

series=dum.textdata(1,2:end);

% Transformation numbers
tcode=dum.data(1,:);

% Raw data
rawdata=dum.data(2:end,:);

% Quarter/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_quarter=ceil(final_datevec(2)/3);  % Convert month to quarter
final_year=final_datevec(1);

% Dates (quarterly) are of the form YEAR+QUARTER/4
% e.g. Q1 1970 is represented as 1970+1/4, Q3 1970 is 1970+3/4
% Dates go from 1959:Q1 to final_year:final_quarter (see above)
dates = (1959+1/4:1/4:final_year+final_quarter/4)';

% T = number of quarters in sample
T=size(dates,1);
rawdata=rawdata(1:T,:);

% =========================================================================
% PART 2: PROCESS DATA
% Transform raw data to be stationary using auxiliary function
% prepare_missing()
yt=prepare_missing(rawdata,tcode);
% Reduce sample to usable dates: remove first two quarters because some
% series have been first differenced
yt=yt(3:T,:);
dates=dates(3:T,:);

