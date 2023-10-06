clear;

% Load data into Matlab
load matlab_data.txt

aftershocks = matlab_data;

% The variable aftershocks is a 458x4 matrix that contains 458 aftershock events 
% each with its time (in seconds from main earthquake) and position in latitude, 
% longitude, and depth respectively
% The last entry in the matrix is the main earthquake event (Magnitude 7.8)
% that caused the aftershocks.

%% Section 1 - Getting the dependance of aftershock frequency on days. p = 1

% Convert seconds to days
days = aftershocks(:,1)./(60*60*24);
% Flip, starting from the first event and finishing with the last
days = flip(days);

% Plot histogram to look at frequency (in days) of earthquakes
histogram(days); xlabel('Days'); ylabel('Number of aftershocks');
title('Frequency of Aftershocks from Feb. 06, 2023 to Apr. 19, 2023');

% Get number of aftershocks for each day. Total 73 days.
num_shocks = zeros(73,1);    % Frequency of aftershocks
for i = 1:length(days)
    j = floor(days(i));      % Get integer days
    num_shocks(j + 1) = num_shocks(j + 1) + 1;
end

% Plot number of aftershocks as a function of days
t = 1:1:73;
figure; plot(t, num_shocks, 'ko', 'markerfacecolor', 'k', 'markersize', 4);
grid("on"); xlabel('Days'); ylabel('Number of Aftershocks'); hold on;

% Do L2 fit (assuming no uncertainty)
G = [ones(73, 1) 1./transpose(t)];
m = linsolve(G, num_shocks);
num_predicted = m(1) + m(2)./t;     % L2 fit

% Plot predicted values of L2
plot(t, num_predicted, 'b*', 'markerfacecolor', 'b', 'markersize', 4);

% Do L1 fit:
misfitL1 = @(m)sum(abs(G*m - num_shocks));   % Objective function
m0 = [0.1; 163];   % Guess the model's parameters
% Get minimum of function
m1 = fminsearch(misfitL1, m0);

% Plot predicted values of L1
num_predicted1 = m1(1) + m1(2)./t;     % L1 fit
plot(t, num_predicted1, 'r*', 'markerfacecolor', 'r', 'markersize', 4);

% Do minmax Linf:
misfitLinf = @(m)max(abs(G*m - num_shocks));   % Objective function
m0 = [0.1; 163];   % Guess the model's parameters
% Get minimum of function
minf = fminsearch(misfitLinf, m0);

% Plot predicted values of Linf
num_predictedinf = minf(1) + minf(2)./t;     % L1 fit
plot(t, num_predictedinf, 'm*', 'markerfacecolor', 'm', 'markersize', 4);
legend('Measured Data', 'Predicted Data L2','Predicted Data L1', 'Predicted Data L infinity');

%% Section 1.1 p = 0.5

% Plot number of aftershocks as a function of days
t = 1:1:73;
figure; plot(t, num_shocks, 'ko', 'markerfacecolor', 'k', 'markersize', 4);
grid("on"); xlabel('Days'); ylabel('Number of Aftershocks'); hold on;

% Do inversion (assuming no uncertainty)
G = [ones(73, 1) 1./transpose(t).^0.5];
m = linsolve(G, num_shocks);
num_predicted = m(1) + m(2)./(t.^0.5);     % L2 fit

% Do L-1 fit assuming outliers don't matter as much:
misfitL1 = @(m)sum(abs(G*m - num_shocks));   % Objective function
m0 = [0.1; 163];   % Guess the model's parameters
% Get minimum of function
m1 = fminsearch(misfitL1, m0);

% Plot predicted values of L2
plot(t, num_predicted, 'b*', 'markerfacecolor', 'b', 'markersize', 4);

% Plot predicted values of L1
num_predicted1 = m1(1) + m1(2)./(t.^0.5);     % L1 fit
plot(t, num_predicted1, 'r*', 'markerfacecolor', 'r', 'markersize', 4);
legend('Measured Data', 'Predicted Data L2','Predicted Data L1');

% Do minmax
misfitLinf = @(m)max(abs(G*m - num_shocks));   % Objective function
m0 = [0.1; 163];   % Guess the model's parameters
% Get minimum of function
minf = fminsearch(misfitLinf, m0);

% Plot predicted values of Linf
num_predictedinf = minf(1) + minf(2)./(t.^0.5);     % L1 fit
plot(t, num_predictedinf, 'm*', 'markerfacecolor', 'm', 'markersize', 4);
legend('Measured Data', 'Predicted Data L2','Predicted Data L1', 'Predicted Data L infinity');
title('p = 0.5');

%% Section 1.2 p = 2

% Plot number of aftershocks as a function of days
t = 1:1:73;
figure; plot(t, num_shocks, 'ko', 'markerfacecolor', 'k', 'markersize', 4);
grid("on"); xlabel('Days'); ylabel('Number of Aftershocks'); hold on;

% Do inversion (assuming no uncertainty)
G = [ones(73, 1) 1./transpose(t).^2];
m = linsolve(G, num_shocks);
num_predicted = m(1) + m(2)./(t.^2);     % L2 fit

% Do L-1 fit assuming outliers don't matter as much:
misfitL1 = @(m)sum(abs(G*m - num_shocks));   % Objective function
m0 = [0.1; 163];   % Guess the model's parameters
% Get minimum of function
m1 = fminsearch(misfitL1, m0);

% Plot predicted values of L2
plot(t, num_predicted, 'b*', 'markerfacecolor', 'b', 'markersize', 4);

% Plot predicted values of L1
num_predicted1 = m1(1) + m1(2)./(t.^2);     % L1 fit
plot(t, num_predicted1, 'r*', 'markerfacecolor', 'r', 'markersize', 4);
legend('Measured Data', 'Predicted Data L2','Predicted Data L1');

% Do minmax
misfitLinf = @(m)max(abs(G*m - num_shocks));   % Objective function
m0 = [0.1; 163];   % Guess the model's parameters
% Get minimum of function
minf = fminsearch(misfitLinf, m0);

% Plot predicted values of Linf
num_predictedinf = minf(1) + minf(2)./(t.^2);     % L1 fit
plot(t, num_predictedinf, 'm*', 'markerfacecolor', 'm', 'markersize', 4);
legend('Measured Data', 'Predicted Data L2','Predicted Data L1', 'Predicted Data L infinity');
title('p = 2');

%% Section 2 - Getting the distribution and spread of aftershocks

% Get average of latitude
lat_avg = mean(aftershocks(:,2));   % Gives 37.6929 degrees

% Convert latitude to km (at the equator each degree is 111 km
lats = aftershocks(:,2).*111;       % Latitudes in km
lats0 = lats - mean(lats);          % Latitudes with center 0

% Convert longitude to km (at the average latitude)
long_dist = 6371*cosd(lat_avg)*pi/180;  % Distance in km for each degree
longs = aftershocks(:,3).*long_dist;    % Longitudes in km
longs0 = longs - mean(longs);           % Longitudes with center 0

% Remove average from depth:
depth0 = aftershocks(:,4) - mean(aftershocks(:,4));     % Depth with center 0

% Compute the covariance matrix:
A = [longs0 lats0 depth0];     % Matrix 39x3
covar = cov(A);

% Get eigenvalues and vectors
[eigen_vectors, eigen_values] = eigs(covar,3);

% Get standard deviations:
std1 = sqrt(covar(1,1));     % 95% of aftershocks are within 220 km from center longitudinally
std2 = sqrt(covar(2,2));     % 95% of aftershocks are within 152 km from center
std3 = sqrt(covar(3,3));     % 95% of aftershocks are within 7.4*2 = 14.8 km depth

% Plot
plot3(longs0, lats0, depth0, 'k.'); grid(); hold on;
xlabel('Longitude (km)'); ylabel('Latitude (km)'); 
zlabel('Depth (km)'); title('Distribution of Aftershocks');

% Plot direction of spread
v1 = eigen_vectors(:,1);        % Eigen vector 1 (largest spread)
v2 = eigen_vectors(:,2);        % Eigen vector 2
v3 = eigen_vectors(:,3);        % Eigen vector 3 (smallest spread)

% Magnify each direction with its standard deviation
S1 = 2*std1; E1 = -S1;         % Magnifying factors: S for start point and E for end point
S2 = 2*std2; E2 = -S2;
S3 = 2*std3; E3 = -S3;

% Plot
plot3([v1(1)*S1 v1(1)*E1], [v1(2)*S1 v1(2)*E1], ...
    [v1(3)*S1 v1(3)*E1], 'r-','linewidth', 2); hold on;
plot3([v2(1).*S2; v2(1).*E2], [v2(2).*S2; v2(2).*E2], ...
    [v2(3).*S2; v2(3).*E2], 'b-','linewidth', 2);
plot3([v3(1).*S3; v3(1).*E3], [v3(2).*S3; v3(2).*E3], ...
    [v3(3).*S3; v3(3).*E3], 'g-','linewidth', 2);

a1 = 1/long_dist; a2 = 1/111; a3 = 1;     % In terms of 1 km
daspect([a1 a2 a3]);

%% Section 2.1 Plot distribution with spread on map

xc = mean(aftershocks(:,3));     % Center points
yc = mean(aftershocks(:,2));

% Plot direction of spread
v1 = eigen_vectors(:,1);        % Eigen vector 1 (largest spread)
v2 = eigen_vectors(:,2);        % Eigen vector 2

% Magnify each direction with its standard deviation
S1 = 2*std1; E1 = -S1;         % Magnifying factors: S for start point and E for end point
S2 = 2*std2; E2 = -S2;

% Convert from km to degrees
v1 = v1./long_dist;
v2 = v2./111;

M1 = [v1(1)*S1 v1(1)*E1;
    v1(2)*S1 v1(2)*E1];
M2 = [v2(1).*S2 v2(1).*E2;
    v2(2).*S2 v2(2).*E2];

% Plot longitude and latitude with map
figure; plot(aftershocks(:,3), aftershocks(:,2), 'k.'); hold on;
xlabel('Longitude (km)'); ylabel('Latitude (km)'); 
title('Distribution of Aftershocks with Map');
load coastline.mat; plot(long, lat, 'k', 'linewidth', 2);            % Plot coastline on top

% Plot spread
plot( M1(1,:) + xc,  M1(2,:)+ yc, ...
     'r-','linewidth', 2); hold on;
plot(M2(1,:) + yc, M2(2,:) + yc, ...
    'b-','linewidth', 2);
legend('Aftershocks', 'Coastline', 'Spread 1', 'Spread 2');
