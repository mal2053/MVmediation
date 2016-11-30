%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters
 
N = 26;     % Number of subjects
B = 35;     % Number of components in the original PVD 

path = '/Users/martinlindquist/Dropbox/nsf/';

W = [];             % Reduced Weights
Wfull = [];         % Full Weights
Theta = [];         % Theta parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Population Value Decomposition
 
% Files should be mat-file containing X (voxels-by-trial), Temp, and Rate 
files = filenames('/Users/martinlindquist/Dropbox/nsf/*.mat');
[x, y, M_tilde, D] = PVD(N, B, files); 
Dt = D';    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute first direction of mediation

[w_1, theta_1, ~]= DirectionsMediationN(x,y,M_tilde,[],[]);

W{1} = w_1;
Theta{1} = theta_1;
Wfull{1} = pinv(Dt)*w_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute second direction of mediation

[w_2, theta_2, ~]= DirectionsMediationN(x,y,M_tilde,W,Theta);

W{2} = w_2;
Theta{2} = theta_2;
Wfull{2} = pinv(Dt)*w_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute third direction of mediation

[w_3, theta_3, ~]= DirectionsMediationN(x,y,M_tilde,W,Theta);

W{3} = w_3;
Theta{3} = theta_3;
Wfull{3} = pinv(Dt)*w_3;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap Dms

Bsamp = 200;

% Bootstrap p-values for first DM
[t1,p1,Wboot1] = BootDM(x, y, M_tilde, [], [], Dt, Bsamp);

% Bootstrap p-values for first DM
[t2,p2,Wboot2] = BootDM(x, y, M_tilde, W{1}, Theta{1}, Dt, Bsamp);

% Bootstrap p-values for first DM
[t3,p3,Wboot3] = BootDM(x, y, M_tilde, W{1:2}, Theta{1:2}, Dt, Bsamp);
