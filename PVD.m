function [x, y, M_tilde, D] = PVD(N, B, files)
% Perform Population Value Decomposition prior to estimating Degrees of
% Mediation
%
%
% INPUT:
%
% N     - Number of subjects
% B     - Number of components (should be lesss than the min number of
%               trials per subject)
% path  - data directory
%
% OUTPUT:
%
% x      - treatment (N X 1 vector)
% y      - outcome (N X 1 vector)
% Mtilde - mediator after PVD (N X B matrix)
% D      - Weight projection matrix (B x voxels)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Read in input variables.
x = [];
y = [];
M = [];


V = [];
Res = []; 
for i=1:N

    % Load data - should be mat-file containing X (voxels-by-trial), Temp, and Rate
    load(files{i});
    
    x = [x; Temp];  % Trial-specific temperatures
    y = [y; Rate];  % Trial-specific pain ratings

    % Compute and store trial-specific SVDs
   
    [Vi,Si,Ui] =svd(X,0);
    Res{i}.U = Ui;
    Res{i}.S = Si;
    Res{i}.V = Vi;
    disp(i)
   
%    V = [V Vi(:,1:B)];

    V = [V Vi];
end


%%    

% Compute D matrix

VV = V'*V;

% Perfom SVD on VV 

[Ct, Bt, ~] = svd(VV);
Bt_inv= diag(1./sqrt(diag(Bt)));
CB_inv = Ct*Bt_inv;

% Obtain A_k tilde, i.e. A_k tilde = V_k tilde C_tilde B_tilde ^ {-1}

A_tilde = V*CB_inv;

D = A_tilde(:,1:B);



%%

% Compute M_tilde

M_tilde = [];

for i = 1:N,
    Ui = Res{i}.U;
    Si = Res{i}.S;
    Vi = Res{i}.V;
    
 %   tmp = Ui(:,1:B)*Si(1:B,1:B)*Vi(:,1:B)'*D;
    tmp = Ui*Si(:,1:B)*Vi(:,1:B)'*D;
    M_tilde = [M_tilde; tmp];
end

