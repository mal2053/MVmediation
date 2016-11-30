function [t,pv,Wboot] = BootDM(x, y, M_tilde, W, Theta, Dt, Bsamp)
% Perform bootstrap on Degrees of Mediation (DMs)
%
%
% INPUT:
%
% x      - treatment (N X 1 vector)
% y      - outcome (N X 1 vector)
% Mtilde - mediator after PVD (N X B matrix)
% Dt     - inverse weight projection matrix (B x voxels)
% W      - weights from previous directions
% Theta  - parametersfrom previous directions
% Bsamp  - number of bootstrap samples
% numdir - number of DMs to perform bootstarp
%
% OUTPUT:
%
% Wboot  - Bootstrapped weights (Voxels X Bsamp X numdir)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bsamp = 200;
% numdir = 1;
p = size(Dt,2);
len = length(y);
Wboot = zeros(p,Bsamp);

for i=1:Bsamp,
    
    disp(i)
    
    % Bootstrap samples of the data
    ind = ceil(unifrnd(0,len,len,1));
    xB = x(ind);
    yB = y(ind);
    MB = M_tilde(ind,:);
    
    
    % Estimate the nth DM
    [w_n, ~, ~]= DirectionsMediationN(xB,yB,MB,W,Theta);
    Wboot(:,i) = pinv(Dt)*w_n;
          
end

% Compute p-values
t = zeros(p,1);
pv = zeros(p,1);

opt = statset('MaxIter',1000);

for i=1:p,
    try
        res = fitgmdist(Wboot(i,:)',2, 'Options', opt) ;
        t(i) = min(abs(res.mu(1)/sqrt(res.Sigma(1))) ,  abs(res.mu(2)/sqrt(res.Sigma(2))));     
    catch
        t(i) = 0;
    end
    pv(i) = 2*(1 - tcdf(t(i),len-1));
end


