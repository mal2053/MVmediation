function [w_N, theta_N, lambda]= DirectionsMediationN(x,y,m, W, Theta, varargin)
% Compute the Nth Direction of Mediation
%
% This code can be used iteratively to compute each direction of mediation
%
% INPUT:
%
% x     - treatment (N X 1 vector)
% m     - mediator (N X p matrix)
% y     - outcome (N X 1 vector)
% W     - weights of previous directions (cell structure). Use [] if first direction.
% Theta - parameters of previous directions (cell structure). Use [] if first direction.
%
% OUTPUT:
%
% w_N     - weights for the Nth direction of mediation
% theta_N - parameters for Nth direction of mediation
% lambda  - lambda value
%
% Example:
%
% First 3 directions:
%
% [w_1, theta_1, lambda]= DirectionsMediationN(x1,y1,m1, [],[]);
% W{1} = w_1;
% Theta{1} = theta_1;
% [w_2, theta_2, lambda]= DirectionsMediationN(x1,y1,m1, W, Theta);
% W{2} = w_2;
% Theta{2} = theta_2;
% [w_3, theta_3, lambda]= DirectionsMediationN(x1,y1,m1, W, Theta);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set-up
imax = 1000;             % Number of iterations
tol = 10^-6;

len = size(m,2);
N1 = length(W);          % Number of directions previously computed.

I = eye(len);
lambda = 0;

if (N1 >0)

    theta = Theta{N1};
    tn = length(theta);
    theta_N = [theta(1:tn-1) 1 theta(tn)];   
    tn = tn+1;

    H = I;
    WM = zeros(len,N1);      % Create matrix of weight values
    for i=1:N1,
        ww = W{i};
        H = H - ww*ww';
        WM(:,i) = ww;
    end
    WM(:,(N1+1)) = ones(len,1);
else
    
    WM = ones(len,1);
    tn =5;
    theta_N = ones(1,tn);
    H = I;

end


w_N = WM(:,end);

count=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start iterative procedure

for i=1:imax,
    
    MW = m*WM;
    mw_N = MW(:,N1+1);
    J = ones(size(m,1),1);
  
    X1 = [J x MW];
    X2 = [J x];
    
    
    % Compute sample variances
    s1 = sqrt(sum((y - X1*((X1'*X1)\X1')*y).^2)/(size(m,1)- N1 - 3));
    s2 = sqrt(sum((mw_N - X2*((X2'*X2)\X2')*mw_N).^2)/(size(m,1)- N1 - 2));
    
    % Compute psi and phi
    HmmH = H*(m'*m)*H;
    psi = (1/s1^2)* HmmH*theta_N(tn-1)^2 + (1/s2^2)*HmmH;
    if (N1 >0)
        phi = (theta_N(tn-1)/s1^2)* H* m'*(y - theta_N(3) - x*theta_N(tn) -MW(:,1:N1)*theta_N(4:(tn-2))') + (1/s2^2)*H*m'*(theta_N(1)+x*theta_N(2));
    else
        phi = (theta_N(tn-1)/s1^2)* H* m'*(y - theta_N(3) - x*theta_N(tn)) + (1/s2^2)*H*m'*(theta_N(1)+x*theta_N(2));
    end
    
    % Estimate lambda_i
    lambda0 = lambda;
    fun = @(lambda)((pinv(lambda*H'*H+ psi) * phi)'* (H'*H) * pinv(lambda*H'*H+ psi)*phi - 1);
    try
        lambda = fzero(fun,0); 
    catch
        lambda = -100000;
    end
    
    % Compute W_i (new weights)
    w_N_new = (pinv(lambda*(H'*H) + psi)*phi);
    w_N_new = w_N_new./sqrt(sum(w_N_new.^2));


    % Estimate theta_i
    fun = @(theta_N)(sum((y-theta_N(3)-x*theta_N(tn)-MW*theta_N(4:(tn-1))').^2) / (1/s1^2));    
    theta_N_1 = fminsearch(fun,theta_N);
    
    fun = @(theta_N)(sum((mw_N-theta_N(1)-x*theta_N(2)).^2) / (1/s2^2));
    theta_N_2 = fminsearch(fun,theta_N);
 
    theta_N_new=[theta_N_2(1:2),theta_N_1(3:tn)];


    % Compute likelihood  
    L = - sum((y-theta_N_new(3)-x*theta_N_new(tn)-MW*theta_N_new(4:(tn-1))').^2) / (2*s1^2) - sum((mw_N-theta_N_new(1)-x*theta_N_new(2)).^2) / (2*s2^2);


    % Print results
    if (mod(i,100) == 1)
        fprintf('Calculating PDM number: %d \n',N1+1);
        fprintf('log-likelihood = %d \n',L);
        fprintf('lambda = %d \n',lambda);
        fprintf('norm = %d \n',w_N_new'*w_N_new);
    end  
    
    if (abs(w_N_new'*w_N_new - 1) > 1000*eps), 
        error('Something wrong')
    end
        

    % Determine whether to exit iterative procedure
    if (abs( sum((theta_N_new-theta_N).^2) )<tol && abs( sum((w_N_new - w_N).^2) )<tol)
        theta_N = theta_N_new;
        w_N = w_N_new; 
        WM(:,end)= w_N;
        break;
    else
        theta_N = theta_N_new;
        w_N = w_N_new; 
        WM(:,end)= w_N;
        count = count + 1;
    end
   
   

end
