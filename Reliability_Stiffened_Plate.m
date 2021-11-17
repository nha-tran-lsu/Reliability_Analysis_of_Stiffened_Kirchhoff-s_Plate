clc;
clear all;
clf;

%% Random variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
n_random_var = 2;                     % number of random variables
x            = zeros(n_random_var,1); % random variables (cirticle displacment, 
                                      % stiffness [N/mm^2], force [N])
mu           = zeros(n_random_var,1); % expected values
sigma        = zeros(n_random_var,1); % standard deviation

%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Criticle displacement
u0 = 14;
Emodule=1092000;      % Young elastic modulus  
Load=1;               % Uniform load

% Stiffness
mu(1)    = 1092000;
sigma(1) = 109200;

% Forces
mu(2)    = 1;
sigma(2) =  0.1;

%% Finding design point using FORM Method 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Variable declaration
eps   = 0.001; % stopping creterium for the fix point iteration
alpha = zeros(n_random_var,1);
dg    = zeros(n_random_var,1);
dx    = 0.001*sigma;
dbeta = 100;

% step 2
x = mu;

% step 3
y = (x-mu)./sigma;
beta = 9999;

step=0;
while (abs(dbeta) > 0.001)
    step=step+1     
    % step 4 and 5
   
    tmp_u1 = cal_SP(x(1), x(2));
    tmp_u2 = cal_SP(x(1)+dx(1), x(2));

    dg(1)  = (tmp_u2 - tmp_u1)/dx(1) * sigma(1);
    tmp_u2 = cal_SP(x(1), x(2)+dx(2));
    dg(2)  = (tmp_u2 - tmp_u1)/dx(2) * sigma(2);
   
    % step 6
    g = -u0 + tmp_u1;
    
    y = 1./(dg'*dg)*(dg'*y - g)*dg;
    
    % step 7
    tmp = beta;
    beta = sqrt(y'*y)
    dbeta = beta - tmp;
    
    % step 8
    x = mu + sigma.*y;
    
end

pf = 1 - normcdf(beta);
sprintf('Probability of failure: %g\n', pf)
pf;
% % % % % % % % % % % % % % % % % % % % 
% u0=10
% step =     1
% beta =   1.299989776988072
% step =     2
% beta =   1.424389879545184
% step =     3
% beta =   1.424228404569439
% ans  = Probability of failure: 0.0771902
% % % % % % % % % % % % % % % % % % % 

% u0=11
% step =     1
% beta =   0.722846633012398
% step =     2
% beta =   0.760648877360251
% step =     3
% beta =   0.760640830412670
% ans = Probability of failure: 0.223436
% % % % % % % % % % % % % % % % % % % % %

% u0=12
% step =     1
% beta =   0.145703489036724
% step =     2
% beta =   0.147204861176776
% step =     3
% beta =   0.147204858905609
% ans = Probability of failure: 0.441485
% % % % % % % % % % % % % % % % % % % % %

% u0=12.5
% step =     1
% beta =   0.142868082951114
% step =     2
% beta =   0.141425152878370
% step =     3
% beta =   0.141425151032474
% ans = Probability of failure: 0.443767
% % % % % % % % % % % % % % % % % % 
% 
% u0=13.5
% step =     1
% beta =   0.720011226926788
% step =     2
% beta =   0.684301082594260
% step =     3
% beta =   0.684294433686579
% ans = Probability of failure: 0.246895
% % % % % % % % % % % % % % % % % % 

% u0=14
% step =     1
% beta =   1.008582798914625
% step =     2
% beta =   0.939353129088658
% step =     3
% beta =   0.939318500388323
% ans =Probability of failure: 0.173784
% % % % % % % % % % % % % % % % % % % 
