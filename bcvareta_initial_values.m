function [sigmae2,epsilon,p,Ip,SigmaVV,SigmaVVinv,maxiter1,maxiter11,K,Kt,indms,WJJ,Svv,lambda2,rho,mask1,mask2] = bcvareta_initial_values(Svv,Kfull,m,vertices,faces)
%% Screening Structured Sparse Bayesian Learning
[p,qfull]      = size(Kfull);
groups         = [];
for ii = 1:qfull
    groups{ii} = ii;
end
indana         = 1:qfull;
[miu]          = screening_ssbl(Svv,Kfull,m,groups);
[indms]        = screening(miu,0.01*qfull,vertices,faces,indana);
K              = Kfull(:,indms);
%% Initialization of variables and tunning parameters
maxiter1       = 300;                                                        % Maximum number of outer EM loop iterations
maxiter11      = 30;                                                        % Maximum number of inner EM loop iterations
epsilon        = 5e0;                                                       % Hyperparmeter of the data nuisance prior
sigmae2        = 1E0;                                                       % Hyperparameter of nuisance initial value
q              = size(K,2);                                                 % Number of cortical generators
Ip             = eye(p);                                                    % Identity matrix on the sensors space
Iq             = eye(q);                                                    % Identity matrix on the cortical generators space
rho            = sqrt(log(q)/m);                                            % Regularization parameter
rho_diag       = 0;                                                         % Regularization mask diagonal
rho_ndiag      = 1;                                                         % Regularization mask nondiagonal
lambda2        = (rho*(rho_diag*eye(q)+rho_ndiag*(ones(q)-eye(q)))).^2;     % Regularization mask squared
mask1          = [];                                                        % "ThetaJJ == 0" Thresholding mask                                             
mask2          = [];                                                        % "ThetaJJ < 0.005varThetaJJ" Thresholding mask
WJJ            = Iq;                                                        % Covariance matrix of cortical activity
SigmaVV        = Ip;                                                        % Covariance matrix of sensors activity
SigmaVVinv     = Ip/SigmaVV;                                                % Precision matrix of sensors activity
%% Generators Covariance matrix initialization and Lead Field and Data scaling
Kt             = K';
scale_K        = sqrt(trace(K*Kt)/p);
K              = K/scale_K;
Kt             = Kt/scale_K;
%% Transference operator and posterior covariance
WJJKt          = WJJ*Kt;
SigmaJJ        = WJJ-(WJJKt/(K*WJJKt+sigmae2*SigmaVV))*K*WJJ;
T              = SigmaJJ*Kt*SigmaVVinv*(1/sigmae2);
%% Computation of sources activity calibration empirical covariance
SJJ            = T*Svv*T';
%% Scaling the sources activity empirical covariance and data empirical covariance
scale_SJJ      = (trace(SJJ)/q)/max(diag(SigmaJJ));
Svv            = Svv/scale_SJJ;
end