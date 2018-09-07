function [ThetaJJ,SJJ,indms] = bcvareta(Svv,Kfull,m,vertices,faces)
% BC-VARETA  computes the Source Activity and Connectivity from the frequency 
% domain representation of MEEG stationary time series as the solution of embeded
% Gaussian Graphical Models by the EM algorithm with Sources Hermitian Graphical LASSO
%
% inputs:
%    Svv           : MEEG time series cross-spectra
%    K             : MEEG Lead Fields
%    Nsamp         : Sample number of MEEG cross-spectra
%
% outputs:
%    WJJ           : source covariance matrix
%    ThetaJJ       : source precision matrix
%    PsiJJ         : source empirical covariance matrix
%
% Pedro Valdes-Sosa, September 2018
% Deirel Paz Linares, September 2018
% Eduardo Gonzalez-Moreira, September 2018
%**************************************************************************
%% step 1 Initialization EM algorithm
[sigmae2,epsilon,p,Ip,SigmaVV,SigmaVVinv,maxiter1,maxiter11,K,Kt,indms,WJJ,Svv,lambda2,rho,mask1,mask2] = bcvareta_initial_values(Svv,Kfull,m,vertices,faces);
%% Outer loop EM algorithm
for cont1 = 1:maxiter1
    disp(['outer loop iteration # ',num2str(cont1)])
    %% Source Posterior Covariance Matrix
    WJJKt            = WJJ*Kt; % Auxiliar Matrix (Sources Covariance Matrix)*(Lead Field)'
    SigmaJJ          = WJJ-WJJKt*pinv(K*WJJKt+sigmae2*SigmaVV)*K*WJJ; % Woodbury formula
    SigmaJJ(mask2)   = 0;
    %%
    %% Data to Source Transference Operator
    T                = SigmaJJ*Kt*SigmaVVinv*(1/sigmae2); % Transference Operator
    Tt               = T'; % Tranconjugated Transference Operator
    A                = (Ip - K*T); % Auxiliar Matrix  I - (Lead Field)*(Transference Operator)
    %%
    %% Effective Source Empirical Covariance Matrix
    SJJ              = T*Svv*Tt; % Source Empirical Covariance Matrix
    PsiJJ            = SigmaJJ + SJJ; % Effective Source Empirical Covariance Matrix
    WJJ              = PsiJJ; % Naive Source Covariance Matrix estimator 
    %%
    %% Source Graphical Model Local Quadratic Approximation
    [ThetaJJ]        = sggm_lqa(PsiJJ,m,lambda2,rho,maxiter11); % Source Precision Matrix estimator
    %%
    %% Unbiased Source Precision Matrix estimator
    ThetaJJ          = 2*ThetaJJ - ThetaJJ*PsiJJ*ThetaJJ; % Unbiased Source Precision Matrix estimator
    varThetaJJ       = abs(diag(ThetaJJ))*abs(diag(ThetaJJ))' + abs(ThetaJJ).^2; % Consistent variances Unbiased Source Precision Matrix estimator    
    ThetaJJ(mask2)   = 0;
    %% Threshold mask given the Unbiased Source Precision Matrix Normal tendency
    [mask2]          = find(abs(ThetaJJ).^2 < (0.05/log(m))*(varThetaJJ - diag(diag(varThetaJJ)))); % "ThetaJJ < 0.05varThetaJJ" Thresholding mask
    ThetaJJ(mask2)   = 0;
    %%
    %% Noise Variance estimator
    sigmae2                = trace(SigmaVVinv*A*Svv*A')/p + trace(Kt*SigmaVVinv*K*SigmaJJ)/p + epsilon;
    %%
end %iterations outer loop
end

%% Source Graphical Model Local Quadratic Approximation
function [PM] = sggm_lqa(EC,m,lambda2,rho,maxiter)
%% sggm-lqa initial values
m2          = m^2;                                                   
m12         = m^(1/2);                                               
nu          = 2.5e0;
EC          = EC + 1e-4*max(abs(diag(EC)))*eye(length(EC));
ECinv       = pinv(EC);
ECinv_ph    = exp(1i*angle(ECinv));
idx         = (lambda2 > 0);
idx0        = (lambda2 == 0);
gamma       = zeros(length(lambda2));
gamma(idx0) = rho*m*abs(ECinv(idx0)).^2;
PM          = ECinv;
for cont = 1:maxiter
    %% Estimation of variances Gamma of Gaussian Mixtures prior
    DET           = 1 + 4*m2*lambda2(idx).*abs(PM(idx)).^2;
    gamma(idx)    = (sqrt(DET) - 1)./(2*m*lambda2(idx));
    ninf          = max(gamma(:));
    gamman        = gamma/ninf;
    %%
    %% Standarization of the Empirical Covariance Matrix     
    st_factor1    = ninf^(-1/2)*ECinv.*gamman.^(1/2)*(nu*m12) + (1/(m12))*ECinv_ph; % Consistent factor1
    st_factor2    = (1 + gamman*(nu*m12)); % Consitent factor2
    ECst_inv      = st_factor1./st_factor2; % Inverted Standard Empirical Covariance Matrix
    ECst          = pinv(ECst_inv); % Standard Empirical Covariance Matrix
    %%
    %% Standard Precision Matrix estimator
    [U,D,V]       = svd(ECst);
    D             = abs(diag(D));
    PMst          = (1/2)*U*diag(sqrt(D.^2 + 4) - D)*V';
    PMst          = (PMst + PMst')/2;
    %%
    %% Unstandarized Precision Matrix estimator
    PM          = gamma.^(1/2).*PMst;
    %%
end
end
