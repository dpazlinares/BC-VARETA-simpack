function [miu,alpha1,alpha2,h] = screening_ssbl(Svv,LeadField,Nsamp,nonovgroups,vertices,faces)
% Elastic Net_Sparse Bayesian Learning
%% Loading Dimensions
Nv          = size(LeadField,2);                                           % Number of nodes
Ngroups     = length(nonovgroups);                                         % Number of spatial groups
sigma       = ones(Nv,1);                                                  % Prior Variances
sigma_post  = ones(Nv,1);                                                  % Posterior Variances
h           = ones(Nv,1);                                                  % Hypeparameter for Gamma
%% User defined parameters
maxiter1    = 20;                                                          % Number of Iterations of outer cycle
maxiter11   = 3;                                                           % Number of Iterations of inner cycle
e1          = 1E-5;                                                        % Rate parameter of the gamma pdf of alpha (higher -> smoother)
e2          = 1E-6;                                                        % Rate parameter of the gamma pdf of k (higher -> smoother)
beta        = 1;                                                           % Data variance
%% Initialization of parameters
K           = LeadField;
Kt          = transpose(K);
Ne          = size(K,1);
Ine         = spdiags(ones(Ne,1),0,Ne,Ne);                                 % Identity matrix with Ne size
scale_K     = sqrt(trace(K*Kt)/Ne);                                        % Lead Field Scale Norm inf
K           = K/scale_K;                                                   % Lead Field scaling
Kt          = Kt/scale_K;
%% Calibration Inverse Solution
sigmaKt     = spdiags(sigma,0,Nv,Nv)*Kt;
sigma_post1 = sigmaKt/(K*sigmaKt+beta*Ine);
sigma_post2 = K*spdiags(sigma,0,Nv,Nv);
for jj=1:Nv
    sigma_post(jj) = sigma(jj)-sigma_post1(jj,:)*sigma_post2(:,jj);
end
TrOp        = (1/beta).*(sigmaKt - sigma_post1*(sigma_post2*Kt));
miu_cal     = diag(TrOp*Svv*TrOp');
scale_V     = mean(abs(miu_cal))/max(sigma_post);
%% Data Scaling
Svv         = Svv/scale_V;
%% Initialization of Hyperparameters
alpha1      = 1E0;                                                         % Hyperparameter of the L2 norm
k           = 1E0;                                                         % Hyperparameter of the Truncate Gamma pdf
%% Main Cycle
for cont1 = 1:maxiter1
    for cont11 = 1:maxiter11
        %% Update Posterior Mean and Covariance matrix
        sigmaKt     = spdiags(sigma,0,Nv,Nv)*Kt;
        sigma_post1 = sigmaKt/(K*sigmaKt+beta*Ine);
        sigma_post2 = K*spdiags(sigma,0,Nv,Nv);
        % Only save the diagonals of the Posterior Covariance
        for jj=1:Nv
            sigma_post(jj) = sigma(jj)-sigma_post1(jj,:)*sigma_post2(:,jj);
        end
        % Iterative Transference Operator
        TrOp        = (1/beta).*(sigmaKt-sigma_post1*(sigma_post2*Kt));
        % Compute 'miu'
        miu2        = abs(diag(TrOp*Svv*TrOp'));
        %% Update Gammas
        for group = 1:Ngroups
            idx_group    = nonovgroups(group);
            idx_group    = idx_group{1};
            h(idx_group) = sqrt((1./4).^2+alpha1.*k.*Nsamp.*sum(miu2(idx_group) + sigma_post(idx_group)))-1./4;
        end
        index_h        = find((miu2 + sigma_post)<0);
        h(index_h)     = 0;
        gamma          = k + h;
        sigma_bar      = h./gamma;
        sigma          = (1/(2*alpha1))*sigma_bar;
    end
    %% Update alpha_1
    index_alpha_1  = find(sigma_bar>0);
    alpha1         = (length(index_alpha_1)/2 + Nv*Nsamp)/(Nsamp*sum((miu2(index_alpha_1) + sigma_post(index_alpha_1))./(sigma_bar(index_alpha_1))) + e1*Nv*Nsamp);
    %% Update alpha_2
    f_aux          = @(k_aux) e2*Nsamp + sum(ones(Nv,1)./(1-sigma_bar))/Nv - (Nsamp - 1/2)/k_aux - trascendent_term(k_aux);
    k              = fzero(f_aux,[0.000001 700000]);
    alpha2         = (4*alpha1*k)^(1/2);
    sigma          = (1/(2*alpha1))*sigma_bar;
end
miu = (miu2*scale_V/scale_K);
% Plot_sources_Haufe(miu2/max(abs(miu2)),vertices,faces,'simple')
% caxis([0 0.1])
end