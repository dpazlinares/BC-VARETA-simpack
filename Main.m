%% BC-VARETA main

%% cleaning
clear all;
clc;
close all;

%% loading data
load('data/CUBAN_DATA_6k_15k_INFO.mat');
load('data/mycolormap_brain_basic_conn.mat');
% load('data/mycolormap.mat');
[Ne1,Np1] = size(K_6k);
Nseg      = 5e2;
snr       = 1;
snr_ch    = 1;

%% three dipoles
Sjj  = zeros(Np1);
% dipole1 = 19;
% dipole2 = 209;
% dipole3 = 205;

% dipole1 = 209;
% dipole2 = 19;
% dipole3 = 205;

dipole1 = 19;
dipole2 = 1956;
dipole3 = 1287;
dipole4 = 2945; 

Thetajj = [1 1 0 0;1 1 0 0;0 0 1 1;0 0 1 1];
Wjj     = pinv(Thetajj);

dipoles = [dipole1 dipole2 dipole3 dipole4];

Sjj(dipoles,dipoles) = Wjj;

title_name = 'Three active dipoles with connection between 1 and 2';

%% active indices for brain
indact  = [];
Neibo   = 14;
d0      = 20;

[index] = surfpatch_v1(dipole1,S_6k.Vertices,S_6k.Faces,d0,Neibo);
if length(index) >= Neibo
    index = index(1:Neibo);
end
index1 = index([14 12 10 8 6 2 4 1 3 5 7 9 11 13]);
indact  = [indact; index1];

[index] = surfpatch_v1(dipole2,S_6k.Vertices,S_6k.Faces,d0,Neibo);
index = setdiff(index,indact,'stable');
if length(index) >= Neibo
    index = index(1:Neibo);
end
index2 = index([14 12 10 8 6 2 4 1 3 5 7 9 11 13]);
indact  = [indact; index2];

[index] = surfpatch_v1(dipole3,S_6k.Vertices,S_6k.Faces,d0,Neibo);
index = setdiff(index,indact,'stable');
if length(index) >= Neibo
    index = index(1:Neibo);
end
index3 = index([14 10 8 2 4 12 6 1 3 5 7 9 11 13]);
indact = [indact; index3];

[index] = surfpatch_v1(dipole4,S_6k.Vertices,S_6k.Faces,d0,Neibo);
index = setdiff(index,indact,'stable');
if length(index) >= Neibo
    index = index(1:Neibo);
end
index4 = index([14 10 8 2 4 12 6 1 3 5 7 9 11 13]);
indact = [indact; index4];

indlab = [];
for ii = 1:length(S_6k.Vertices)
    indlab{ii} = ii;
end

%% active indices for scalp
d = zeros(length(ASA_343.Channel),1);
s_fix = S_6k.Vertices(dipole1,:);
for ii = 1:length(ASA_343.Channel)
    d(ii) = (sum([(s_fix(1)-ASA_343.Channel(ii).Loc(1))^2,...
                  (s_fix(2)-ASA_343.Channel(ii).Loc(2))^2,...
                  (s_fix(3)-ASA_343.Channel(ii).Loc(3))^2]))^0.5;
end
d = d/max(d);
figure; hold on;
plot(d);
[B,I1] = sort(d);

d = zeros(length(ASA_343.Channel),1);
s_fix = S_6k.Vertices(dipole2,:);
for ii = 1:length(ASA_343.Channel)
    d(ii) = (sum([(s_fix(1)-ASA_343.Channel(ii).Loc(1))^2,...
                  (s_fix(2)-ASA_343.Channel(ii).Loc(2))^2,...
                  (s_fix(3)-ASA_343.Channel(ii).Loc(3))^2]))^0.5;
end
d = d/max(d);
plot(d);
[B,I2] = sort(d);

d = zeros(length(ASA_343.Channel),1);
s_fix = S_6k.Vertices(dipole3,:);
for ii = 1:length(ASA_343.Channel)
     d(ii) = (sum([(s_fix(1)-ASA_343.Channel(ii).Loc(1))^2,...
                   (s_fix(2)-ASA_343.Channel(ii).Loc(2))^2,...
                   (s_fix(3)-ASA_343.Channel(ii).Loc(3))^2]))^0.5;
end
d = d/max(d);
plot(d);
[B,I3] = sort(d);

d = zeros(length(ASA_343.Channel),1);
s_fix = S_6k.Vertices(dipole4,:);
for ii = 1:length(ASA_343.Channel)
     d(ii) = (sum([(s_fix(1)-ASA_343.Channel(ii).Loc(1))^2,...
                   (s_fix(2)-ASA_343.Channel(ii).Loc(2))^2,...
                   (s_fix(3)-ASA_343.Channel(ii).Loc(3))^2]))^0.5;
end
d = d/max(d);
plot(d);
[B,I4] = sort(d);


% It = [I1([10 8 6 4 2 1 3 5 7 9 11]);...
%       I2([10 8 6 4 2 1 3 5 7 9 11]);...
%       I3([10 8 6 4 2 1 3 5 7 9 11])];
It = [I1([6 4 2 1 3 5 7]);...
      I2([4 5 3 1 2 7 6]);...
      I3([4 2 1 3 5 6 7]);
      I4([4 2 1 3 5 6 7])];
indele = unique(It,'stable');
% Itt = It;

%% sources space

sources_iv = abs(diag(Sjj));
connect_iv = abs(Sjj);
figure('Color','k'); hold on;
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',sources_iv,'FaceColor','interp','EdgeColor','none','FaceAlpha',.85);
set(gca,'Color','k','XColor','k','YColor','k','ZColor','k');
colormap(gca,cmap);
title('Ground truth','Color','w','FontSize',16);
figure('Color','k');
imagesc((abs(connect_iv(indact,indact))/max(max(abs(connect_iv(indact,indact))))).^0.99);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(indact),'YTick',1:length(indact),...
    'XTickLabel',indlab(indact),'XTickLabelRotation',90,...
    'YTickLabel',indlab(indact),'YTickLabelRotation',0); 
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colormap(gca,cmap);
axis square;
title('Ground truth','Color','w','FontSize',16);

%% forward problem
Svv = K_6k*Sjj*K_6k';

%% biological noise, mixture of independent pink noise sources...
% n_noise_sources   = round(0.03*Np1);
% ind_noise_rand    = randperm(Np1);
% ind_noise         = ind_noise_rand(1:n_noise_sources)';
% pn                = mkpinknoise(Nseg, n_noise_sources)';
% norm_brain_noise  = norm(pn, 'fro');
% Sjj_brain_noise   = K_6k(:, ind_noise)*pn;
% Sjj_brain_noise   = Sjj_brain_noise ./ norm_brain_noise;
% [Svv_brain_noise] = xspectrum(Sjj_brain_noise,[],[],[],1);
% Svv_brain_noise   = Svv_brain_noise(:,:,1);
% Svv_brain_noise   = Svv_brain_noise ./ norm(Svv_brain_noise(:), 'fro');

%% white electrode noise...
% Svv_channel_noise   = randn(Ne1,Nseg);
% Svv_channel_noise   = Svv_channel_noise ./ norm(Svv_channel_noise, 'fro');
% [Svv_channel_noise] = xspectrum(Svv_channel_noise,[],[],[],1);
% Svv_channel_noise   = Svv_channel_noise(:,:,1);
% Svv_channel_noise   = Svv_channel_noise ./ norm(Svv_channel_noise(:), 'fro');

%% updating Svv with noise
% Svv                 = Svv ./ norm(Svv(:), 'fro');
% Svv                 = snr*Svv+(1-snr)*Svv_brain_noise;
% Svv                 = Svv ./ norm(Svv, 'fro');
% Svv                 = snr_ch*Svv + (1-snr_ch)*Svv_channel_noise;
% Svv                 = Svv ./ norm(Svv, 'fro');

%% inverse covariance matrix 
Svv_inv = pinv(Svv);

%% electrodes space
X = zeros(length(ASA_343.Channel),1);
Y = zeros(length(ASA_343.Channel),1);
Z = zeros(length(ASA_343.Channel),1);
label = [];
for ii = 1:length(ASA_343.Channel)
    X(ii) = ASA_343.Channel(ii).Loc(1);
    Y(ii) = ASA_343.Channel(ii).Loc(2);
    Z(ii) = ASA_343.Channel(ii).Loc(3);
    label{ii} = ASA_343.Channel(ii).Name;
end
C = abs(diag(Svv));
C = C/max(C); 
C(C<0.01) = 0;

figure('Color','k'); hold on; set(gca,'Color','k');
scatter3(X,Y,Z,100,C.^1,'filled');
patch('Faces',S_h.Faces,'Vertices',S_h.Vertices,'FaceVertexCData',0.01*(ones(length(S_h.Vertices),1)),'FaceColor','interp','EdgeColor','none','FaceAlpha',.35);
colormap(gca,cmap);
title('Scalp','Color','w','FontSize',16);

figure('Color','k');
imagesc((abs(Svv_inv(indele,indele))/max(max(abs(Svv_inv(indele,indele))))).^0.5);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length({ASA_343.Channel(indele).Name}),'YTick',1:length({ASA_343.Channel(indele).Name}),...
    'XTickLabel',{ASA_343.Channel(indele).Name},'XTickLabelRotation',90,...
    'YTickLabel',{ASA_343.Channel(indele).Name},'YTickLabelRotation',0); 
xlabel('electrodes','Color','w');
ylabel('electrodes','Color','w');
colormap(gca,'jet');
axis square;
title('Scalp','Color','w','FontSize',16);

%% eloreta
P          = mkfilt_eloreta(K_6k,1);
P          = P';
connect_iv = P*Svv*P';
sources_iv = abs(diag(connect_iv));
sources_iv = sources_iv/max(sources_iv(:));
connect_iv = pinv(connect_iv);
% connect_iv = connect_iv-diag(diag(connect_iv));
connect_iv = abs(connect_iv)/max(abs(connect_iv(:)));
sources_iv(sources_iv < 0.01) = 0;
connect_iv(connect_iv < 0.01) = 0;
figure('Color','k'); hold on;
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',sources_iv,'FaceColor','interp','EdgeColor','none','FaceAlpha',.85);
set(gca,'Color','k');
colormap(gca,cmap);
title('eLORETA','Color','w','FontSize',16);
figure('Color','k');
imagesc((abs(connect_iv(indact,indact))/max(max(abs(connect_iv(indact,indact))))).^0.5);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(indact),'YTick',1:length(indact),...
    'XTickLabel',indlab(indact),'XTickLabelRotation',90,...
    'YTickLabel',indlab(indact),'YTickLabelRotation',0); 
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colormap(gca,'jet');
axis square;
title('eLORETA','Color','w','FontSize',16);

%% lcmv
P          = mkfilt_lcmv(K_6k,Svv,3e4);
P          = P';
connect_iv = P*Svv*P';
sources_iv = abs(diag(connect_iv));
sources_iv = sources_iv/max(sources_iv(:));
connect_iv = pinv(connect_iv);
% connect_iv = connect_iv-diag(diag(connect_iv));
connect_iv = abs(connect_iv)/max(abs(connect_iv(:)));
sources_iv(sources_iv < 0.01) = 0;
connect_iv(connect_iv < 0.01) = 0;

figure('Color','k'); hold on;
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',sources_iv,'FaceColor','interp','EdgeColor','none','FaceAlpha',.85);
set(gca,'Color','k');
colormap(gca,cmap);
title('LCMV','Color','w','FontSize',16);
figure('Color','k');
imagesc((abs(connect_iv(indact,indact))/max(max(abs(connect_iv(indact,indact))))).^0.5);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(indact),'YTick',1:length(indact),...
    'XTickLabel',indlab(indact),'XTickLabelRotation',90,...
    'YTickLabel',indlab(indact),'YTickLabelRotation',0); 
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colormap(gca,'jet');
axis square;
title('LCMV','Color','w','FontSize',16);

%% bc-vareta
[ThetaJJ,SJJ,indms] = bcvareta(Svv,K_6k,Nseg,S_6k.Vertices,S_6k.Faces);

sources_iv         = zeros(length(K_6k),1);
sources_iv(indms)  = abs(diag(SJJ));
sources_iv         = sources_iv/max(sources_iv(:));
ind_zr = sources_iv < 0.01;
sources_iv(ind_zr) = 0;

figure('Color','k'); hold on;
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',sources_iv.^0.5,'FaceColor','interp','EdgeColor','none','FaceAlpha',.95);
set(gca,'Color','k');
colormap(gca,cmap);
title('SJJ diagonal','Color','w','FontSize',16);


connect_iv              = zeros(length(K_6k));
connect_iv(indms,indms) = abs(ThetaJJ);
connect_iv              = connect_iv/max(connect_iv(:));
% connect_iv(ind_zr,:) = 0;
% connect_iv(:,ind_zr) = 0;

figure('Color','k');
imagesc(abs(connect_iv(indact,indact)).^0.5);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(indact),'YTick',1:length(indact),...
    'XTickLabel',indlab(indact),'XTickLabelRotation',90,...
    'YTickLabel',indlab(indact),'YTickLabelRotation',0); 
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colormap(gca,cmap);
axis square;
title('ThetaJJ','Color','w','FontSize',16);
