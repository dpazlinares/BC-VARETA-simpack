%% cleaning
clear all;
clc;
close all;

%% loading data
load('data/CUBAN_DATA_6k_15k_INFO.mat');
load('data/mycolormap_brain_basic_conn.mat');
[Ne1,Np1] = size(K_6k);
Nseg      = 1e5;
snr       = 0.90;
snr_ch    = 0.90;

%% one dipole
temp_c = zeros(Np1);
temp_c(1660,1660) = 1;
title_name = 'One active dipole';

%% two dipoles
temp_c = zeros(Np1);
dipole1 = 3920;
dipole2 = 5309;
% dipole3 = 2803;
dipole3 = 2713;
temp_c(dipole1,dipole1) = 1;
temp_c(dipole2,dipole2) = 1;
title_name = 'Two active dipoles';

d = zeros(length(ASA_343.Channel),1);
s_fix = S_6k.Vertices(3920,:);
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
s_fix = S_6k.Vertices(5309,:);
for ii = 1:length(ASA_343.Channel)
    d(ii) = (sum([(s_fix(1)-ASA_343.Channel(ii).Loc(1))^2,...
                 (s_fix(2)-ASA_343.Channel(ii).Loc(2))^2,...
                 (s_fix(3)-ASA_343.Channel(ii).Loc(3))^2]))^0.5;
end
d = d/max(d);
plot(d);
[B,I2] = sort(d);

It = [I1(1:12);  I2(1:12)];
Itt = unique(It,'stable');


%% three dipoles
temp_c = zeros(Np1);
dipole1 = 84;
dipole2 = 527;
% dipole3 = 2803;
dipole3 = 2713;
temp_c(dipole1,dipole1) = 1;
temp_c(dipole2,dipole2) = 1;
temp_c(dipole3,dipole3) = 1;
% temp_c(dipole1,dipole2) = 1;
% temp_c(dipole2,dipole1) = 1;
% title_name = 'Three active dipoles with connection between 2 and 3';
title_name = 'Three active dipoles with no connection';

%%
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

% It = [I1([10 8 6 4 2 1 3 5 7 9 11]);...
%       I2([10 8 6 4 2 1 3 5 7 9 11]);...
%       I3([10 8 6 4 2 1 3 5 7 9 11])];
It = [I3([6 4 2 1 3 5]);...
      I2([4 5 3 1 2 7]);...
      I1([4 2 1 3 5 6 7])];
Itt = unique(It,'stable');
% Itt = It;


%% forward problem
Svv = K_6k*temp_c*K_6k';
Svv_inv = pinv(Svv);

%% biological noise, mixture of independent pink noise sources...
n_noise_sources   = round(0.03*Np1);
ind_noise_rand    = randperm(Np1);
ind_noise         = ind_noise_rand(1:n_noise_sources)';
pn                = mkpinknoise(Nseg, n_noise_sources)';
norm_brain_noise  = norm(pn, 'fro');
Sjj_brain_noise   = K_6k(:, ind_noise)*pn;
Sjj_brain_noise   = Sjj_brain_noise ./ norm_brain_noise;
[Svv_brain_noise] = xspectrum(Sjj_brain_noise,[],[],[],1);
Svv_brain_noise   = Svv_brain_noise(:,:,1);
Svv_brain_noise   = Svv_brain_noise ./ norm(Svv_brain_noise(:), 'fro');

%% white electrode noise...
Svv_channel_noise   = randn(Ne1,Nseg);
Svv_channel_noise   = Svv_channel_noise ./ norm(Svv_channel_noise, 'fro');
[Svv_channel_noise] = xspectrum(Svv_channel_noise,[],[],[],1);
Svv_channel_noise   = Svv_channel_noise(:,:,1);
Svv_channel_noise   = Svv_channel_noise ./ norm(Svv_channel_noise(:), 'fro');
Svv                 = Svv ./ norm(Svv(:), 'fro');
Svv                 = snr*Svv + (1-snr)*Svv_brain_noise;
Svv                 = Svv ./ norm(Svv, 'fro');
Svv                 = snr_ch*Svv + (1-snr_ch)*Svv_channel_noise;
Svv                 = Svv ./ norm(Svv, 'fro');

%% topo
cfg.layout     = 'EEG1005.lay';
cfg.channel    = 'eeg';
eeg.label      = label;
eeg.dimord     = 'chan_freq';
eeg.freq       = 1;
eeg.powspctrm  = sqrt(C);
figure('Color','k'); 
ft_topoplotER(cfg, eeg);
colormap(gca,cmap);
set(gca,'Color',[0.5 0.5 0.5]);

%% ploting cross-spectrum

% figure('Color','k');
% imagesc((abs(Svv(Itt,Itt))/max(max(abs(Svv(Itt,Itt))))).^0.4);
% set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
%     'XTickLabel',{ASA_343.Channel(Itt).Name},'XTickLabelRotation',45,...
%     'YTickLabel',{ASA_343.Channel(Itt).Name},'YTickLabelRotation',45); 
% xlabel('electrodes','Color','w');
% ylabel('electrodes','Color','w');
% colormap(gca,'jet');
% title('Synthetic Svv','Color','w','FontSize',16);

figure('Color','k');
imagesc((abs(Svv_inv(Itt,Itt))/max(max(abs(Svv_inv(Itt,Itt))))).^0.4);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTickLabel',{ASA_343.Channel(Itt).Name},'XTickLabelRotation',45,...
    'YTickLabel',{ASA_343.Channel(Itt).Name},'YTickLabelRotation',45); 
xlabel('electrodes','Color','w');
ylabel('electrodes','Color','w');
colormap(gca,'jet');
title('Inverse Svv','Color','w','FontSize',16);

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
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',0.1*rand(length(S_6k.Vertices),1)+diag(temp_c),'FaceColor','interp','EdgeColor','none','FaceAlpha',.89);
colormap(gca,cmap);
title(title_name,'Color','w','FontSize',16);

%% eloreta
A = mkfilt_eloreta2(K_h, 4.05);
h = A'*Svv*A;
temp_h = diag(abs(h));
temp_h = temp_h/max(temp_h);
temp_h(temp_h < 0.1) = 0;
figure('Color','k'); hold on; set(gca,'Color','k');
patch('Faces',S_h.Faces,'Vertices',S_h.Vertices,'FaceVertexCData',temp_h,...
    'FaceColor','interp','EdgeColor','none','FaceAlpha',.75);
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',diag(abs(temp_c)),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
colormap(gca,cmap);
title(title_name,'Color','w','FontSize',16);

%% bc-vareta
Shh = {Svv};
K = {K_h};
[Ne,Np] = size(K{1});
groups    = [];
for ii = 1:Np
    groups{ii} = ii;
end
[miu,alpha1,alpha2,h] = cross_nonovgrouped_enet_ssbl(Shh,K,Nseg,groups);
temp_h = h;
temp_h = temp_h/max(temp_h);
temp_h(temp_h < 0.1) = 0;
figure('Color','k'); hold on; set(gca,'Color','k');
patch('Faces',S_h.Faces,'Vertices',S_h.Vertices,'FaceVertexCData',(abs(temp_h)),...
    'FaceColor','interp','EdgeColor','none','FaceAlpha',.75);
patch('Faces',S_6k.Faces,'Vertices',S_6k.Vertices,'FaceVertexCData',temp_c,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
colormap(gca,cmap);
title(title_name,'Color','w','FontSize',16);