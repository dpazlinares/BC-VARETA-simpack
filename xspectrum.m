function [Svv,F,Ns] = xspectrum(data,Fs,Fm,deltaf,plotPSD)
% xspectrum estimates the Cross Spectrum of the input M/EEG data
% Inputs:
%    data     = M/EEG data matrix, in which every row is a channel
%    Fs       = sampling frequency (in Hz), default = 200
%    Fm       = maximun frequency (in Hz) in the estimated spectrum, default = 19
%    deltaf   = frequency resolution, default = 0.3906
%    plotPSD  = plotting flag, 1 plots estimated PSD, 0 doesn't (default)
% Outputs:
%    PSD      = estimated power spectral density of input EEG data
%    Svv      = estimated cross spectrum of input EEG data
%    Ns       = number of segments in which the EEG signal is wrapped
%
%
%% 
% =============================================================================
% This function is part of the BC-VARETA toolbox:
% https://github.com/egmoreira/BC-VARETA-toolbox
% =============================================================================@
%
% Authors:
% Pedro A. Valdes-Sosa, 2010-2018
% Alberto Taboada-Crispi, 2016
% Deirel Paz-Linares, 2017-2018
% Eduardo Gonzalez-Moreira, 2017-2018
%
%**************************************************************************
%% Initialization oF variables...
if (nargin < 5) || isempty(plotPSD)
    plotPSD = 0;                                        % plotting flag
end
if (nargin < 4) || isempty(deltaf)
    deltaf   = 0.3906;                                  % frequency resolution
end
if (nargin < 3) || isempty(Fm)
    Fm       = 19;                                      % maximun frequency in the estimated spectrum, in Hz
end
if (nargin<2) || isempty(Fs)
    Fs = 200;                                           % sampling frequency in Hz
end
NFFT     = round(Fs/deltaf);                            % number of time points per window
Nw       = 2;                                           % number of windows for Thomson spectral estimate
F        = 0:deltaf:Fm;                                 % frequency vector
%% Estimation of the Cross Spectrum...
e       = dpss(NFFT,Nw);                                % discrete prolate spheroidal (Slepian) sequences
e       = reshape(e,[1,NFFT,2*Nw]);
[Nc,Ns] = size(data);                                   % number of channels (rows) and samples (columns)
Ns      = fix(Ns/NFFT);                                 % number of segments in which the EEG signal is wrapped
Ns      = max(1,Ns);
if NFFT > Ns
    data = [data zeros(Nc,NFFT-Ns)];                    % zero padding
end
data(:,Ns*NFFT+1:end) = [];                             % discards samples after Ns*NFFT
data = reshape(data,Nc,NFFT,Ns);                        % 'resized' EEG data
lf = length(F);                                         % length of F vector
Svv = zeros(Nc,Nc,lf);                                  % allocated matrix for the cross spectrum
for k = 1:Ns
    w = data(:,:,k);                                    % k-th window
    W = repmat(w,[1,1,2*Nw]).*repmat(e,[Nc,1,1]);       % multiplied by Slepian seq
    W = fft(W,[],2);                                    % FFT
    W = W(:,1:lf,:);                                    % pruning values of the FFT
    for i=1:lf
        Svv(:,:,i) = Svv(:,:,i)+cov(squeeze(W(:,i,:)).',1);
    end
end
Svv = Svv/Ns;                                           % normalizing
%% Estimation of Power Spectral Density (PSD)...
PSD = zeros(Nc,lf);                                     
for freq = 1:lf
    PSD(:,freq) = diag(squeeze(real(Svv(:,:,freq))));
end
%% Plotting PSD...
if plotPSD
    figure('Color','k');
    plot(F,10*log10(abs(PSD)));
    set(gca,'Color','k','XColor','w','YColor','w');
    ylabel('PSD (dB)','Color','w');
    xlabel('Freq. (Hz)','Color','w');
    title('Power Spectral Density','Color','w');
end
end