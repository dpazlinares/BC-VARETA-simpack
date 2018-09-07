function noise=mkpinknoise(n,m);
% makes m channels of pink noise
% of length n. Each column is one 
% channel
%
% Guido Nolte, 2012-2015
% g.nolte@uke.de

% If you use this code for a publication, please ask Guido Nolte for the
% correct reference to cite.

% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.



% randn('state',sum(100*clock))
n1=2*ceil((n-1)/2)+1;
scal=sqrt(1./[1:(n1-3)/2]');
ff=zeros(n1,m);
ff(2:(n1-1)/2,:)=repmat(scal,1,m);
noise=fft(randn(n1,m));
noise=2*real(ifft(noise.*ff));
noise=noise(1:n,:);

return;


