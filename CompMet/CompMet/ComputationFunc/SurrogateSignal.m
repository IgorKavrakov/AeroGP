function [Sur,X,Y,f] = SurrogateSignal(x,dt,n,FFTflag)
% Generate Fourier surrogates
% This file is part of CompMet.
% CompMet is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%  CompMet is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with CompMet.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020

%% Input
% x is the original signal
% n is the number of suruogate signals to be crated
% dt is the time step
% FFTflag is a flag wheter to compute the ffts
% beta is the padding operation (NSig/beta) is padded on both sides
%% Test
% clear all;close all;clc;
% dt=0.01;
% NSamples=10001;
% t=(1:1:NSamples)*dt;
% XNoise=1.*cos(2.*pi.*2.*t)+1.4.*cos(2.*pi.*3.*t);
% noiseSamp=random('norm',0,std(XNoise)./5,[1 NSamples]);
% x=1.*cos(2.*pi.*2.*t)+1.4.*cos(2.*pi.*3.*t+pi/3)+noiseSamp;
% n=5;
% FFTflag=1;
%% Prepare
if exist('n')~=1||n<=0; n=1; end
if exist('dt')~=1||dt<=0; dt=1; end
if exist('FFTflag')~=1||FFTflag<=0; FFTflag=0; end

NSig=length(x);
Sur=zeros(NSig,n);

if mod(length(x),2)==1
    x=x(1:end-1);
    NSig=NSig-1;
end
Nfft=NSig/2;
x=x(:);

%% Generate surogates
MN=mean(x);
X_FFT = fft(x);
X_FFT = abs(X_FFT(2:Nfft+1));
Phi_rand = exp(1j.*2*pi*rand(Nfft,n));
B=zeros(2*Nfft,n);
B(1:Nfft,:)=bsxfun(@times,X_FFT,Phi_rand);
B=2.*ifft(B);
Sur(1:NSig,:)=real(bsxfun(@times,B,exp(1j.*(1:1:NSig)'.*2.*pi/NSig)));

%% FFTS
if FFTflag
Fs=1/(dt.*NSig);
f=Fs.*(1:1:NSig/2);

X=fft(x);
X=abs(X(1:NSig/2))./NSig;
X(1:end-1)=X(1:end-1).*2;

Y=fft(Sur(1:NSig,:),[],1);
Y=abs(Y(1:NSig/2,:))./NSig;
Y(1:end-1,:)=Y(1:end-1,:).*2;

else
f=[];
X=[];
Y=[];
end


end

