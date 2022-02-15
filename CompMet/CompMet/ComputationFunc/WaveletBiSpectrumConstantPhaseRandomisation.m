function [f1,WavletAutoBiSpectrum] = WaveletBiSpectrumConstantPhaseRandomisation(WT,FreqBins,dt,T,NIntervals,Flag,Fact,eps)
% Calcuate wavelet bi-spectrum using Phase Randomisation 
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

%%
% WT - Wavlet coefficients;
% f - Frequencies
% t -  Time vector
% T - time for averaging [TStart TEnd]
% NIntervals - Number of intervals to be splitted
% Rand - Randomization variable for the phase (0 if we don't want);
% Fact - for phase randomization of the wavelet bicoherence
if nargin==3; T(1)=0; T(2)=dt*size(WT,2);NIntervals=1; Flag=1; Fact=10*pi; eps=0.2; end
if nargin==4; NIntervals=1; Flag=1;Fact=10*pi;eps=0.2; end
if nargin==5; Flag=1;Fact=10*pi;eps=0.2; end
if nargin==6; Fact=10*pi; eps=0.2;end
if nargin==7; eps=0.2;end

if size(T)==1; T(1)=0; T(2)=dt*size(WT,2); end


TStart=T(1);TEnd=T(2);
StepStart=max(1,floor(TStart/dt));StepEnd=min(floor(TEnd/dt),size(WT,2));
NSteps=StepEnd-StepStart+1;
WT=WT(:,StepStart:1:StepEnd);
IntervalSteps=floor(NSteps/NIntervals);


BinsConsidered=find(FreqBins<FreqBins(end)/2);
BinsConsidered=BinsConsidered(end);
WavletAutoBiSpectrum=zeros(BinsConsidered,BinsConsidered,NIntervals);
if Flag==1; NIntervals=1; end

Mplyer=1;
f1=FreqBins(1:BinsConsidered);
indx=hankel([1:BinsConsidered],[BinsConsidered,1:BinsConsidered-1])+1;
indx=reshape(indx,[BinsConsidered.^2,1]);
ConstRandom = (pi-eps).*rand(BinsConsidered,BinsConsidered) + eps;
Sec=randi([0 1],BinsConsidered,BinsConsidered);Sec(Sec==0)=-1;
ConstRandom=ConstRandom.*Sec; clear Sec;

for j=1:NIntervals

   for k=1:IntervalSteps-1
            CC=reshape(WT(indx,k+(j-1)*IntervalSteps),[BinsConsidered BinsConsidered]);
            CC2=WT(1:BinsConsidered,k+(j-1)*IntervalSteps)*(WT(1:BinsConsidered,k+(j-1)*IntervalSteps).');
            Random=(rand(size(ConstRandom))-0.5).*Fact*2;
            CohOneStep=abs(CC2.*conj(CC)).*exp(1j.*Random.*ConstRandom);
            WavletAutoBiSpectrum(:,:,j)=WavletAutoBiSpectrum(:,:,j)+Mplyer.*CohOneStep.*dt./2;
            Mplyer=2;
    end

        Mplyer=1;
        CC=reshape(WT(indx,IntervalSteps+(j-1)*IntervalSteps),[BinsConsidered BinsConsidered]);
        CC2=WT(1:BinsConsidered,IntervalSteps+(j-1)*IntervalSteps)*(WT(1:BinsConsidered,IntervalSteps+(j-1)*IntervalSteps).');
        Random=(rand(size(ConstRandom))-0.5).*Fact*2;
        CohOneStep=abs(CC2.*conj(CC)).*exp(1j.*Random.*ConstRandom);
        WavletAutoBiSpectrum(:,:,j)=WavletAutoBiSpectrum(:,:,j)+Mplyer.*CohOneStep.*dt./2;
        a=ones(BinsConsidered,BinsConsidered);
        WavletAutoBiSpectrum(:,:,j)=WavletAutoBiSpectrum(:,:,j).*triu(a).*flip(tril(a));
 
    
end


