function [lik,Der] = GP_Process_Opt(x,y,hyp,jitter,noiseFlag)
%Obtain the likelihood for the GP

% By Igor Kavrakov

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of AeroGP.
%  AeroGP is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  AeroGP is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with AeroGP.  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Allan McRobie, Guido Morgenthal 2022

%%
if noiseFlag~=-1
    noise=log(0);
else
    noise=hyp(end);    
end

Kern=Kern_Exp(x,x,hyp);   %Exponential kernel
Noise=eye(size(x,1)).*exp(noise)+eye(size(x,1)).*jitter; %Noise+jitter term

[L_Kern,p]=chol(Kern+Noise,'lower');
if p>0; lik=inf; Der=hyp*Inf; return; else; end
alpha=L_Kern'\(L_Kern\y);
lik=1/2*y'*alpha+sum(log(diag(L_Kern)))+length(y)/2*log(2*pi);

if nargout==2
Q  =  L_Kern'\(L_Kern\eye(size(Kern))) - alpha*alpha';
Der = Kern_Exp_Der(x,x,hyp,Kern,Q);        %Get the derivative of the kernels
    if noiseFlag==-1
    Der(end+1)=exp(noise)*trace(Q)/2;%Get the derivative of the noise
    end
end

end

