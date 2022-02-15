function [Der] = Kern_Exp_Der(x1,x2,hyp,Kern,Q)
% Get Exponential Kernel Derivatives

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
D=size(x1,2); % Dimensionality;
Der=zeros(1+D,1); %Derivatives

%Derivative w.r.t. A
Der(1) = trace(Q*Kern)/2; 

%Derivative w.r.t. length scale
for i=1:D
Kern_D=1/2*exp(-hyp(i+1)).*(bsxfun(@minus,x1(:,i),x2(:,i)')).^2.*Kern;
Der(1+i)=trace(Q*Kern_D)/2; 
end

end

