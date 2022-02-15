function [Kern] = Kern_Exp(x1,x2,hyp)
% Get Exponential Kernel

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
n=size(x1,1); % Number observation points;
p=size(x2,1); %Number of prediction points;

Kern=zeros(n,p);
for i=1:D
Kern=Kern+exp(-hyp(i+1))*(bsxfun(@minus,x1(:,i),x2(:,i)')).^2;
end
Kern=exp(hyp(1)).*exp(-1/2*Kern);

end

