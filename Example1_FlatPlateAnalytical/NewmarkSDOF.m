function [u1,udot1,u2dot1] = NewmarkSDOF(p,u,udot,u2dot,beta,gamma,dt,m,c,k)
% Newmark-beta Method for SDOF Linear Systems 
% In book P177 Structural dynamics - Chopra

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

%% Implementation
% p: Force for i time step;
% u: Displacement from i time step;
% udot: Velocity from i time step;
% u2dot: Acceleration from i time step;
% beta,gamma: Time integration parameters;
% dt: Time step;
% m:  Mass;
% c:  Damping;
% k:  Stifness;

 %1.3
a1=1/(beta.*dt^2).*m+gamma./(beta.*dt).*c;
a2=1/(beta.*dt).*m+(gamma./beta-1).*c;
a3=(1/(2.*beta)-1).*m+dt.*(gamma./(2.*beta)-1).*c;
 %1.4
kbar=k+a1;
%2.0 Calculation for for next time step
%2.1
pbar=p+a1.*u+a2.*udot+a3.*u2dot;
%2.2
u1=pbar./kbar;
%2.3
udot1=gamma./(beta.*dt).*(u1-u)+(1-gamma./beta).*udot+dt.*(1-gamma./(2.*beta)).*u2dot;
%2.4      
u2dot1=1./(beta.*dt.^2)*(u1-u)-1./(beta.*dt).*udot-(1./(2.*beta)-1)*u2dot;

end