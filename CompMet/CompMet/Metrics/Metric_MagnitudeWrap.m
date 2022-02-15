function [ OutStruct ] = Metric_MagnitudeWrap(X1,X2,t,Plots)
% Calcuate Magnitude Metric based on Dynamic Time Warping
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


OutStruct.Name='Magnitude Metric - Wraped (MW)';
OutStruct.Value=exp(-norm(X1-X2)./norm(X1));



%% Plots
if isfield(Plots,'WarpedSignals')&&Plots.WarpedSignals==1
Blue=[0 155/255 180/255];Red=[183/255 26/255 73/255];Green=[0.4660 0.6740 0.1880];Purp=[0.4940;0.1840;0.5560];Orange=[0.8500 0.3250  0.0980]; Gray=[0.5 0.5 0.5];
dt_Wraped=t(end)/(length(X1));
t_Wraped=0:dt_Wraped:t(end)-dt_Wraped;

[positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth,1.2,0.5,1,0.6,1,1,1,1);
f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);
ax=axes('position',positions{1,1},'Layer','top');
hold on
plot(t_Wraped,X1,'color',Blue);
plot(t_Wraped,X2,'color',Red);

l=legend({'$X_1$-Warp','$X_2$-Warp'},'location','south');  set(l, 'Interpreter', 'latex','FontSize',8)
xlabel('$t$','Interpreter','latex','FontSize',9);ylabel('Amplitude','Interpreter', 'latex','FontSize',9);
set(gca,'FontSize',8,'FontName','Times');
title('Warped Signals');
box on; grid on;
print(gcf,[Plots.Name '\'  Plots.Name '_Signals_Warp'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_Signals_Wrap']);
end



end

