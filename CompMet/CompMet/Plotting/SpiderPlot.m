function SpiderPlot(Metric,Names,Color,LegendTitles)
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
TextLocRho=ones(1,length(Metric))*1.087;
TextLocTheta=deg2rad(ones(1,length(Metric))*-4);
LineType={'.-','s-','o-','v-','^-','+-','-','-','-','-','-','-'};  % Cell array of line type.
MarkerSize={8,3,3,3,3,3,3,3,3,3,3};         % MarkerSize

if nargin==3; LegendTitles=[]; end

NMetrics=max(size(Names));
[~,Col]=size(Metric); if NMetrics==Col; Metric=Metric'; end;
[~,NSets]=size(Metric);
Theta=(0:1/NMetrics:1-1/NMetrics).*2.*pi+pi/2;
Radius=[0:1/5:1];
% Generate iso lines
clear x_axes y_axes
for i=1:NMetrics
[x_axes(i,:), y_axes(i,:)] = pol2cart(Theta(i).*ones(1,length(Radius)), Radius);
end
x_axes(end+1,:)=x_axes(1,:);
y_axes(end+1,:)=y_axes(1,:);

grey = [1, 1, 1] * 0.25;
h = line(x_axes, y_axes,...
    'LineWidth', 0.2,...
    'Color', grey);

% Radial lines
RadialLines=[0 1.1];
for i=1:NMetrics
[RadialLines_X(i,:), RadialLines_Y(i,:)] = pol2cart(Theta(i).*ones(1,length(RadialLines)), RadialLines);
end
hold on
h = line(RadialLines_X', RadialLines_Y','LineWidth', 0.2,'Color', grey);


%Write axis
for i=1:length(Radius)
text(Radius(2)/10,Radius(i)+Radius(2)/3.5,num2str(Radius(i),'%.1f'),  'FontSize',8,'FontName','Times')
end

if length(TextLocRho)==1
TextLocRho=TextLocRho.*ones(NMetrics,1);    
end
if length(TextLocTheta)==1
TextLocTheta=TextLocTheta.*ones(NMetrics,1);    
end

%Write labels
for i=1:NMetrics
    if Theta(i)>=pi/2&&Theta(i)<3/2*pi; TextLocRho(i)=TextLocRho(i)*1.165; end
[x_label, y_label] = pol2cart(Theta(i)+TextLocTheta(i),TextLocRho(i));
text(x_label,y_label,Names{i}, 'Interpreter', 'LaTeX',  'FontSize',9)
end
clear x_axes y_axes

%Generate Metrics
for j=1:NSets
    for i=1:NMetrics
    [x_axes(i), y_axes(i)] = pol2cart(Theta(i),[Metric(i,j)]);
       
    end
    x_axes(end+1)=x_axes(1);y_axes(end+1)=y_axes(1);
    
        for i=1:NMetrics
            if i==1
        legendhandle(j) = plot([x_axes(i) x_axes(i+1)], [y_axes(i) y_axes(i+1)],LineType{j},'Markersize',MarkerSize{j},'Color', Color{j});
        else
        plot([x_axes(i) x_axes(i+1)], [y_axes(i) y_axes(i+1)],LineType{j},'Markersize',MarkerSize{j},'Color', Color{j});        
        end
        end
    clear x_axes y_axes
end

if ~isempty(LegendTitles)
l=legend(legendhandle,LegendTitles,'location','eastoutside');
set(l, 'FontSize',7,'interpreter','latex')
end

xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis square
axis off
set(gcf,'color','w');

end

