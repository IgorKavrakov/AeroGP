function Example1b_Plots(GP,FP)
% This is a speciic function for plots in the Flat Plate Example

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
%% Plot properties
close all;
Colors=[57  106 177
        218 124 48
        62  150 81
        204 37  41
        83  81  84
        107 76  154
        146 36  40
        148 139 61]./255;


FigWidth=18.48;nby=2; nbx=2;spacey=1; spacex=2.18; leftmargin=1.2; rightmargin=0.5; topmargin=0.2; bottommargin=0.9;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
FontAxis=7; FontLabel=9; FontLegend=7;
set(groot,'defaultLineLineWidth',0.6)

%% Plot 1  Flutter
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

FP(1).Pred.tau=(0:FP(1).Pred.Samp)*FP(1).Par.ds;
ax=axes('position',positions{1,2},'Layer','top');hold on;
h3=plot(FP(1).Pred.tau,atan(FP(3).Pred.u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(3,:));
h2=plot(FP(1).Pred.tau,atan(FP(2).Pred.u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(4,:));
h1=plot(FP(1).Pred.tau,atan(FP(1).Pred.u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_h$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ylim([-0.2 0.2]);yticks([-0.2:0.1:0.2]); ytickformat('%.1f');%xtickformat('%.1f');
l=legend([h1,h2,h3],{'Analytical Damped $V_r$=13.30','Analytical Critical $V_{r,cr}$=13.33','Analytical Divergent $V_r$=13.40'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2,2},'Layer','top');hold on;
plot(FP(1).Pred.tau,FP(3).Pred.u(:,2)*180/pi,'-','color',Colors(3,:));
plot(FP(1).Pred.tau,FP(2).Pred.u(:,2)*180/pi,'-','color',Colors(4,:));
plot(FP(1).Pred.tau,FP(1).Pred.u(:,2)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');%xtickformat('%.1f');
% l=legend({'Analytical Damped $U$=77 m/s','Analytical critical $U_{cr}$=78 m/s','Analytical divergent $U$=79 m/s'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{1,1},'Layer','top');hold on;
h3=plot(FP(1).Pred.tau,atan(GP(4).u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(3,:));
h2=plot(FP(1).Pred.tau,atan(GP(3).u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(4,:));
h1=plot(FP(1).Pred.tau,atan(GP(2).u_d(:,1)/FP(1).Par.B)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_h$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ylim([-0.2 0.2]);yticks([-0.2:0.1:0.2]);ytickformat('%.1f');%xtickformat('%.1f');
l=legend([h1,h2,h3],{'GP Damped $V_r$=13.33','GP Critical $V_{r,cr}$=13.40','GP Divergent $V_r$=13.45'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2,1},'Layer','top');hold on;
plot(FP(1).Pred.tau,GP(4).u(:,2)*180/pi,'-','color',Colors(3,:));
plot(FP(1).Pred.tau,GP(3).u(:,2)*180/pi,'-','color',Colors(4,:));
plot(FP(1).Pred.tau,GP(2).u(:,2)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');%xtickformat('%.1f');
% l=legend({'Analytical Damped $U$=77 m/s','Analytical critical $U_{cr}$=78 m/s','Analytical divergent $U$=79 m/s'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

print(gcf,'Example1_FlatPlateAnalytical/Ex1b_Flut','-dpdf')


end

%Plot function
function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey)
    subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
    subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
    for i=1:nbx
       for j=1:nby
           xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
           yfirst=bottommargin+(j-1.0)*(subysize+spacey);
           positions{i,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
       end
    end
end