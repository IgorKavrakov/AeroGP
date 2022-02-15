function Example1a_Plots(GP,FP)
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

%% Plot 1 - Training set input
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

GP.Train(1).x=GP.Train(1).x*GP.Train(1).x_max;% Denormalize
GP.Train(2).x=GP.Train(2).x*GP.Train(2).x_max;% Denormalize

ax=axes('position',positions{3},'Layer','top');hold on;
h10=fill([50,100,100,50]',[-0.4,-0.4,0.4,0.4]',[0.8 0.8 0.8],'LineStyle','none','facealpha',.5);
plot(FP.Train.tau,GP.Train(1).x(:,1)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_h$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');
xlim([0 280]);xticks(0:40:280);
grid on; box on;

ax=axes('position',positions{4},'Layer','top');hold on;
h10=fill([50,100,100,50]',[-0.4,-0.4,0.4,0.4]',[0.8 0.8 0.8],'LineStyle','none','facealpha',.5);
plot(FP.Train.tau,GP.Train(1).x(:,2+GP.Par.Lag)*180/pi,'-','color',Colors(4,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha_a$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');
xlim([0 280]);xticks(0:40:280);
grid on; box on;

GP.Pred(1).y=GP.Pred(1).y*GP.Train(1).y_max;% Denormalize
GP.Pred(2).y=GP.Pred(2).y*GP.Train(2).y_max;% Denormalize
GP.Train(1).y=GP.Train(1).y*GP.Train(1).y_max;% Denormalize
GP.Train(2).y=GP.Train(2).y*GP.Train(2).y_max;% Denormalize
GP.Pred(1).y_targ=GP.Pred(1).y_targ*GP.Train(1).y_max;% Denormalize
GP.Pred(2).y_targ=GP.Pred(2).y_targ*GP.Train(2).y_max;% Denormalize
GP.Pred(1).y_std=sqrt(diag(GP.Pred(1).y_var))*GP.Train(1).y_max; %Denormalize
GP.Pred(2).y_std=sqrt(diag(GP.Pred(2).y_var))*GP.Train(2).y_max; %Denormalize

ax=axes('position',positions{1},'Layer','top');hold on;
h10=fill([FP.Train.tau; flip(FP.Train.tau)], [GP.Pred(1).y+2.56.*GP.Pred(1).y_std; flip(GP.Pred(1).y-2.56*GP.Pred(1).y_std)],Colors(1,:),'LineStyle','none','facealpha',.3);
h1=plot(FP.Train.tau,GP.Pred(1).y_targ,'-','color',Colors(4,:));
h3=plot(FP.Train.tau,GP.Pred(1).y,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_L$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.2f');%xtickformat('%.1f');
l=legend([h1 h3],{'Analytical','GP Train'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([50 100]);
grid on; box on;

ax=axes('position',positions{2},'Layer','top');hold on;
h10=fill([FP.Train.tau; flip(FP.Train.tau)], [GP.Pred(2).y+2.56.*GP.Pred(2).y_std; flip(GP.Pred(2).y-2.56*GP.Pred(2).y_std)],Colors(1,:),'LineStyle','none','facealpha',.3);
h1=plot(FP.Train.tau,GP.Pred(2).y_targ,'-','color',Colors(4,:));
h3=plot(FP.Train.tau,GP.Pred(2).y,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.3f');%xtickformat('%.1f');
l=legend([h1 h3],{'Analytical','GP Train'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([50 100])
grid on; box on;
f.Renderer='Painters';
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_Train','-dpdf')

%% Plot 2 - Random prediction output
nby=1;FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

GP.Pred(3).x=GP.Pred(3).x*GP.Train(3).x_max;% Denormalize
GP.Pred(3).y=GP.Pred(3).y*GP.Train(3).y_max;% Denormalize
GP.Pred(3).y_std=sqrt(diag(GP.Pred(3).y_var))*GP.Train(3).y_max; %Denormalize
GP.Pred(3).y_targ=GP.Pred(3).y_targ*GP.Train(3).y_max;% Denormalize
GP.Pred(3).tau=(0:length(GP.Pred(3).y)-1)'*FP.Par.ds;

ax=axes('position',positions{1,1},'Layer','top');hold on;
plot(GP.Pred(3).tau,GP.Pred(3).x(:,1)*180/pi,'-','color',Colors(4,:));
plot(GP.Pred(3).tau,GP.Pred(3).x(:,2+GP.Par.Lag)*180/pi,'-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$\alpha$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');%xtickformat('%.1f');
l=legend({'$\alpha_h$','$\alpha_a$'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2},'Layer','top');hold on;
h10=fill([GP.Pred(3).tau; flip(GP.Pred(3).tau)], [GP.Pred(3).y+2.56.*GP.Pred(3).y_std; flip(GP.Pred(3).y-2.56*GP.Pred(3).y_std)],Colors(1,:),'LineStyle','none','facealpha',.3);
h1=plot(GP.Pred(3).tau,GP.Pred(3).y_targ,'-','color',Colors(4,:));
h3=plot(GP.Pred(3).tau,GP.Pred(3).y,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis);
xlabel('$\tau$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.3f');%xtickformat('%.1f');
l=legend([h1 h3],{'Analytical','GP Prediction'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_RandOut','-dpdf')

%% Plot 3 - Harmonic prediction output: Lift and Moment at Vr=6 due to alpha_a
topmargin=0.4; FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

StepsVr=cumsum([1 FP.Pred.Cycl.*FP.Pred.Vr./FP.Par.ds+1]);StepsVr(end)=StepsVr(end);
Vr=6; %To predict this;
VrInd=find(FP.Pred.Vr==Vr);
VrInd=[StepsVr(VrInd),StepsVr(VrInd+1)-1];

GP.Pred(6).x=GP.Pred(6).x*GP.Train(6).x_max;% Denormalize
GP.Pred(6).y=GP.Pred(6).y*GP.Train(6).y_max;% Denormalize
GP.Pred(6).y_targ=GP.Pred(6).y_targ*GP.Train(6).y_max;% Denormalize
GP.Pred(6).tau=(0:length(GP.Pred(6).y)-1)'*FP.Par.ds;

GP.Pred(7).x=GP.Pred(7).x*GP.Train(7).x_max;% Denormalize
GP.Pred(7).y=GP.Pred(7).y*GP.Train(7).y_max;% Denormalize
GP.Pred(7).y_targ=GP.Pred(7).y_targ*GP.Train(7).y_max;% Denormalize
tau=(0:VrInd(2)-VrInd(1))'*FP.Par.ds./Vr;

ax=axes('position',positions{1},'Layer','top');hold on;
h1=plot(tau,GP.Pred(6).y_targ(VrInd(1):VrInd(2))*5,'-','color',Colors(4,:));
h3=plot(tau,GP.Pred(6).y(VrInd(1):VrInd(2))*5,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis); hAxes=gca;hAxes.YAxis.Exponent = 0;
xlabel('$\tau/V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_L$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.2f');
l=legend([h1 h3],{'Analytical','GP Prediction'},'location','best');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
grid on; box on;

ax=axes('position',positions{2},'Layer','top');hold on;
h1=plot(tau,GP.Pred(7).y_targ(VrInd(1):VrInd(2))*5,'-','color',Colors(4,:));
h3=plot(tau,GP.Pred(7).y(VrInd(1):VrInd(2))*5,'--','color',Colors(1,:),'LineWidth',0.8);
set(gca, 'FontSize',FontAxis); hAxes=gca;hAxes.YAxis.Exponent = 0;
xlabel('$\tau/V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$C_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.3f');
grid on; box on;
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_HarmonicOut','-dpdf')

%% Plot 4 - Flutter Derivatives
%Preprocess Lift and moment due to alpha_h
GP.Pred(4).x=GP.Pred(4).x*GP.Train(4).x_max;% Denormalize
GP.Pred(4).y=GP.Pred(4).y*GP.Train(4).y_max;% Denormalize
GP.Pred(4).y_targ=GP.Pred(4).y_targ*GP.Train(4).y_max;% Denormalize

GP.Pred(5).x=GP.Pred(5).x*GP.Train(5).x_max;% Denormalize
GP.Pred(5).y=GP.Pred(5).y*GP.Train(5).y_max;% Denormalize
GP.Pred(5).y_targ=GP.Pred(5).y_targ*GP.Train(5).y_max;% Denormalize

for i=1:length(FP.Pred.Vr)
 S=StepsVr(i);
 E=StepsVr(i+1)-1;
    [H1(i),H4(i)] = DerFit(FP.Pred.Vr(i),atan(GP.Pred(4).x(S:E,1)),atan(GP.Pred(4).x(S:E,3+2*GP.Par.Lag)),GP.Pred(4).y(S:E),1);
    [A1(i),A4(i)] = DerFit(FP.Pred.Vr(i),atan(GP.Pred(5).x(S:E,1)),atan(GP.Pred(5).x(S:E,3+2*GP.Par.Lag)),GP.Pred(5).y(S:E),1);
    [H2(i),H3(i)] = DerFit(FP.Pred.Vr(i),GP.Pred(6).x(S:E,2+GP.Par.Lag),GP.Pred(6).x(S:E,4+2*GP.Par.Lag),GP.Pred(6).y(S:E),0);
    [A2(i),A3(i)] = DerFit(FP.Pred.Vr(i),GP.Pred(7).x(S:E,2+GP.Par.Lag),GP.Pred(7).x(S:E,4+2*GP.Par.Lag),GP.Pred(7).y(S:E),0);
end
K=2*pi./FP.Pred.Vr;
C=conj(besselh(1,K/2))./(conj(besselh(1,K/2))+ 1j.*conj(besselh(0,K/2)));
G=imag(C); F=real(C);

H1_ana=-2.*pi.*F./K;              
H2_ana=-pi./2*(1+4.*G./K+F)./K;     
H3_ana=-pi.*(2.*F-0.5.*G.*K)./(K.*K); 
H4_ana=pi.*0.5.*(1+4.*G./K);        
A1_ana=0.5.*pi.*F./K;              
A2_ana=-pi./(2.*K.*K).*(K./4-G-K.*F./4); 
A3_ana=0.5.*pi.*(K.*K./32+F-K.*G./4)./(K.*K); 
A4_ana=-0.5.*pi.*G./K; 

%Plot
nby=2;
FigDepth=nby*4+topmargin+bottommargin+(nby-1)*spacey;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

ax=axes('position',positions{1,2},'Layer','top');hold on;
h1=plot(FP.Pred.Vr,H1_ana,'-','color',Colors(4,:));
h2=plot(FP.Pred.Vr,H1,'.-','color',Colors(4,:));
h3=plot(FP.Pred.Vr,H4_ana,'-','color',Colors(1,:));
h4=plot(FP.Pred.Vr,H4,'.-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$H_1^*,H_4^*$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylimit=get(gca,'ylim'); plot([14 14],ylimit,'k'); plot([2 2],ylimit,'k');
l=legend([h1 h2 h3 h4],{'$H_1^*$ Analytical','$H_1^*$ GP','$H_4^*$ Analytical','$H_4^*$ GP'},'location','southwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 16]);xticks(0:2:16);
grid on; box on;

ax=axes('position',positions{2,2},'Layer','top');hold on;
h1=plot(FP.Pred.Vr,H2_ana,'-','color',Colors(4,:));
h2=plot(FP.Pred.Vr,H2,'.-','color',Colors(4,:));
h3=plot(FP.Pred.Vr,H3_ana/2,'-','color',Colors(1,:));
h4=plot(FP.Pred.Vr,H3/2,'.-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$H_2^*,H_3^*/2$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylimit=get(gca,'ylim'); 
plot([14 14],ylimit,'k'); plot([2 2],ylimit,'k');
l=legend([h1 h2 h3 h4],{'$H_2^*$ Analytical','$H_2^*$ GP','$H_3^*$ Analytical','$H_3^*$ GP'},'location','southwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 16]);xticks(0:2:16);
grid on; box on;

ax=axes('position',positions{1,1},'Layer','top');hold on;
h1=plot(FP.Pred.Vr,A1_ana,'-','color',Colors(4,:));
h2=plot(FP.Pred.Vr,A1,'.-','color',Colors(4,:));
h3=plot(FP.Pred.Vr,A4_ana,'-','color',Colors(1,:));
h4=plot(FP.Pred.Vr,A4,'.-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$A_1^*,A_4^*$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ytickformat('%.1f');%xtickformat('%.1f');
ylimit=get(gca,'ylim'); plot([14 14],ylimit,'k'); plot([2 2],ylimit,'k');
l=legend([h1 h2 h3 h4],{'$A_1^*$ Analytical','$A_4^*$ GP','$A_1^*$ Analytical','$A_4^*$ GP'},'location','northwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 16]);xticks(0:2:16);
grid on; box on;

ax=axes('position',positions{2,1},'Layer','top');hold on;
h1=plot(FP.Pred.Vr,A2_ana,'-','color',Colors(4,:));
h2=plot(FP.Pred.Vr,A2,'.-','color',Colors(4,:));
h3=plot(FP.Pred.Vr,A3_ana/2,'-','color',Colors(1,:));
h4=plot(FP.Pred.Vr,A3/2,'.-','color',Colors(1,:));
set(gca, 'FontSize',FontAxis);
xlabel('$V_r$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$A_2^*,A_3^*/2$ [-]','Interpreter', 'latex','FontSize',FontLabel)
ylimit=get(gca,'ylim'); plot([14 14],ylimit,'k'); plot([2 2],ylimit,'k');
l=legend([h1 h2 h3 h4],{'$A_2^*$ Analytical','$A_3^*$ GP','$A_2^*$ Analytical','$A_3^*$ GP'},'location','northwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([0 16]);xticks(0:2:16);
grid on; box on;
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_FlutterDer','-dpdf')

%% Plot 5 - Comparison Metrics: Random prediction output
addpath(genpath('CompMet'));
fmin=0.01;fmax=0.5;f0=0.8;Tc=1/2/0.08;            %Metric properties values (we work in nondimensional units! *2pi)
Prop = CompMetProp(FP.Par.ds,fmin,fmax,f0,Tc);  %Metric properties 

X1=GP.Pred(3).y_targ;     %Signals 
X2=GP.Pred(3).y;          %Signals
M=CompMet(X1,X2,Prop); %Calculate (y_targ - Reference!)
Metric(1,1)=M.Phase.Value; Metric(2,1)=M.RMS.Value;  Metric(3,1)=M.MagnitudeMetricWrap.Value;  
Metric(4,1)=M.Peak.Value;  Metric(5,1)=M.PDF.Value;  Metric(6,1)=M.Wavelet.Value(1);          Metric(7,1)=M.Wavelet.Value(2);  

%Plot
FigWidth=5.5;nby=1; nbx=1;leftmargin=0.2; rightmargin=0.3; topmargin=0.25; bottommargin=0;
FigDepth=FigWidth;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

ax=axes('position',positions{1},'Layer','top');hold on;
SpidColor{1}=Colors(1,:);
SpdName={'$\mathcal{M}_{\varphi}$','$\mathcal{M}_\mathrm{rms}$','$\mathcal{M}_m$','$\mathcal{M}_p$','$\mathcal{M}_\mathrm{pdf}$','$\mathcal{M}_w$','$\mathcal{M}_{wf}$'};
SpiderPlot(Metric,SpdName,SpidColor)
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_MetricRand','-dpdf')

%% Plot 6 - Comparison Metrics: Sinusoidal prediction output
fmin=0.01;fmax=0.5;f0=6;Tc=1/2/6;                     %Metric properties values (note the f0 to be reduced frequency, a multiplication of 2pi is needed)
Prop = CompMetProp(FP.Par.ds,fmin,fmax,f0,Tc);  %Metric properties 

X1=GP.Pred(6).y_targ(VrInd(1):VrInd(2));
X2=GP.Pred(6).y(VrInd(1):VrInd(2));     
M=CompMet(X1,X2,Prop); %Calculate (y_targ - Reference!) %Lift
Metric(1,1)=M.Phase.Value; Metric(2,1)=M.RMS.Value;  Metric(3,1)=M.MagnitudeMetricWrap.Value;  
Metric(4,1)=M.Peak.Value;  Metric(5,1)=M.PDF.Value;  Metric(6,1)=M.Wavelet.Value(1);          Metric(7,1)=M.Wavelet.Value(2);  

X1=GP.Pred(7).y_targ(VrInd(1):VrInd(2));
X2=GP.Pred(7).y(VrInd(1):VrInd(2));    
M=CompMet(X1,X2,Prop); %Calculate (y_targ - Reference!) %Moment
Metric(1,2)=M.Phase.Value; Metric(2,2)=M.RMS.Value;   Metric(3,2)=M.MagnitudeMetricWrap.Value;  
Metric(4,2)=M.Peak.Value;  Metric(5,2)=M.PDF.Value;   Metric(6,2)=M.Wavelet.Value(1);         Metric(7,2)=M.Wavelet.Value(2);  

%Plot
FigWidth=7.5;nby=1; nbx=1;leftmargin=0.2; rightmargin=0.3; topmargin=0.25; bottommargin=0;
FigDepth=5.5;
f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

ax=axes('position',positions{1},'Layer','top');hold on;
SpidColor{1}=Colors(1,:);
SpidColor{2}=Colors(4,:);
LegendTitles={'$C_L$','$C_M$'};
SpdName={'$\mathcal{M}_{\varphi}$','$\mathcal{M}_\mathrm{rms}$','$\mathcal{M}_m$','$\mathcal{M}_p$','$\mathcal{M}_\mathrm{pdf}$','$\mathcal{M}_w$','$\mathcal{M}_{wf}$'};
SpiderPlot(Metric,SpdName,SpidColor,LegendTitles)
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_MetricSin','-dpdf')


%% Plot 7 - Static wind coefficients
% %Preprocess Lift and moment due to alpha
GP.Pred(8).x=GP.Pred(8).x*GP.Train(8).x_max;% Denormalize
GP.Pred(8).y=GP.Pred(8).y*GP.Train(8).y_max;% Denormalize
GP.Pred(8).y_targ=GP.Pred(8).y_targ*GP.Train(8).y_max;% Denormalize

GP.Pred(9).x=GP.Pred(9).x*GP.Train(9).x_max;% Denormalize
GP.Pred(9).y=GP.Pred(9).y*GP.Train(9).y_max;% Denormalize
GP.Pred(9).y_targ=GP.Pred(9).y_targ*GP.Train(9).y_max;% Denormalize

FP.Pred.StatTau=20;          %P Cycles - number of cycles for prediction
Inc=(FP.Pred.StatTau./FP.Par.ds+1); %Initiate

CL=zeros(length(FP.Pred.StatAlpha),1);
CM=zeros(length(FP.Pred.StatAlpha),1);
for i=1:length(FP.Pred.StatAlpha)
CL(i)=mean(GP.Pred(8).y(1+(i-1)*Inc:i*Inc))*-1;  %-4 as CM=CL/-4
CM(i)=mean(GP.Pred(9).y(1+(i-1)*Inc:i*Inc))*4;     
end
C_Analytical=FP.Pred.StatAlpha*pi/180*pi*2; %Analytical CL

%Plot
leftmargin=1.2; rightmargin=0.5; topmargin=0.1; bottommargin=0.9;
FigWidth=9; nbx=1; nby=1;FigDepth=3.5+1;

f=figure('visible','on');clf(f);set(gcf, 'PaperUnits', 'centimeters');set(gcf, 'PaperSize', [FigWidth FigDepth]);set(gcf, 'PaperPositionMode', 'manual');set(gcf, 'PaperPosition', [0 0 FigWidth FigDepth]);
[ positions ] = subplot_pos(FigWidth,FigDepth,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey);

ax=axes('position',positions{1},'Layer','top');hold on;
h1=plot(FP.Pred.StatAlpha,CL,'.-','color',Colors(1,:));
h2=plot(FP.Pred.StatAlpha,CM,'.-','color',Colors(4,:));
h3=plot(FP.Pred.StatAlpha,C_Analytical,'-','color','k');
set(gca, 'FontSize',FontAxis);
xlabel('$\alpha$ [deg]','Interpreter', 'latex','FontSize',FontLabel)
ylabel('$-C_L,4C_M$ [-]','Interpreter', 'latex','FontSize',FontLabel)
l=legend([h1 h2 h3],{'$C_L$ GP','$C_M$ GP','Analytical'},'location','northwest');
set(l,'FontSize',FontLegend,'Interpreter', 'latex');
xlim([-1 1]);xticks(-1:0.5:1);
xtickformat('%.1f');ytickformat('%.2f');

grid on; box on;
print(gcf,'Example1_FlatPlateAnalytical/Ex1a_StatWind','-dpdf')


end
%% Help functions
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

%Derivative fit function
function [Der1,Der2] = DerFit(Vr,alpha,alpha_d,F,Type)
K=2*pi/Vr;
if Type %Vertical Displacement Related Derivatives
    Ss(:,1)=K.*alpha;
    Ss(:,2)=-alpha_d;
else    %Rotation Related Derivatives
    Ss(:,1)=K*alpha_d;
    Ss(:,2)=K.^2.*alpha;
end
StF=Ss'*F; %Perform least square fit
StS=Ss'*Ss;
Der=StS\StF;

Der1=Der(1); %Velocity related derivative
Der2=Der(2); %Displacement related derivative
end

%Comparison metrics standard values
function [Prop] = CompMetProp(dt,fmin,fmax,f0,Tc)
%For details, see: 
%Kavrakov, I., Kareem, A., and Morgenthal, G. 2020. Comparison Metrics for Time-histories: Application to Bridge Aerodynamics. J. Eng. Mech., 146 (9), 040200093.
Prop.dt=dt;             %Signal time-step

%Metric activision (Which metrics to be activated: 1- Active; 0 - Inactive)
Prop.Metrics.Phase=1;           %Phase & Cross correlation 
Prop.Metrics.Peak=1;            %Peak 
Prop.Metrics.RMS=1;             %RMS 
Prop.Metrics.MagnitudeWarp=1;   %Magnitude Warped 
Prop.Metrics.Magnitude=0;       %Magnitude Unwarped 
Prop.Metrics.Wavelet=1;         %Wavelet 
Prop.Metrics.PDF=1;             %Probability Distribution metric 
Prop.Metrics.Stationarity=0;    %Stationarity 
Prop.Metrics.WaveletBicoherence=0; %Bicoherence (nonlinear coupling)

%Phase metric property
Prop.TPhase=Tc;       %Considered significant delay (cf. Eq. 2)

%Wavelet metric properties
Prop.WavletProperties.nLevel=1200;   %Frequency scales
Prop.WavletProperties.beta=3;         %cf. Kijewski & Kareem - beta factor for wavelet (keeping 3 is fine)  
Prop.WavletProperties.fmax=fmax;       %(* Maximum frrequency (Strongly dependent on the signals compared)
Prop.WavletProperties.fmin=fmin;       %(* Minimum frequency  (Strongly dependent on the signals compared)
Prop.WavletProperties.f0=f0;           %(* Wavelet central frequency (cf. Eq. 10) (Strongly dependent on the signals compared)
Prop.WavletProperties.Padding=0;       %(Recommended 0) Should the signal be padded for the end effects? Generally no

%Probability Distribution Metric Properties. Using Kernel method according to Botev (cf. Botev, Z., Grotowski, J., and Kroese, D. (2010). Kernel density estimation via diffusion. Ann.763 Stat., 38, 2916–2957)
Prop.PDFProperties.KernelType='Adaptive';  % Adaptive kernel fit (Cf. Botev)
Prop.PDFProperties.Kerneldiscretization=2^16; %Discretization of the kernel result (should be power of 2) this is discretization of the integral
Prop.PDFProperties.StandardScore=1;   %Comparison of the Standard Score of the PDFs (i.e. zero mean with normalized standard deviation). Otherwise =0 (full PDF). 

%No plots
Prop.Plots.General=0;               %Plot time histories & Bar Metrics

end

