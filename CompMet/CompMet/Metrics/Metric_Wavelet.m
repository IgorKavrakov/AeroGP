function [OutStruct] = Metric_Wavelet(X1,X2,dt,WavletProperties,t,Plots)
% Calcuate Wavelet Metric 
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
if isfield(WavletProperties,'BetaCone')&&WavletProperties.BetaCone>0; BetaCone=WavletProperties.BetaCone; else BetaCone=1; end

OutStruct.Name='Wavlet Metric (W); Local frequency normalization (WI); Wavelet Mean Spectrum (WM)';    
[OutStruct.WT1,OutStruct.FreqBins1,OutStruct.Scales1] = MorletWavletTransform(X1, WavletProperties.nLevel, WavletProperties.fmax, WavletProperties.fmin, WavletProperties.f0, 1/dt,WavletProperties.beta,WavletProperties.Padding);
[OutStruct.WT2,OutStruct.FreqBins2,OutStruct.Scales2] = MorletWavletTransform(X2, WavletProperties.nLevel, WavletProperties.fmax, WavletProperties.fmin, WavletProperties.f0, 1/dt,WavletProperties.beta,WavletProperties.Padding);

%Computing cone of influence
OutStruct.ConeInfluence=zeros(size(OutStruct.WT1));
Cone=1/sqrt(2).*OutStruct.Scales1*BetaCone;

for i=1:length(OutStruct.Scales1)
    OutStruct.ConeInfluence(i,find((Cone(i)<t).*(t<t(end)-Cone(i))>0))=1;
end
 
% delta_phi=atan(imag(OutStruct.WT1)./imag(OutStruct.WT2))./pi;

%Total
OutStruct.Value=exp(-sqrt(sum(sum((abs(OutStruct.WT1.*(OutStruct.ConeInfluence))-abs(OutStruct.WT2.*(OutStruct.ConeInfluence))).^2)))/sqrt(sum(sum(abs(OutStruct.WT1.*(OutStruct.ConeInfluence)).^2))));

% Frequency Marginal
WT1_FreqM=bsxfun(@rdivide,abs(OutStruct.WT1.*(OutStruct.ConeInfluence)),max(abs(OutStruct.WT1.*(OutStruct.ConeInfluence)),[],1));
WT2_FreqM=bsxfun(@rdivide,abs(OutStruct.WT2.*(OutStruct.ConeInfluence)),max(abs(OutStruct.WT2.*(OutStruct.ConeInfluence)),[],1));
TAverage=0;A=zeros(length(t),1);
for i=1:length(t)
Val=find(OutStruct.ConeInfluence(:,i)==1);
if any(Val)
A(i)=sqrt(sum((WT1_FreqM(Val,i)-WT2_FreqM(Val,i)).^2,1))./sqrt(sum(WT1_FreqM(Val,i).^2,1));
TAverage=TAverage+1;
end
end
OutStruct.Value(2)=exp(-1./TAverage.*sum(A));

% Time Marginal
% WT1_TimeM=bsxfun(@rdivide,abs(OutStruct.WT1),max(abs(OutStruct.WT1),[],2));
% WT2_TimeM=bsxfun(@rdivide,abs(OutStruct.WT2),max(abs(OutStruct.WT2),[],2));
% OutStruct.value_tm=exp(-mean(sqrt(sum((WT1_TimeM-WT2_TimeM).^2,2))./sqrt(sum(WT1_TimeM.^2,2))));

% Frequency AverageSpectra
OutStruct.WT1_FreqMean=mean(abs(OutStruct.WT1.*(OutStruct.ConeInfluence)),2);
OutStruct.WT2_FreqMean=mean(abs(OutStruct.WT2.*(OutStruct.ConeInfluence)),2);

OutStruct.ValueScalogram=exp(-sqrt(sum(sum(abs(abs(OutStruct.WT1.*(OutStruct.ConeInfluence)).^2-abs(OutStruct.WT2.*(OutStruct.ConeInfluence)).^2).^2)))/sqrt(sum(sum(abs(OutStruct.WT1.*(OutStruct.ConeInfluence)).^2))));

%% Plots
if ~isfield(WavletProperties,'MaxNormPlot'); MaxNormPlot=1.5; else MaxNormPlot=WavletProperties.MaxNormPlot; end
    maximumX1=max(max(abs(OutStruct.WT1)));

if isfield(Plots,'Scalogram')&&Plots.Scalogram==1; 

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth, 1.2, 0.5, 1.7, 0.6, 2, 1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(t,OutStruct.FreqBins1,abs(OutStruct.WT1)./maximumX1,1024/20, 'edgecolor','none');
    shading flat
    line(BetaCone.*1/sqrt(2).*OutStruct.Scales1,OutStruct.FreqBins1,max(max(abs(OutStruct.WT1)./maximumX1))*ones(1,length(OutStruct.Scales1)),'color','r','linewidth',2)     
    line(t(end)-BetaCone.*1/sqrt(2).*OutStruct.Scales1,OutStruct.FreqBins1,max(max(abs(OutStruct.WT1)./maximumX1))*ones(1,length(OutStruct.Scales1)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Scalogram - X1 (normalized w.r.t. X1)');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(t,OutStruct.FreqBins2,abs(OutStruct.WT2)./maximumX1,1024/20, 'edgecolor','none');
    shading flat
    line(BetaCone.*1/sqrt(2).*OutStruct.Scales2,OutStruct.FreqBins2,max(max(abs(OutStruct.WT2)./maximumX1))*ones(1,length(OutStruct.Scales2)),'color','r','linewidth',2)     
    line(t(end)-BetaCone.*1/sqrt(2).*OutStruct.Scales2,OutStruct.FreqBins2,max(max(abs(OutStruct.WT2)./maximumX1))*ones(1,length(OutStruct.Scales2)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Scalogram - X2 (normalized w.r.t. X1)');
    box on; grid on; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_Scalogram'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_Scalogram']);
end

if isfield(Plots,'ScalogramDiff')&&Plots.ScalogramDiff==1; 
    Delta_a=abs(abs(OutStruct.WT1)-abs(OutStruct.WT2))./repmat(max(abs(OutStruct.WT1),[],1),[length(OutStruct.FreqBins1) 1]);

    [positions] = subplot_pos(Plots.FigWidth/2,Plots.FigDepth,1.2,0.5,1.7,0.6,1,1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(t,OutStruct.FreqBins1,Delta_a,1024/20, 'edgecolor','none');
    shading flat
    line(BetaCone.*1/sqrt(2).*OutStruct.Scales1,OutStruct.FreqBins1,max(max(abs(OutStruct.WT1)./maximumX1))*ones(1,length(OutStruct.Scales1)),'color','r','linewidth',2)     
    line(t(end)-BetaCone.*1/sqrt(2).*OutStruct.Scales1,OutStruct.FreqBins1,max(max(abs(OutStruct.WT1)./maximumX1))*ones(1,length(OutStruct.Scales1)),'color','r','linewidth',2)
    caxis([0 MaxNormPlot]);
    h=colorbar('location','southoutside','position',[0.1286    0.0525    0.8214    0.0435]);
    xlabel('$t$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Frequency norm. diff Scalogram - abs(W(X1(t,f))-W(X2(t,f))/max_f(W(X1(t,f)))','Interpreter', 'none');
    box on; grid on; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_ScalogramDiff'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_ScalogramDiff']);
end

if isfield(Plots,'ScalogramNorm')&&Plots.ScalogramNorm==1; 

    [positions] = subplot_pos(Plots.FigWidth,Plots.FigDepth,1.2,0.5,1.7,0.6,2,1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    hold on
    contourf(t,OutStruct.FreqBins1,WT1_FreqM,1024/20, 'edgecolor','none');
    shading flat
    line(BetaCone.*1/sqrt(2).*OutStruct.Scales1,OutStruct.FreqBins1,max(max(abs(OutStruct.WT1)./maximumX1))*ones(1,length(OutStruct.Scales1)),'color','r','linewidth',2)     
    line(t(end)-BetaCone.*1/sqrt(2).*OutStruct.Scales1,OutStruct.FreqBins1,max(max(abs(OutStruct.WT1)./maximumX1))*ones(1,length(OutStruct.Scales1)),'color','r','linewidth',2)
    caxis([0 1]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Scalogram - X1 (norm. w.r.t. instantanious freq)');
    box on; grid on; hold off

    ax=axes('position',positions{2,1},'Layer','top');
    hold on
    contourf(t,OutStruct.FreqBins2,WT2_FreqM,1024/20, 'edgecolor','none');
    shading flat
    line(1/sqrt(2).*OutStruct.Scales2,OutStruct.FreqBins2,max(max(abs(OutStruct.WT2)./maximumX1))*ones(1,length(OutStruct.Scales2)),'color','r','linewidth',2)     
    line(t(end)-1/sqrt(2).*OutStruct.Scales2,OutStruct.FreqBins2,max(max(abs(OutStruct.WT2)./maximumX1))*ones(1,length(OutStruct.Scales2)),'color','r','linewidth',2)
    caxis([0 1]);
    h=colorbar('location','southoutside','position',[0.0816    0.0525    0.8214    0.0435]);
    xlabel('$t$ ','Interpreter', 'latex','FontSize',9);
    ylabel('$f$ ','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    title('Scalogram - X2 (norm. w.r.t. instantanious freq)');
    box on; grid on; hold off

    print(gcf,[Plots.Name '\'  Plots.Name '_ScalogramFreqMarginal'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_ScalogramFreqMarginal']);
end

if isfield(Plots,'MeanPSD')&&Plots.MeanPSD==1; 
    Blue=[0 155/255 180/255];Red=[183/255 26/255 73/255];Green=[0.4660 0.6740 0.1880];Purp=[0.4940;0.1840;0.5560];Orange=[0.8500 0.3250  0.0980]; Gray=[0.5 0.5 0.5];

    [positions] = subplot_pos(Plots.FigWidth/2,Plots.FigDepth,1.2,0.5,1,0.6,1,1,1.5,1.5);
    f=figure(); clf(f); set(gcf, 'PaperUnits', 'centimeters'); set(gcf, 'PaperSize', [Plots.FigWidth Plots.FigDepth]); set(gcf, 'PaperPositionMode', 'manual'); set(gcf, 'PaperPosition', [0 0 Plots.FigWidth Plots.FigDepth]);

    ax=axes('position',positions{1,1},'Layer','top');
    loglog(OutStruct.FreqBins1,OutStruct.WT1_FreqMean,'color',Blue);
    hold on
    loglog(OutStruct.FreqBins2,OutStruct.WT2_FreqMean,'color',Red);
    xlabel('$f$ ','Interpreter', 'latex','FontSize',9);
    ylabel('PSD','Interpreter', 'latex','FontSize',9);
    set(gca,'FontSize',8,'FontName','Times');
    l=legend({'$X_1$','$X_2$'},'location','south');  set(l, 'Interpreter', 'latex','FontSize',8)
    title('PSD Based on Wavlet Analysis','Interpreter', 'none');
    box on; grid on; hold off
    
    print(gcf,[Plots.Name '\'  Plots.Name '_PSD_Wavelet'],'-dpdf');%savefig(gcf,[Plots.Name '\'  Plots.Name '_PSD_Wavelet']);
end


end

