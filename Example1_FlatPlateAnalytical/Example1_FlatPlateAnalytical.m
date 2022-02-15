function [FP] = Example1_FlatPlateAnalytical(FP,Lag)
% This function gives the training and prediction data based on analytical
% Flat Plate solution

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

%% Time step & Initiate prediction
ds=FP.Par.ds;

alpha_h_pr=[];
alpha_h_d_pr=[];
alpha_a_pr=[];
alpha_a_d_pr=[];
CL_pr=[];
CM_pr=[];
PPoint=[];
    

%% Prediction & Training
    if ~isfield(FP,'Pred')||strcmp(FP.Pred.Excitation,'Train') %Training and prediction set are same
    VrR=FP.Train.Vr;
    Alpha_STD=FP.Train.Alpha_STD*pi/180;
    Alpha_AmpR=FP.Train.Alpha_AmpRel;

    Cycl=FP.Train.Cycl;
    TPoint=Cycl*VrR(end)./ds;
    tau=ds*(0:1:TPoint-1).'; %Reduced time

    Vr=ds*TPoint./(1:TPoint)';
    Vr(TPoint/2:end)=-inf;
    Vals= Vr>=VrR(1) & Vr<=VrR(2);

    X_FFT=(1-Alpha_AmpR).*rand(TPoint,1) + Alpha_AmpR;
    Phi_rand = exp(1j.*2*pi*rand(TPoint,1));
    X_FFT(~Vals)=0;

    alpha_h=real(ifft(X_FFT.*Phi_rand));
    alpha_h_d=real(ifft(1j*2*pi./Vr.*X_FFT.*Phi_rand));

    X_FFT=(1-Alpha_AmpR).*rand(TPoint,1) + Alpha_AmpR;
    Phi_rand = exp(1j.*2*pi*rand(TPoint,1));
    X_FFT(~Vals)=0;

    alpha_a=real(ifft(X_FFT.*Phi_rand));
    alpha_a_d=real(ifft(1j*2*pi./Vr.*X_FFT.*Phi_rand));
    alpha_a_2d=real(ifft((1j*2*pi./Vr).^2.*X_FFT.*Phi_rand));

    FactAmp_h=Alpha_STD./std(alpha_h);
    FactAmp_a=Alpha_STD./std(alpha_a);

    alpha_h=(alpha_h-mean(alpha_h)).*FactAmp_h; 
    alpha_h_d=(alpha_h_d-mean(alpha_h_d)).*FactAmp_h; 

    alpha_a=(alpha_a-mean(alpha_a)).*FactAmp_a; 
    alpha_a_d=(alpha_a_d-mean(alpha_a_d)).*FactAmp_a; 
    alpha_a_2d=(alpha_a_2d-mean(alpha_a_2d)).*FactAmp_a; 

    Phi_se=1-0.165.*exp(-0.089.*tau)-0.335.*exp(-0.6.*tau);% Wagner function        
    h_primeB=sec(alpha_h).^2.*alpha_h_d;     %H_dd/B. Relationship 1+tan^2(alpha)=sec^2(alpha) is used.
    h_conv=ConvFFT(Phi_se,h_primeB,ds);      %Rise due to vert vel
    a_conv=ConvFFT(Phi_se,alpha_a_d,ds);      %Rise due to rotation 
    a_d_conv=ConvFFT(Phi_se,alpha_a_2d,ds/4); %Rise due to angular velocity

    CLh=2*pi*(-h_conv-1/4*h_primeB);                                      %Lift due to h
    CMh=pi/2*( h_conv);                                                       %Moment due to h
    CLa=2*pi*(-a_conv-a_d_conv-1/4*alpha_a_d);                   %Lift due to a
    CMa=pi/2*( a_conv+a_d_conv-1/4*alpha_a_d-alpha_a_2d/32);  %Moment due to a

    CL=CLh-mean(CLh)+CLa-mean(CLa);
    CM=CMh-mean(CMh)+CMa-mean(CMa);

    alpha_h(:,2:Lag)=0; 
    alpha_a(:,2:Lag)=0; 

    for ll=1:Lag %Regressors
        alpha_h(1+ll:end,1+ll)=alpha_h(1:TPoint-ll,1);   
        alpha_a(1+ll:end,1+ll)=alpha_a(1:TPoint-ll,1); 
    end
        
    alpha_h_pr=alpha_h;
    alpha_h_d_pr=alpha_h_d;
    alpha_a_pr=alpha_a;
    alpha_a_d_pr=alpha_a_d;
    CL_pr=CL;
    CM_pr=CM;
    PPoint=TPoint;
    
    FP.Train.alpha_h=alpha_h;
    FP.Train.alpha_h_d=alpha_h_d;
    FP.Train.alpha_a=alpha_a;
    FP.Train.alpha_a_d=alpha_a_d;
    FP.Train.CL=CL;
    FP.Train.CM=CM;
    FP.Train.Samp=TPoint;
    FP.Train.tau=tau;
       
    elseif strcmp(FP.Pred.Excitation,'Rand')%Random excitation
        load('FP_RandInp.mat','a_pr','a_d_pr','a_2d_pr','h_pr','h_d_pr','ds'); %Loads the input for the prediction ,alpha_pr, alpha_d_pr, and alpha_2d_pr, ds(Reduced time step). These are used for prediction    
        PPoint=length(a_pr); 
        tau_predict=(0:ds:ds*(PPoint-1))'; %Time & Reduced time
        Phi_se=1-0.165.*exp(-0.089.*tau_predict)-0.335.*exp(-0.6.*tau_predict);% Wagner function
        h_conv=ConvFFT(Phi_se,h_d_pr,ds);      %Rise due to vert vel
        a_conv=ConvFFT(Phi_se,a_d_pr,ds);  %Rise lift due to rotation 
        a_d_conv=ConvFFT(Phi_se,a_2d_pr,ds/4); %Rise lift due to angular velocity

        alpha_h_pr=zeros(PPoint,1+Lag); 
        alpha_h_d_pr=1./(1+h_pr.^2).*h_d_pr;
        alpha_a_pr=zeros(PPoint,1+Lag); 
        alpha_a_d_pr=a_d_pr;
        
        CMh_pr=pi/2*( h_conv);                                                     
        CMa_pr=pi/2*( a_conv+a_d_conv-1/4*a_d_pr-a_2d_pr/32);
              
        CL_pr=zeros(PPoint,1);
        CM_pr=CMh_pr+CMa_pr;

        for ll=0:Lag% Add the lag, leaving out 0 initial conditions for the regressors
            alpha_h_pr(1+ll:end,1+ll)  =atan(h_pr(1:end-ll,1));
            alpha_a_pr(1+ll:end,1+ll)  =a_pr(1:end-ll,1);              
        end
    elseif strcmp(FP.Pred.Excitation,'SinH')|| strcmp(FP.Pred.Excitation,'SinA')  %Sinusoidal excitation
        Vr=FP.Pred.Vr;
        Alpha_Amp=FP.Pred.Alpha_Amp*pi/180;
        Cycl=FP.Pred.Cycl;

        PPoint=sum(Cycl.*Vr./ds+1)*length(Alpha_Amp);
        
        alpha_h_pr=zeros(PPoint,1+Lag); 
        alpha_h_d_pr=zeros(PPoint,1);

        alpha_a_pr=zeros(PPoint,1+Lag); 
        alpha_a_d_pr=zeros(PPoint,1);
        
        CLh_pr=zeros(PPoint,1); 
        CMh_pr=zeros(PPoint,1); 
        CLa_pr=zeros(PPoint,1); 
        CMa_pr=zeros(PPoint,1); 

         Per=0;Num=0;Inc=sum((Cycl.*Vr./ds+1)); %Initiate
        for i=1:length(Alpha_Amp)
            for j=1:length(Vr)

            [alphaVr,alphaVr_d,CLhVr,CMhVr,CLaVr,CMaVr] = FP_Analy_SE (Vr(j),ds,Cycl(j),Alpha_Amp(i));
            NElemVR=length(alphaVr);
            
            if strcmp(FP.Pred.Excitation,'SinH')
            alpha_h_d_pr(1+Num+Per:Per+Num+NElemVR,1) =1./(1+alphaVr.^2).*alphaVr_d;
            CLh_pr(1+Num+Per:Per+Num+NElemVR)=CLhVr;%Lift coefficient due to vert disp
            CMh_pr(1+Num+Per:Per+Num+NElemVR)=CMhVr;%Moment coefficient due to rot 
            elseif strcmp(FP.Pred.Excitation,'SinA')     
            alpha_a_d_pr(1+Num+Per:Per+Num+NElemVR,1) =alphaVr_d;               
            CLa_pr(1+Num+Per:Per+Num+NElemVR)=CLaVr;%Lift coefficient due to vert disp
            CMa_pr(1+Num+Per:Per+Num+NElemVR)=CMaVr;%Moment coefficient due to rot
            end
              for ll=0:Lag %Regressors
                  if strcmp(FP.Pred.Excitation,'SinH')
                  alpha_h_pr(1+ll+Num+Per:Per+Num+NElemVR,1+ll)   =atan(alphaVr(1:NElemVR-ll,1));
                  alpha_h_pr(1+Num+Per:1+ll+Num+Per,1+ll)         =atan(alphaVr(NElemVR-ll:end,1));
                  elseif strcmp(FP.Pred.Excitation,'SinA')
                  alpha_a_pr(1+ll+Num+Per:Per+Num+NElemVR,1+ll)   =alphaVr(1:NElemVR-ll,1);
                  alpha_a_pr(1+Num+Per:1+ll+Num+Per,1+ll)         =alphaVr(NElemVR-ll:end,1);
                  end
              end
              Num=Num+NElemVR; %Reset
            end 
            Per=Per+Inc;Num=0; %Reset
        end
      
        CL_pr=CLh_pr+CLa_pr;
        CM_pr=CMh_pr+CMa_pr;
        
    elseif strcmp(FP.Pred.Excitation,'Stat') %Static excitation
        Alpha_Amp=FP.Pred.StatAlpha*pi/180;
        PPoint=(FP.Pred.StatTau./ds+1)*length(Alpha_Amp);
        
        alpha_h_pr=zeros(PPoint,1+Lag); 
        alpha_h_d_pr=zeros(PPoint,1);

        alpha_a_pr=zeros(PPoint,1+Lag); 
        alpha_a_d_pr=zeros(PPoint,1);
        
        CL_pr=zeros(PPoint,1);
        CM_pr=zeros(PPoint,1);    
        
        Inc=(FP.Pred.StatTau./ds+1); %Initiate
        for i=1:length(Alpha_Amp)
           alpha_a_pr(1+(i-1)*Inc:i*Inc,:)=Alpha_Amp(i);
           CL_pr(1+(i-1)*Inc:i*Inc)=Alpha_Amp(i)*2*pi; %Slope 2*pi
           CM_pr(1+(i-1)*Inc:i*Inc)=Alpha_Amp(i)*pi/2; %Slope pi/2            
        end
    
        
    elseif strcmp(FP.Pred.Excitation,'Aeroelastic')  %Flutter       
        M=[FP.Par.m FP.Par.I];
        C=4*pi*[FP.Par.fh FP.Par.fa].*FP.Par.psi.*M;
        K=4*pi^2*[FP.Par.fh FP.Par.fa].^2.*M;

        M=M.*FP.Par.U^2/FP.Par.B^2; % Analysis w.r.t. nondimensional time
        C=C.*FP.Par.U/FP.Par.B;     % Analysis w.r.t. nondlimensional time
        
        PPoint=floor(FP.Par.Redtime./ds);
        tau_predict=(0:ds:ds*(PPoint-1))'; %Time & Reduced time
        Phi_se=1-0.165.*exp(-0.089.*tau_predict)-0.335.*exp(-0.6.*tau_predict);% Wagner function
        
        u=zeros(PPoint,2); %Dispalcements
        u_d=zeros(PPoint,2); %Velocity
        u_2d=zeros(PPoint,2); %Acceleration
        
        u(1,1)=0.5; %Initial conditions
        
        CL_pr=zeros(PPoint,1);
        CM_pr=zeros(PPoint,1);
        Pres=1/2*FP.Par.rho*FP.Par.U^2*FP.Par.B.*[1 FP.Par.B];
        B=FP.Par.B;
        
       for i=1:PPoint %Time integration loop
        [~,h_conv]=ConvFFT(Phi_se(1:i),u_2d(1:i,1)./B,ds);
        [~,a_conv]=ConvFFT(Phi_se(1:i),u_d(1:i,2),ds);
        [~,a_d_conv]=ConvFFT(Phi_se(1:i),u_2d(1:i,2).*1/4,ds);
        Conv=h_conv+a_conv+a_d_conv;
                
        CL_pr(i)=2*pi*(-1/4*u_2d(i,1)/B-1/4*u_d(i,2)-Conv);
        CM_pr(i)=pi/2*(-1/4*u_d(i,2)-u_2d(i,2)/32   +Conv);

        p=Pres.*[CL_pr(i) CM_pr(i)];
        [u(i+1,:),u_d(i+1,:),u_2d(i+1,:)] = NewmarkSDOF(p,u(i,:),u_d(i,:),u_2d(i,:),1/4,1/2,ds,M,C,K);    
       end
       %Store extra
        FP.Pred.u=u; FP.Pred.u_d=u_d; FP.Pred.u_2d=u_2d;
        FP.Par.K=K; FP.Par.C=C;       FP.Par.M=M;
        FP.Par.Pres=Pres;
        
        alpha_h_pr=[];  
        alpha_h_d_pr=[];
        alpha_a_pr=[]; 
        alpha_a_d_pr=[]; 
    end

%% Store
%Prediction
FP.Pred.alpha_h=alpha_h_pr;
FP.Pred.alpha_h_d=alpha_h_d_pr;
FP.Pred.alpha_a=alpha_a_pr;
FP.Pred.alpha_a_d=alpha_a_d_pr;
FP.Pred.CL=CL_pr;
FP.Pred.CM=CM_pr;
FP.Pred.Samp=PPoint;
    
end

%% Self-excited forces due to sinusoidal effective angle
function [alpha,alpha_d,CLh,CMh,CLa,CMa] = FP_Analy_SE (Vr,ds,Cycl,Amp)
          K=2*pi/Vr; k=K/2; %Reduced frequency & Halfreduced frequency
          tau=(0:ds:(Vr*Cycl)).'; %Total reduced time for number of Cycl
          alpha=Amp*exp(1j.*K.*tau-1j*pi/2); %Sinusoidal alpha          
          alpha_d=(1j.*K.*Amp*exp(1j.*K.*tau-1j*pi/2));     %Angular velocity 
          alpha_2d=((1j.*K).^2.*Amp*exp(1j.*K.*tau-1j*pi/2));%Angular acceleration 
          
          Theo=conj(besselh(1,k))./(conj(besselh(1,k))+ 1j.*conj(besselh(0,k)));%Theodorsen function
          
          CLh=2*pi*real(-Theo*alpha-1/4*alpha_d);                     %Lift coefficient due to vert disp
          CMh=pi/2*real( Theo*alpha);                                      %Moment coefficient due to rot  
          CLa=2*pi*real(-Theo*(alpha+alpha_d*1/4)-1/4*alpha_d);            %Lift coefficient due to vert disp
          CMa=pi/2*real( Theo*(alpha+alpha_d*1/4)-1/4*alpha_d-alpha_2d/32);%Moment coefficient due to rot         
         
          alpha=real(alpha);
          alpha_d=real(alpha_d);
end

%% Fast convolution
function [Conv,Cent] = ConvFFT(M1,M2,dt)
M1 = [M1;zeros(length(M1),1)];
M2 = [M2;zeros(length(M2),1)];
  
Conv = ifft(fft(M1,[],1).*fft(M2,[],1),[],1);
Conv = Conv(1:length(M1)/2,:)*dt;
Cent=Conv(end,:);
end



