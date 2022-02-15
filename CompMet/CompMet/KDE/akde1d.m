function [pdf,grid]=akde1d(X,grid,gam)
%% fast adaptive kernel density estimation in one-dimension;
%  provides optimal accuracy/speed tradeoff, controlled with parameter "gam";
% INPUTS:   X  - data as a 'n' by '1' vector;
%
%         grid - (optional) mesh over which density is to be computed;
%                default mesh uses 2^12 points over range of data;
%
%          gam - (optional) cost/accuracy tradeoff parameter, where gam<n;
%                default value is gam=ceil(n^(1/3))+20; larger values
%                may result in better accuracy, but always reduce speed;
%                to speedup the code, reduce the value of "gam"; 
%
% OUTPUT: pdf - the value of the estimated density at 'grid'
%
%%  EXAMPLE:
%   data=[exp(randn(10^3,1))]; % log-normal sample
%   [pdf,grid]=akde1d(data); plot(grid,pdf)
%
% Note: If you need a very fast estimator use my "kde.m" function.
% This routine is more adaptive at the expense of speed. Use "gam"
% to control a speed/accuracy tradeoff. 
%
%%  Reference:
%  Kernel density estimation via diffusion
%  Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%  Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
[n,d]=size(X);
% begin scaling preprocessing
MAX=max(X,[],1);MIN=min(X,[],1);scaling=MAX-MIN;
MAX=MAX+scaling/10;MIN=MIN-scaling/10;scaling=MAX-MIN;
X=bsxfun(@minus,X,MIN);X=bsxfun(@rdivide,X,scaling);
if (nargin<2)|isempty(grid) % failing to provide grid
    grid=(MIN:scaling/(2^12-1):MAX)';
end
mesh=bsxfun(@minus,grid,MIN);mesh=bsxfun(@rdivide,mesh,scaling);
if nargin<3 % failing to provide speed/accuracy tradeoff
    gam=ceil(n^(1/3))+20;
end
% end preprocessing
% algorithm initialization
del=.2/n^(d/(d+4));perm=randperm(n);mu=X(perm(1:gam),:);
w=rand(1,gam);w=w/sum(w);Sig=del^2*rand(gam,d);
ent=-Inf;
for iter=1:1500 % begin algorithm
    Eold=ent;
    [w,mu,Sig,del,ent]=regEM(w,mu,Sig,del,X); % update parameters
    err=abs((ent-Eold)/ent); % stopping condition
%     fprintf('Iter.    Tol.      Bandwidth \n');
%     fprintf('%4i    %8.2e   %8.2e\n',iter,err,del);
%     fprintf('----------------------------\n');
    if (err<10^-5)|iter>200, break, end
end
% now output density values at grid
pdf = probfun(mesh,w,mu,Sig)/prod(scaling); % evaluate density
del=del*scaling; % adjust bandwidth for scaling
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,mu,Sig,del,ent]=regEM(w,mu,Sig,del,X)
[gam,d]=size(mu);[n,d]=size(X);
log_lh=zeros(n,gam); log_sig=log_lh;
for i=1:gam
    s=Sig(i,:);
    Xcentered = bsxfun(@minus, X, mu(i,:));
    xRinv = bsxfun(@rdivide, Xcentered.^2, s);
    xSig = sum(bsxfun(@rdivide, xRinv, s),2)+eps;
    log_lh(:,i)=-.5*sum(xRinv, 2)-.5*sum(log(s))...
        +log(w(i))-d*log(2*pi)/2-.5*del^2*sum(1./s);
    log_sig(:,i)=log_lh(:,i)+log(xSig);
end
maxll = max (log_lh,[],2); maxlsig = max (log_sig,[],2);
p= exp(bsxfun(@minus, log_lh, maxll));
psig=exp(bsxfun(@minus, log_sig, maxlsig));
density = sum(p,2); psigd=sum(psig,2);
logpdf=log(density)+maxll; logpsigd=log(psigd)+maxlsig;
p = bsxfun(@rdivide, p, density);% normalize classification prob.
ent=sum(logpdf); w=sum(p,1);
for i=find(w>0)
    mu(i,:)=p(:,i)'*X/w(i);  %compute mu's
    Xcentered = bsxfun(@minus, X,mu(i,:));
    Sig(i,:)=p(:,i)'*(Xcentered.^2)/w(i)+del^2; % compute sigmas
end
w=w/sum(w);
curv=mean(exp(logpsigd-logpdf));
del=1/(4*n*(4*pi)^(d/2)*curv)^(1/(d+2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=probfun(x,w,mu,Sig)
[gam,d]=size(mu);
out=0;
for k=1:gam
    S=Sig(k,:);
    xx=bsxfun(@minus, x,mu(k,:));
    xx=bsxfun(@rdivide,xx.^2,S);
    out=out+exp(-.5*sum(xx,2)+log(w(k))-.5*sum(log(S))-d*log(2*pi)/2);
end
end


