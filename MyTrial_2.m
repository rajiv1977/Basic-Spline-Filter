%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      Spline Filter                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   S.Rajiv, B. Balaji, R.Tharmarasa,  and T.Kirubarajan                    %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%           sithirr@mcmaster.ca, tharman@grads.ece.mcmaster.ca, kiruba@mcmaster.ca          %
%                                                                                           %
%                                 B. Balaji and M.McDonald                                  %
%              Defence R&D Canada, 3701 Carling Avenue, Ottawa, ON K1A 0Z4, Canada.         %
%						   bhashyam.balaji@drdc-rddc.gc.ca                                  %
%                           mike.mcdonald@drdc-rddc.gc.ca                                   %
%                                                                                           %
%                                       M.Pelletier                                         %
%                           FLIR - Radars, Laval, QC, Canada.                               %
%                               Michel.Pelletier@flir.com                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MyTrial()

clear all;
close all;

%Signal Part

USE_ONLY_SIGNIFICANYT_POINTS = 0;
dt=0.01;
X0=-0.50;
NTsteps=10;
X1(1:NTsteps+1)=0;
X1(1)=X0;
sigma_x=0.3;
for ii=1:NTsteps
    Winc=sqrt(dt)*randn;
    fx1=1.2*cos(3*X1(ii));
    gx1=sigma_x;
    X1(ii+1)=X1(ii)+dt*fx1+Winc*gx1;
end
figure(1);
plot(X1)
%Measurement Part: Continuous-Discrete Model

cdY0=0;
cdY1(1:NTsteps+1)=0;
cdY1(1)=cdY0;
cdsigma_y1=0.05; %0.05, 0.0125 in C-C
for ii=1:NTsteps
    cdWincy=sqrt(dt)*randn;
    %cdf1=X1(ii+1);
    cdgy1=cdsigma_y1;
    cdY1(ii+1)=atan(X1(ii+1))+cdgy1*cdWincy;
end
figure(3);
plot(cdY1)

% Compute the transition probability density

dx1=0.01;
X1lims=-1.5:dx1:1.5;

lx1=length(X1lims);

r=0.5; %Symmetric discretization

for ii=1:lx1
    for jj=1:lx1
        x=X1lims(ii);
        y=X1lims(jj);
        X1r=X1lims(jj)+r*(X1lims(ii)-X1lims(jj));
        P(jj,ii)=(1./(sqrt(2*pi*dt)*sigma_x)).*...
            exp(-0.5*dt*(1/sigma_x^2)*( (x-y)./dt-1.2*cos(3*X1r)).^2+...
            r*dt*3.6*sin(3*X1r));
    end
end
% P(jj,ii) is the probability of going from jj to ii.
% Coverting Transition Probability to a 2D Spline
% The values of x are stored in X1lims vector and the transiiton
% probabilites are stored in P matrix.
% Points ={X1lims,X1lims};

%definingtransitionsplines;
for ii = 1: lx1
    Pnorm(ii,:) = P(ii,:)/sum(P(ii,:));
end
for tt = 1:lx1
    if USE_ONLY_SIGNIFICANYT_POINTS == 1
        [sign_points] = find_significant_points(Pnorm(:,tt), lx1);
    else
        sign_points = 1:lx1;
    end
    Probability(tt) = spaps(X1lims(sign_points), Pnorm(sign_points,tt), 1e-10);    
end
%Here we are essentially trying to multiply the Transition probabilty
%spline with the prior.

%Initialize pdf
%This part is spline version.
t(1:lx1)=1;
% the prior is the first rox of t matrix as it represents
prior = spaps(X1lims, t(1,:), 1e-10);
%for ss = 1:NTsteps
%  multiplyingpriorandtransition;
predicted = prior;
for kk = 1:NTsteps
    for uu = 1: lx1
        tpred(uu)= fncmb( predicted ,'*', Probability(uu));
        predsum(uu) = fnval(fnint(tpred(uu)), tpred(uu).breaks(end));
    end

    predsumforplot = predsum.';

    predictedval = predsumforplot/sum(predsumforplot*dx1);
    if USE_ONLY_SIGNIFICANYT_POINTS == 1
        sign_points = find_significant_points(predictedval', lx1);
    else
        sign_points = 1:lx1;
    end
    predicted = spaps(X1lims(sign_points),predictedval(sign_points), 1e-20);

    Corr=(1/(sqrt(2*pi)*cdsigma_y1))*exp(-(1/(2*cdsigma_y1^2))*(cdY1(kk+1)-atan(X1lims)).^2);
    if USE_ONLY_SIGNIFICANYT_POINTS == 1
        sign_points = find_significant_points(Corr, lx1);
    else
        sign_points = 1:lx1;
    end
    likelihoodDistribution = spaps(X1lims(sign_points), Corr(sign_points), 1e-10);
    updatedPDF(kk) = fncmb(predicted, '*', likelihoodDistribution);
    cum_pro = fnval(fnint(updatedPDF(kk)), updatedPDF(kk).breaks(end));    
    updatedPDF(kk).coefs = updatedPDF(kk).coefs/cum_pro;
    predicted = updatedPDF(kk);
end

%The correction part has to come in here.



% To normalize the predicted values to accomodate for the transition
% probability matrix having values greater than one.


%Original from Balaji paper
u(1:NTsteps,1:lx1)=1;

for kk=1:NTsteps
    upred=P(:,:)*u(kk,:).'*dx1; %Prediction
    Corr=(1/(sqrt(2*pi)*cdsigma_y1))*exp(-(1/(2*cdsigma_y1^2))*(cdY1(kk+1)-atan(X1lims)).^2);
    ucorr=Corr.*upred.'; %Correction
    u(kk+1,:)=ucorr/(sum(ucorr*dx1)); %Normalized
end

figure
fnplt(updatedPDF(NTsteps));
hold on
plot(X1lims,u(kk+1,:),'r');
legend('Spline', 'Numerical')

%Conditional Mean

for ll=1:NTsteps
X1mean(ll)=dot(squeeze(u(ll,:)),X1lims)*dx1;
X2mean(ll) = 0;
for ii = 1:lx1
    X2mean(ll) = X2mean(ll) + dot(fnval(updatedPDF(ll), updatedPDF(ll).breaks(ii)),X1lims(ii));
end
X2mean(ll) = X2mean(ll)*dx1;
end

figure(4);
plot(1:NTsteps,X1(2:end),'-',1:NTsteps,X1mean,'o',1:NTsteps,X2mean,'r');
legend('Truth','Numerical','Spline')

end 

function [sign_points] = find_significant_points(P, lx1)

sign_points = find(P > 0.01);
if sign_points(1) > 2
    sign_points = [floor(sign_points(1)/2) sign_points];
end
if isempty(sign_points == 1)
    sign_points = [1 sign_points];
end
if sign_points(end) < lx1-2
    sign_points = [sign_points ceil((lx1 + sign_points(end))/2)];
end
if isempty(sign_points == lx1)
    sign_points = [sign_points lx1];
end

end
