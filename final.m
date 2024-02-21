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

spOrder = 5;
multi = 1; 

USE_ONLY_SIGNIFICANYT_POINTS = 0;
dt=0.01;
X0=-0.50;
NTsteps=200;
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
x = X1lims; 
y = x; 
[xx,yy] = ndgrid(x,y);
z = P;
sp_Transition = spap2({length(x) - spOrder + 1, length(x) - spOrder + 1},[spOrder spOrder],{x,y},z);
tran_coefs = sp_Transition.coefs;

%Here we are essentially trying to multiply the Transition probabilty
%spline with the prior.

%Initialize pdf
%This part is spline version.
t(1:lx1) = 1;
% the prior is the first rox of t matrix as it represents
sp_Prior = spap2(length(x) - spOrder + 1, spOrder, X1lims, t(1,:));

t(1:lx1) = 1;
sp_Basic = spap2(length(x) - spOrder + 1, spOrder, X1lims, t(1,:));

r = length(sp_Basic.coefs);
IntOfBasicSP = zeros(r,r);
coefs = sp_Basic.coefs;
L = length(coefs);

for i = 1:r
    coef_temp_i = zeros(1,L);
    coef_temp_i(i) = coefs(i); 
    sp_temp_i = sp_Basic;
    sp_temp_i.coefs = coef_temp_i;
    for j = i:min(i+spOrder,r)
        coef_temp_j = zeros(1,L);
        coef_temp_j(j) = coefs(j);
        sp_temp_j = sp_Basic;
        sp_temp_j.coefs = coef_temp_j;
        product = fncmb(sp_temp_i ,'*', sp_temp_j);
        IntOfBasicSP(i,j) = fnval(fnint(product), product.breaks(end));
        IntOfBasicSP(j,i) = IntOfBasicSP(i,j);
    end
end

V = zeros(r,r);
for i = 1:r
    for j = 1:r
        temp_i(j) = tran_coefs(1,i,j);
    end
    for t = 1:r
        V(i,t) = temp_i * IntOfBasicSP(:,t);
    end
end
%for ss = 1:NTsteps
%  multiplyingpriorandtransition;
prior = sp_Prior;
pri_coefs = prior.coefs;
numOfBasicSpline = length(pri_coefs);
for kk = 1:NTsteps
    predict = prior;
    pri_coefs = prior.coefs;
    pre_coefs = zeros(1,numOfBasicSpline);
    for i = 1:numOfBasicSpline
%         for j = 1:numOfBasicSpline
%             temp_t = pri_coefs * IntOfBasicSP(j,:)';
%             pre_coefs(i) = pre_coefs(i) + temp_t * tran_coefs(1,i,j);
%         end
        pre_coefs(i) = pri_coefs * V(i,:)';
    end
    predict.coefs = pre_coefs;

    Corr=(1/(sqrt(2*pi)*cdsigma_y1))*exp(-(1/(2*cdsigma_y1^2))*(cdY1(kk+1)-atan(X1lims)).^2);
    likelihoodDistribution = spap2(length(X1lims) - spOrder + 1, spOrder, X1lims, Corr);
    updatedPDF(kk) = fncmb(predict, '*', likelihoodDistribution);
    cum_pro = fnval(fnint(updatedPDF(kk)), updatedPDF(kk).breaks(end));    
    updatedPDF(kk).coefs = updatedPDF(kk).coefs/cum_pro;
    updatedDistribution = updatedPDF(kk);
    vals = fnval(updatedDistribution, X1lims);
    %prior = spap2(length(X1lims) - spOrder + 1, spOrder, X1lims, vals);
    prior = predict;
    prior.coefs = slvblk(spcol(predict.knots, spOrder, X1lims','slvblk','noderiv'),vals')';
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

for kk = 1:50:NTsteps
figure
fnplt(updatedPDF(kk));
hold on
plot(X1lims,u(kk+1,:),'r');
legend('Spline', 'Numerical')
end
%Conditional Mean

for ll=1:NTsteps
X1mean(ll)=dot(squeeze(u(ll,:)),X1lims)*dx1;
X2mean(ll) = 0;
for ii = 1:lx1
    X2mean(ll) = X2mean(ll) + dot(fnval(updatedPDF(ll), X1lims(ii)), X1lims(ii));
end
X2mean(ll) = X2mean(ll)*dx1;
end

figure(4);
plot(1:NTsteps,X1(2:end),'-',1:NTsteps,X1mean,'-o',1:NTsteps,X2mean,'r');
legend('Truth','Numerical','Spline');
axis([0 200 min(X1) max(X1)]);
end