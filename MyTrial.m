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
NTsteps=10%200;
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
        X1r=X1lims(ii)+r*(X1lims(jj)-X1lims(ii));
        P(jj,ii)=(1./(sqrt(2*pi*dt)*sigma_x)).*...
            exp(-0.5*dt*(1/sigma_x^2)*( (x-y)./dt-1.2*cos(3*X1r)).^2+...
            r*dt*3.6*sin(3*X1r));
    end
end

%This part is spline version.
% Coverting Transition Probability to a 2D Spline
% The values of x are stored in X1lims vector and the transiiton
% probabilites are stored in P matrix.
Points ={X1lims,X1lims};

%definingtransitionsplines;
for ii = 1: lx1
    Pnorm(ii,:) = P(ii,:)/sum(P(ii,:));
end
for tt = 1:lx1
    if USE_ONLY_SIGNIFICANYT_POINTS == 1
        [sign_points] = find_significant_points(P(tt,:), lx1);
    else
        sign_points = 1:lx1;
    end
    Probability(tt) = spaps(X1lims(sign_points), P(tt,sign_points), 1e-10);    
end
%Here we are essentially trying to multiply the Transition probabilty
%spline with the prior.

%Initialize pdf
%This part is spline version.
t(1:lx1)=1;
% the prior is the first rox of t matrix as it represents
%  P(i,j) means the probability of going from j to i in this case.
prior = spaps(X1lims, t(1,:), 1e-10);
%for ss = 1:NTsteps
%  multiplyingpriorandtransition;
predicted = prior;
for kk = 1:NTsteps
    for uu = 1: lx1
        tpred(uu)= fncmb( predicted ,'*', Probability(uu));
        predsum(uu) = fnval(fnint(tpred(uu)), tpred(uu).breaks(end));
    end
end
if sign_points(end) < lx1-2
    sign_points = [sign_points ceil((lx1 + sign_points(end))/2)];
end
if isempty(sign_points == lx1)
    sign_points = [sign_points lx1];
end

end
