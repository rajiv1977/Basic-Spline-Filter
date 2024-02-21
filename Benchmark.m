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
function Benchmark()
clc;
clear all;
close all;
NTsteps = 20;
RMSE_s = zeros(1,NTsteps);
RMSE_p1 = zeros(1,NTsteps);
RMSE_p2 = zeros(1,NTsteps);
RMSE_p3 = zeros(1,NTsteps);
CRLB_M = zeros(1,NTsteps);
NEES_s = zeros(1,NTsteps);
NEES_p1 = zeros(1,NTsteps);
NEES_p2 = zeros(1,NTsteps);
NEES_p3 = zeros(1,NTsteps);


M = 1000;
DO_PLOT = 0;
SCALE = 1;

F = 1;
Tao = 1;
U = 0;
Yp = 50;

m_noise_sd = 0.5;%1*pi/180;
p_noise_sd = sqrt(10);%0.5;
variance = p_noise_sd^2;
MAX_ERROR = 1e20;

spOrder = 3;
x_min = -30*SCALE;%70;
x_max = 25*SCALE;%110;
NumberOfKnots = 10;
prior_knots = x_min:(x_max - x_min)/NumberOfKnots:x_max;
J_0 = 1/((x_max-x_min)^2/12);
dx= min(prior_knots):.1:max(prior_knots);
sp_Prior = spap2(augknt(prior_knots,spOrder), spOrder, dx, ones(1,length(dx))/(max(dx) - min(dx)));
NofSpline = sp_Prior.number;
%[X1lims V dx_sp IntOfBasicSP sp_Prior] = SP_preconstruction(NofSpline,variance,spOrder,NTsteps);

for m = 1:M
    %>>>>>>>>>>>>>>>> Generation of truths and measurements >>>>>>>>>>>>>>
    X0=0.1;%80;
    X1(1:NTsteps)=0;
    X1(1)=X0;
    for i = 2:NTsteps
        noise = randn * sqrt(variance);
        X1(i) = X1(i-1)/2 + SCALE*25*X1(i-1)/(1 + X1(i-1)^2) + 8*cos(1.2*i) + noise + U;
        %X1(i) = SCALE*X1(i-1) + U + noise;
    end
    
    for ii=1:NTsteps
        noise = m_noise_sd * randn;
        cdY1(ii)=X1(ii)^2/20 + noise;
        %cdY1(ii) = atan(Yp/(X1(ii) - 4 * ii)) + noise;
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



    %>>>>>>>>>>>>> Spline filter
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %updatedPDF = SF(sp_Prior,NTsteps,V,X1lims,spOrder,m_noise_sd,cdY1,IntOfBasicSP);
    updatedPDF = SF_movingKnots(sp_Prior,sp_Prior.knots,NTsteps,spOrder,m_noise_sd,cdY1,variance,NofSpline,SCALE,U,Yp);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    %>>>>>>>>>>>>> Particle filter >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    numOfParticle1 = 100;
    [x_PF1 w_PF1] = PF(cdY1,NTsteps,variance,m_noise_sd,numOfParticle1,SCALE,U,Yp,x_min,x_max);
    
    numOfParticle2 = 1000;
    [x_PF2 w_PF2] = PF(cdY1,NTsteps,variance,m_noise_sd,numOfParticle2,SCALE,U,Yp,x_min,x_max);
    
    numOfParticle3 = 10000;
    [x_PF3 w_PF3] = PF(cdY1,NTsteps,variance,m_noise_sd,numOfParticle3,SCALE,U,Yp,x_min,x_max);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    %>>>>>>>>>>>>> Numerical method >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%     dx1_numerical = 1;
%     [X1lims_numerical u] = numerical(cdY1,NTsteps,variance,m_noise_sd,dx1_numerical);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    %>>>>>>>>>>>>> Plot distributions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if DO_PLOT == 1
        for kk = 1:5
            figure
            fnplt(updatedPDF(kk));
            hold on
            plot(X1lims_numerical,u(kk,:),'.r');
            legend('Spline', 'Numerical');
        end

%         for kk = 1:1:NTsteps
%             figure
%             hold on;
%             hist(x_PF(kk,:),length(unique(x_PF(kk,:))));
%         end
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    %>>>>>>>>>>>>>>>>>>>>>>Find state estimates>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    for ll=1:NTsteps
%         X2mean(ll) = 0;
%         x = min(updatedPDF(ll).knots):.01:max(updatedPDF(ll).knots);
%         sum = 0;
%         for ii = 1:length(x)
%             X2mean(ll) = X2mean(ll) + max(0,fnval(updatedPDF(ll), x(ii))) * x(ii);
%             sum = sum + max(0,fnval(updatedPDF(ll), x(ii)));
%         end
%          X2mean(ll) = X2mean(ll)/sum;
         X2mean(ll) = findMean(updatedPDF(ll).coefs, updatedPDF(ll).knots);
        xPFmean1(ll) = (x_PF1(ll,:)*w_PF1(ll,:).');
        xPFmean2(ll) = (x_PF2(ll,:)*w_PF2(ll,:).');
        xPFmean3(ll) = (x_PF3(ll,:)*w_PF3(ll,:).');
    end
    xSP = X2mean;
    if DO_PLOT == 1
        figure;
        plot(1:NTsteps+1,X1,'-.',1:NTsteps,X1mean,'-o',1:NTsteps,X2mean,'-+',1:NTsteps,xPFmean,'-v');
        legend('Truth','Numerical','Spline','PF');
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    
    CRLB = crlb(X1,m_noise_sd,SCALE*F,Tao,variance,J_0,Yp);
    CRLB_M = CRLB_M + CRLB;
    for k = 1:NTsteps
        RMSE_s(k) = RMSE_s(k) + (min((xSP(1,k) - X1(1,k)), MAX_ERROR))^2;
        RMSE_p1(k) = RMSE_p1(k) + (min((xPFmean1(1,k) - X1(1,k)), MAX_ERROR))^2;
        RMSE_p2(k) = RMSE_p2(k) + (min((xPFmean2(1,k) - X1(1,k)), MAX_ERROR))^2;
        RMSE_p3(k) = RMSE_p3(k) + (min((xPFmean3(1,k) - X1(1,k)), MAX_ERROR))^2;
        NEES_s(k) = NEES_s(k) + ((xSP(1,k) - X1(1,k)))^2/CRLB(k);
        NEES_p1(k) = NEES_p1(k) + ((xPFmean1(1,k) - X1(1,k)))^2/CRLB(k);
        NEES_p2(k) = NEES_p2(k) + ((xPFmean2(1,k) - X1(1,k)))^2/CRLB(k);
        NEES_p3(k) = NEES_p3(k) + ((xPFmean3(1,k) - X1(1,k)))^2/CRLB(k);
    end
    m
end
RMSE_s = (RMSE_s/M).^(.5);
RMSE_p1 = (RMSE_p1/M).^(.5);
RMSE_p2 = (RMSE_p2/M).^(.5);
RMSE_p3 = (RMSE_p3/M).^(.5);
CRLB_M = (CRLB_M/M).^(.5);
figure
plot(RMSE_s,'-ob');
hold on;
plot(RMSE_p1,'-vr');
hold on;
plot(RMSE_p2,'-+g');
hold on;
plot(RMSE_p3,'-sm');
hold on;
plot(CRLB_M,'-dk');
hold on;
legend('SP','PF1','PF2','PF3','PCRLB');
title('RMSE');

NEES_s = NEES_s/M;
NEES_p1 = NEES_p1/M;
NEES_p2 = NEES_p2/M;
NEES_p3 = NEES_p3/M;
figure
plot(NEES_s,'-ob');
hold on;
plot(NEES_p1,'-vr');
hold on;
plot(NEES_p2,'-+g')
hold on;
plot(NEES_p3,'-sm');
hold on;
legend('SP','PF1','PF2','PF3');
title('NEES');
disp(num2str(sum(RMSE_s)/NTsteps));
disp(num2str(sum(RMSE_p1)/NTsteps));
disp(num2str(sum(RMSE_p2)/NTsteps));
disp(num2str(sum(RMSE_p3)/NTsteps));
save SplineResult1D.mat RMSE_s RMSE_p1 RMSE_p2 RMSE_p3 CRLB_M NEES_s NEES_p1 NEES_p2 NEES_p3
end