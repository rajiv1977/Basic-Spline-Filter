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
function Benchmark_2D_IMM()
clc;
clear all;
close all;
NTsteps1 = 15;
NTsteps2 = 10;

m_t = [ones(1,NTsteps1) zeros(1,NTsteps2)];%true mode prob

NTsteps = NTsteps1 + NTsteps2;
RMSE_s = zeros(1,NTsteps);
RMSE_p1 = zeros(1,NTsteps);
RMSE_p2 = zeros(1,NTsteps);
CRLB_M = zeros(1,NTsteps);
NEES_s = zeros(1,NTsteps);
NEES_p1 = zeros(1,NTsteps);
NEES_p2 = zeros(1,NTsteps);
M = 1;
DO_PLOT = 0;
SCALE = 1;
Yp = 50;
J_0 = 1/100^2;

m_noise_sd = 1*pi/180;
p_noise_sd = .8;
variance = p_noise_sd^2;

spOrder = 3;
prior_knots1{1} = 0:4:40;
prior_knots1{2} = 0:2.5:10;
dx= min(prior_knots1{1}):.1:max(prior_knots1{1});
dv= min(prior_knots1{2}):.01:max(prior_knots1{2});
sp_Prior1 = spap2({augknt(prior_knots1{1},spOrder),augknt(prior_knots1{2},spOrder)}, [spOrder spOrder], {dx dv}, ...
    ones(length(dx),length(dv))/(max(dx) - min(dx))/(max(dv) - min(dv)));
NofSpline = sp_Prior1.number;

prior_knots2{1} = 0:4:40;
prior_knots2{2} = 0:.5:2;
dx= min(prior_knots2{1}):.1:max(prior_knots2{1});
dv= min(prior_knots2{2}):.01:max(prior_knots2{2});
sp_Prior2 = spap2({augknt(prior_knots2{1},spOrder),augknt(prior_knots2{2},spOrder)}, [spOrder spOrder], {dx dv}, ...
    ones(length(dx),length(dv))/(max(dx) - min(dx))/(max(dv) - min(dv)));

%[X1lims V dx_sp IntOfBasicSP sp_Prior] = SP_preconstruction(NofSpline,variance,spOrder,NTsteps);
F1 = [1 1;0 1];
Tao1 = [.5;.01];
G1 = [0;0];

F2 = [1 1;0 1];
Tao2 = [.5;.01];
G2 = [-0.4;0];

for m = 1:M
    %>>>>>>>>>>>>>>>> Generation of truths and measurements >>>>>>>>>>>>>>
    X0=[30;5];
    X1(:,1)=X0;
    for i = 2:NTsteps1
        noise = randn * sqrt(variance)/100;
        X1(:,i) = F1*X1(:,i-1) + G1 * (i - 1) + Tao1 *noise;
    end
    for i = NTsteps1+1:NTsteps1+NTsteps2
        noise = randn * sqrt(variance)/100;
        X1(:,i) = F2*X1(:,i-1) + G2 * (i - 1) + Tao2 *noise;
    end
    
    for ii=1:NTsteps
        noise = (m_noise_sd) * randn/10;
        %cdY1(ii) = atan2(20,(X1(1,ii) - 4 * ii)) + noise;
        cdY1(ii) = atan(20/(X1(1,ii) - 4 * ii)) + noise;
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %>>>>>>>>>>>>>>> EKF >>>>>>>>>>>>>>>>>>>>>>>>>>>
    X1_pri = [30;5];
    P1_pri = [400 0;0 25];
    X2_pri = [30;5];
    P2_pri = [400 0;0 25];
    X1
    [updatedX updatedP modeProb] = EKF_2D_IMM(X1_pri,P1_pri,X2_pri,P2_pri,NTsteps,m_noise_sd,cdY1,variance,F1,Tao1,F2,Tao2,G1,G2);
    modeProb
    t=1:NTsteps;
    plot(t,modeProb(1,:),'-r+')
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  


    %>>>>>>>>>>>>> Spline filter
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %updatedPDF = SF(sp_Prior,NTsteps,V,X1lims,spOrder,m_noise_sd,cdY1,IntOfBasicSP);
    %updatedPDF = SF_movingKnots_2D(sp_Prior,sp_Prior.knots,NTsteps,spOrder,m_noise_sd,cdY1,variance,NofSpline,SCALE,DO_PLOT);
    [updatedPDF modeProb] = SF_movingKnots_2D_IMM(sp_Prior1,sp_Prior1.knots,sp_Prior1,sp_Prior1.knots,NTsteps,spOrder,m_noise_sd,cdY1,variance,...
        NofSpline,SCALE,DO_PLOT,F1,Tao1,G1,F2,Tao2,G2);
    if M == 1
        figure
        plot(1:NTsteps,m_t,'-bo')
        hold on;
        plot(1:NTsteps,modeProb(1,:),'-r+')
        hold on;
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    %>>>>>>>>>>>>> Particle filter >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    numOfParticle1 = 100;
    [x_PF1 w_PF1] = PF_2D(cdY1,NTsteps,variance,m_noise_sd,numOfParticle1,SCALE);
    
    numOfParticle2 = 1000;
    [x_PF2 w_PF2] = PF_2D(cdY1,NTsteps,variance,m_noise_sd,numOfParticle2,SCALE);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    %>>>>>>>>>>>>>>>>>>>>>>Find state estimates>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    xSP = zeros(2,NTsteps);
    for ll=1:NTsteps
        x{1} = min(updatedPDF(ll).knots{1}):1:max(updatedPDF(ll).knots{1});
        x{2} = min(updatedPDF(ll).knots{2}):.1:max(updatedPDF(ll).knots{2});
        sum = 0;
        for ii = 1:length(x{1})
            for jj = 1:length(x{2})
                xSP(:,ll) = xSP(:,ll) + max(0,fnval(updatedPDF(ll), {x{1}(ii) x{2}(jj)})) * [x{1}(ii);x{2}(jj)];
                sum = sum + max(0,fnval(updatedPDF(ll), {x{1}(ii) x{2}(jj)}));
            end
        end
        xSP(:,ll) = xSP(:,ll)/sum;

        PF1 = zeros(length(x_PF1(ll,:,1)),length(x_PF1(ll,1,:)));
        for i = 1:length(x_PF1(ll,:,1))
            for j = 1:length(x_PF1(ll,1,:))
                PF1(i,j) = x_PF1(ll,i,j);
            end
        end
        xPFmean1(:,ll) = zeros(2,1);
        for i = 1:length(w_PF1(ll,:))
            xPFmean1(:,ll) = xPFmean1(:,ll) + PF1(:,i)*w_PF1(ll,i);
        end
        PF2 = zeros(length(x_PF2(ll,:,1)),length(x_PF2(ll,1,:)));
        for i = 1:length(x_PF2(ll,:,1))
            for j = 1:length(x_PF2(ll,1,:))
                PF2(i,j) = x_PF2(ll,i,j);
            end
        end
        xPFmean2(:,ll) = zeros(2,1);
        for i = 1:length(w_PF2(ll,:))
            xPFmean2(:,ll) = xPFmean2(:,ll) + PF2(:,i)*w_PF2(ll,i);
        end
        ll
    end
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if M == 1
        figure
        plot(X1(1,:),'-ob');
        hold on;
        plot(xSP(1,:),'-^m');
        hold on;
        plot(xPFmean1(1,:),'-vr');
        hold on;
        plot(xPFmean1(1,:),'-+g');
        hold on;
        legend('Truth','SP','PF1','PF2');
        title('Position');
        figure
        plot(X1(2,:),'-ob');
        hold on;
        plot(xSP(2,:),'-^m');
        hold on;
        plot(xPFmean1(2,:),'-vr');
        hold on;
        plot(xPFmean1(2,:),'-+g');
        hold on;
        legend('Truth','SP','PF1','PF2');
        title('Velocity');
    end
    CRLB = crlb(X1,m_noise_sd,F1,Tao1,variance,J_0,Yp);
    CRLB_M = CRLB_M + CRLB;
    for k = 1:NTsteps
        RMSE_s(k) = RMSE_s(k) + ((xSP(1,k) - X1(1,k)))^2;
        RMSE_p1(k) = RMSE_p1(k) + ((xPFmean1(1,k) - X1(1,k)))^2;
        RMSE_p2(k) = RMSE_p2(k) + ((xPFmean2(1,k) - X1(1,k)))^2;
        NEES_s(k) = NEES_s(k) + ((xSP(1,k) - X1(1,k)))^2/CRLB(k);
        NEES_p1(k) = NEES_p1(k) + ((xPFmean1(1,k) - X1(1,k)))^2/CRLB(k);
        NEES_p2(k) = NEES_p2(k) + ((xPFmean2(1,k) - X1(1,k)))^2/CRLB(k);
    end
    m
end
RMSE_s = (RMSE_s/M).^(.5);
RMSE_p1 = (RMSE_p1/M).^(.5);
RMSE_p2 = (RMSE_p2/M).^(.5);
CRLB_M = (CRLB_M/M).^(.5);
figure
plot(RMSE_s,'-ob');
hold on;
plot(RMSE_p1,'-vr');
hold on;
plot(RMSE_p2,'-+g');
hold on;
plot(CRLB_M,'-dk');
hold on;
legend('SP','PF1','PF2','CRLB');
title('RMSE');

NEES_s = NEES_s/M;
NEES_p1 = NEES_p1/M;
NEES_p2 = NEES_p2/M;
figure
plot(NEES_s,'-ob');
hold on;
plot(NEES_p1,'-vr');
hold on;
plot(NEES_p2,'-+g');
hold on;
legend('SP','PF1','PF2');
title('NEES');

save SplineResult.mat RMSE_s RMSE_p1 RMSE_p2 CRLB_M NEES_s NEES_p1 NEES_p2
end