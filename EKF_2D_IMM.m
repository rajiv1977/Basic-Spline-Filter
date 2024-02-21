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
function [updatedX updatedP modeProb] = EKF_2D_IMM(X1_pri,P1_pri,X2_pri,P2_pri,NTsteps,m_noise_sd,Z,variance,F1,Tao1,F2,Tao2,G1,G2)
Markov = [0.9 0.1;0.2 0.8];
u_prior1 = 0.5;
u_prior2 = 0.5;
v_max = 10; 

for k = 1:NTsteps
    u11_predict = u_prior1 * Markov(1,1);
    u21_predict = u_prior2 * Markov(2,1);
    u12_predict = u_prior1 * Markov(1,2);
    u22_predict = u_prior2 * Markov(2,2);
    
    sum1 = u11_predict + u21_predict;
    u11_predict = u11_predict/sum1;
    u21_predict = u21_predict/sum1;
    
    sum2 = u12_predict + u22_predict;
    u12_predict = u12_predict/sum2;
    u22_predict = u22_predict/sum2;
    
    X1_0 = X1_pri * u11_predict + X2_pri * u21_predict; 
    X2_0 = X1_pri * u12_predict + X2_pri * u22_predict;
    
%     P2_pri_star = P2_pri;
%     P2_pri_star(2,2) = (v_max/2)^2;
    
    P1_0 = u11_predict * (P1_pri + (X1_pri - X1_0) * (X1_pri - X1_0)') + u21_predict * (P2_pri + (X2_pri - X1_0) * (X2_pri - X1_0)'); 
    P2_0 = u12_predict * (P1_pri + (X1_pri - X2_0) * (X1_pri - X2_0)') + u22_predict * (P2_pri + (X2_pri - X2_0) * (X2_pri - X2_0)');
    
    [X1_pri P1_pri likelihood1] = EKF(X1_0, P1_0, Z(k), F1, Tao1, variance, m_noise_sd^2, k,G1);
    [X2_pri P2_pri likelihood2] = EKF(X2_0, P2_0, Z(k), F2, Tao2, variance, m_noise_sd^2, k,G2);
    
    u_prior1 = (u11_predict + u21_predict) * likelihood1;
    u_prior2 = (u12_predict + u22_predict) * likelihood2;
    sum_u = u_prior1 + u_prior2;
    u_prior1 = u_prior1/sum_u;
    u_prior2 = u_prior2/sum_u;
    modeProb(1,k) = u_prior1;
    modeProb(2,k) = u_prior2;
    
    updatedX(:,k) = X1_pri * u_prior1 + X2_pri * u_prior2;
    updatedP(:,:,k) = u_prior1 * (P1_pri + (X1_pri - updatedX(:,k)) * (X1_pri - updatedX(:,k))')...
        + u_prior2 * (P2_pri + (X2_pri - updatedX(:,k)) * (X2_pri - updatedX(:,k))');
end
end

function [X_hat P_hat likelihood] = EKF(X, P, z, F, Tao, Q, R,k,G)
    X_pre = F * X + G * (k - 1);
    P_pre = F * P * F' + Tao * Q * Tao';
    
    H = find_H(X_pre(1), k);
    z_hat = atan(20/(X_pre(1) - 4*k));
    
    S = R + H * P_pre * H';
    v = z - z_hat;
    W = P_pre * H' * inv(S);
    P_hat = P_pre - W * S * W';
    likelihood = normpdf(v,0,sqrt(S));
    if likelihood < 1e-300
        likelihood = 1e-300;
    end
    X_hat = X_pre + W * v;
end

function H = find_H(x_pre,k)
    H = [-20/(400 + (x_pre - 4*k)^2), 0];
end