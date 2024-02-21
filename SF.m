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
function [updatedPDF] = SF(prior,NTsteps,V,X1lims,spOrder,m_noise_sd,cdY1,IntOfBasicSP)
pri_coefs = prior.coefs;
numOfBasicSpline = length(pri_coefs);
for kk = 1:NTsteps
    predict = prior;
    pri_coefs = prior.coefs;
    pre_coefs = zeros(1,numOfBasicSpline);
    for i = 1:numOfBasicSpline
        for j = 1:numOfBasicSpline
            V_temp(j) = V(kk,i,j);
        end
        pre_coefs(i) = pri_coefs * V_temp';
    end
    predict.coefs = pre_coefs;

    x=X1lims;
    %delta = cdY1(kk) - x.^2/20;
    delta = cdY1(kk) - x;
    Corr = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
    likelihoodDistribution = spap2(length(X1lims) - spOrder + 1, spOrder, X1lims, Corr);
    
    updatedPDF(kk) = fncmb(predict, '*', likelihoodDistribution);
    cum_pro = fnval(fnint(updatedPDF(kk)), updatedPDF(kk).breaks(end));
    
    
    for j = 1:numOfBasicSpline
        temp1(j)=pre_coefs*IntOfBasicSP(:,j);
    end
    cum_pro1 = likelihoodDistribution.coefs * temp1';
    
    
    updatedPDF(kk).coefs = updatedPDF(kk).coefs/cum_pro;
    updatedDistribution = updatedPDF(kk);
    vals = fnval(updatedDistribution, X1lims);
    %prior = spap2(length(X1lims) - spOrder + 1, spOrder, X1lims, vals);
    prior = predict;
    prior.coefs = slvblk(spcol(predict.knots, spOrder, X1lims','slvblk','noderiv'),vals')';
end