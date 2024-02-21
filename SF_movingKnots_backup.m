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
function [updatedPDF] = SF_movingKnots(prior_sp,prior_knots,K,spOrder,m_noise_sd,Z,variance,M)
N = spOrder;
knots = prior_knots;
prior = prior_sp;
dx=.1;
for k = 1:K
    %>>> predict >>>
     cur_knots = knots/2 + 25 * knots./(1 + knots.^2) + 8*cos(1.2*k);
     %cur_knots = knots;
     cur_knots = sort(cur_knots);
    %cur_knots = knots;
%     [x, y] = ndgrid(knots,cur_knots);
%     delta = y - x;
%     P = (1/(sqrt(2*pi*variance))) * exp(-(delta.^2/(2*variance)));
    x_range = min(knots):dx:max(knots);    
    cur_x_range = x_range/2 + 25 * x_range./(1 + x_range.^2) + 8*cos(1.2*k);
    %cur_x_range = x_range;
    cur_x_range = sort(cur_x_range);
    P=zeros(length(cur_x_range),length(cur_x_range));
    for i = 1:length(cur_x_range)
        for j = 1:length(cur_x_range)
            delta =  x_range(j) - cur_x_range(i);
            P(j,i) = (1/(sqrt(2*pi*variance))) * exp(-(delta)^2/(2*variance));
        end
    end
    %spTran = spapi({N,N},{knots,cur_knots},P);
    %spTran = spap2({M - N + 1,M - N + 1},[N, N],{knots,cur_knots},P);
    spTran = spap2({augknt(unique(knots),N), augknt(unique(cur_knots),N)},{N,N},{x_range,cur_x_range},P);
    tran_coefs = spTran.coefs;
    
    
    %>>>>>>>>>>>>>>>>> From donna >>>>>>>>>>>>>>>>>>>>>>>
%     for i = 1:length(cur_x_range)
%         Pnorm(i,:) = P(i,:)/sum(P(i,:));
%     end 
%     for t = 1:length(x_range)
%         Probability(t) = spaps(x_range, Pnorm(:,t), 1e-10);    
%     end
%     
%     
%     for uu = 1:length(x_range)
%         tpred(uu)= fncmb(prior ,'*', Probability(uu));
%         predsum(uu) = fnval(fnint(tpred(uu)), tpred(uu).breaks(end));
%     end
% 
%     predsumforplot = predsum.';
% 
%     predictedval = predsumforplot/sum(predsumforplot*dx);
%     predicted_donna = spaps(x_range,predictedval, 1e-20);
    
    
%     Corr=(1/(sqrt(2*pi)*cdsigma_y1))*exp(-(1/(2*cdsigma_y1^2))*(cdY1(kk+1)-atan(X1lims)).^2);
%     if USE_ONLY_SIGNIFICANYT_POINTS == 1
%         sign_points = find_significant_points(Corr, lx1);
%     else
%         sign_points = 1:lx1;
%     end
%     likelihoodDistribution = spaps(X1lims(sign_points), Corr(sign_points), 1e-10);
%     updatedPDF(kk) = fncmb(predicted, '*', likelihoodDistribution);
%     cum_pro = fnval(fnint(updatedPDF(kk)), updatedPDF(kk).breaks(end));    
%     updatedPDF(kk).coefs = updatedPDF(kk).coefs/cum_pro;
%     predicted = updatedPDF(kk);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    %sp_Basic = spapi(N,knots,ones(1,length(knots)));
    %sp_Basic = spap2(M-N+1,N,knots,ones(1,length(knots)));
    sp_Basic = spap2(augknt(unique(knots),N),N,x_range,ones(1,length(x_range)));
    
    r=M;
    numOfBasicSpline = r;
    IntOfBasicSP = zeros(r,r);
    coefs = sp_Basic.coefs;
    L = length(coefs);
    for i = 1:r
        coef_temp_i = zeros(1,L);
        coef_temp_i(i) = coefs(i); 
        sp_temp_i = sp_Basic;
        sp_temp_i.coefs = coef_temp_i;
        for j = i:min(i+N,r)
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
    
    predict = prior;
    pri_coefs = prior.coefs;
    pre_coefs = zeros(1,numOfBasicSpline);
    for i = 1:numOfBasicSpline
        for j = 1:numOfBasicSpline
            temp_t = pri_coefs * IntOfBasicSP(j,:)';
            pre_coefs(i) = pre_coefs(i) + temp_t * tran_coefs(1,i,j);
        end
        pre_coefs(i) = pri_coefs * V(i,:)';
    end
    
    %temp = spapi(N, cur_knots, ones(1,length(knots)));
    %temp = spap2(M-N+1,N, cur_knots, ones(1,length(knots)));
    %temp = spap2(cur_knots,N, cur_knots, ones(1,length(cur_knots)));
    x_range = min(cur_knots):dx:max(cur_knots); 
    temp = spap2(augknt(unique(cur_knots),N),N,x_range,ones(1,length(x_range)));
    
    predict = temp;
    predict.coefs = pre_coefs;
    %>>>>> update >>>>
    x = min(cur_knots):dx:max(cur_knots);
    delta = Z(k) - x.^2/20;
    %delta = Z(k) - x;
    Corr = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
    %likelihoodDistribution = spapi({N},{cur_knots},Corr);
    %likelihoodDistribution = spap2(M-N+1,N,x,Corr);
    likelihoodDistribution = spap2(x,N,x,Corr);
    
    updatedPDF(k) = fncmb(predict, '*', likelihoodDistribution);
    cum_pro = fnval(fnint(updatedPDF(k)), updatedPDF(k).breaks(end));
    updatedPDF(k).coefs = updatedPDF(k).coefs/cum_pro;
    updatedDistribution = updatedPDF(k);
    vals = fnval(updatedDistribution, x);
    [dump ind] = sort(cur_knots);
    %prior = spap2(augknt(unique(cur_knots),N),N,x,vals);
    prior = spap2(M-N+1,N,x,vals);
    %prior = spap2(M-N+1,N,cur_knots,vals);
    knots = unique(prior.knots);
    fnplt(prior);
end