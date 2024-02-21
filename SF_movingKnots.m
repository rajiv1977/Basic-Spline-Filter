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
function [updatedPDF] = SF_movingKnots(prior_sp,prior_knots,K,spOrder,m_noise_sd,Z,variance,M,SCALE,U,Yp)
SHOW_FIG = 0;


N = spOrder;
knots = prior_knots;
prior = prior_sp;
NN_num = 200;
% >>>>>>>>>> Embedded pf for comparison only >>>>>>>>>
if SHOW_FIG == 1
    N_PF = 1000;
    w_PF = ones(1,N_PF) * 1/N_PF;
    x_min_PF = -30*SCALE;
    x_max_PF = 25*SCALE;
    len_PF = x_max_PF - x_min_PF;
    x_PF = x_min_PF:(len_PF/(N_PF-1)):x_max_PF;
end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


for k = 1:K
    if SHOW_FIG == 1
         figure
         fnplt(prior);
        hold on;
        plot_PF(x_PF,w_PF);
    end
    %>>> predict >>>
    cur_knots = SCALE*knots + U;
    cur_knots = sort(cur_knots);
    x_range = min(knots):((max(knots) - min(knots)))/(NN_num - 1):max(knots);
%     cur_x_range_min = min(SCALE*x_range+U) - 3*sqrt(variance);
%     cur_x_range_max = max(SCALE*x_range+U) + 3*sqrt(variance);
    cur_x_range_min = min(x_range/2 + SCALE*25*x_range./(1+x_range.^2) + 8*cos(1.2*k) + U) - 3*sqrt(variance);
    cur_x_range_max = max(x_range/2 + SCALE*25*x_range./(1+x_range.^2) + 8*cos(1.2*k) + U) + 3*sqrt(variance);
    cur_x_range = cur_x_range_min:((cur_x_range_max - cur_x_range_min))/(NN_num - 1):cur_x_range_max;
    
    cur_x_range = sort(cur_x_range);
    P=zeros(length(cur_x_range),length(cur_x_range));
    for j = 1:length(cur_x_range)
        for i = 1:length(x_range)
            %delta =  cur_x_range(j) - (SCALE*x_range(i) + U);
            delta = cur_x_range(j) - (x_range(i)/2 + SCALE*25*x_range(i)/(1+x_range(i)^2) + 8*cos(1.2*k) + U);
            P(i,j) = (1/(sqrt(2*pi*variance))) * exp(-(delta)^2/(2*variance));
        end
    end
    spTran = spap2({M-N+1,M-N+1},{N,N},{x_range,cur_x_range},P);
    knots = spTran.knots{1};
    cur_knots = spTran.knots{2};
    prior = opt_sp(prior,knots,N,NN_num);
    tran_coefs = spTran.coefs;
    
    sp_Basic = spap2(augknt(unique(knots),N),N,x_range,ones(1,length(x_range)));
    
    r=M;
    numOfBasicSpline = r;
%    IntOfBasicSP = zeros(r,r);
%     coefs = sp_Basic.coefs;
%     L = length(coefs);
%     for i = 1:r
%         coef_temp_i = zeros(1,L);
%         coef_temp_i(i) = coefs(i); 
%         sp_temp_i = sp_Basic;
%         sp_temp_i.coefs = coef_temp_i;
%         for j = i:min(i+N,r)
%             coef_temp_j = zeros(1,L);
%             coef_temp_j(j) = coefs(j);
%             sp_temp_j = sp_Basic;
%             sp_temp_j.coefs = coef_temp_j;
%             product = fncmb(sp_temp_i ,'*', sp_temp_j);
%             IntOfBasicSP(i,j) = fnval(fnint(product), product.breaks(end));
%             IntOfBasicSP(j,i) = IntOfBasicSP(i,j);
%         end
%     end
    
    
    IntOfBasicSP = zeros(r,r);
    for i = 1:r
        kknot1 = sp_Basic.knots(i:i+3);%%%%%%%%%%%%%%%%%%%%
        for j = i:min(i+N,r)
            kknot2 = sp_Basic.knots(j:j+3);%%%%%%%%%%%%%%%%%%
            IntOfBasicSP(i,j) = integrationOfTwoSingleSpline(kknot1, kknot2);%%%%%%%%%%%%%%%%%%%
            IntOfBasicSP(j,i) = IntOfBasicSP(i,j);%%%%%%%%%%%%%%%%
        end
    end
    
    
    
    
    
    V = zeros(r,r);
    for i = 1:r
        for j = 1:r
            temp_i(j) = tran_coefs(1,j,i);
        end
        for t = 1:r
            V(i,t) = temp_i * IntOfBasicSP(:,t);
        end
    end
    
    predict = prior;
    pri_coefs = prior.coefs;
    pre_coefs = zeros(1,numOfBasicSpline);
    for i = 1:numOfBasicSpline
        pre_coefs(i) = pri_coefs * V(i,:)';
    end
    
    x_range = min(cur_knots):((max(cur_knots) - min(cur_knots)))/(NN_num - 1):max(cur_knots); 
    temp = spap2(augknt(unique(cur_knots),N),N,x_range,ones(1,length(x_range)));
    
    predict = temp;
    predict.coefs = pre_coefs;
    cump = fnval(fnint(predict),predict.knots(end));
    predict.coefs = predict.coefs/cump;

% >>>>>>>>>>>>>>>>>>>> PF prediction >>>>>>>>>>>>>
    if SHOW_FIG == 1
        for i = 1:N_PF
            x_PF(i) = x_PF(i)/2 + SCALE*25*x_PF(i)/(1+x_PF(i)^2) +8*cos(1.2*k)+ randn * sqrt(variance);
        end
    end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
    if SHOW_FIG == 1
         figure
         fnplt(predict);
         title(num2str(k));
         hold on;
         plot(cur_knots,zeros(1,length(knots)),'v');
         plot_PF(x_PF,w_PF);
    end
    

    %>>>>> update >>>>
    
    x = min(cur_knots):((max(cur_knots) - min(cur_knots)))/(NN_num - 1):max(cur_knots);
    %delta = Z(k) - atan(Yp./(x - 4*k));
    delta = Z(k) - x.^2/20;
    Corr = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
    
    ind1 = find(Corr>1e-5);
    if length(ind1) > 20
        min_ind = max(min(ind1)-1,1);
        max_ind = min(max(ind1)+1,length(x));
        x = x(min_ind:max_ind);
        Corr = Corr(min_ind:max_ind);
    end
    
    
%     Corr
%     if length(Corr) == 0
%         a=1;
%     end
    likelihoodDistribution = spaps(x,Corr,1e-8);
    
    if SHOW_FIG == 1
        figure
        fnplt(likelihoodDistribution);
    end
    
    updatedDistribution = fncmb(predict, '*', likelihoodDistribution);
    cum_pro = fnval(fnint(updatedDistribution), updatedDistribution.breaks(end));
    updatedDistribution.coefs = updatedDistribution.coefs/cum_pro;
    vals = fnval(updatedDistribution, x);
    ind = find(vals > 1e-5);
    min_ind = min(ind);
    max_ind = max(ind);
    x = x(min_ind):((x(max_ind) - x(min_ind)))/(NN_num - 1):x(max_ind);
    vals = fnval(updatedDistribution, x);
    prior = spap2(M-N+1,N,x,vals);
    updatedPDF(k) = prior;
    knots = unique(prior.knots);
    
% >>>>>>>>>>>>> PF updated >>>>>>>>>>>>>>>>>>>>>>>>
    if SHOW_FIG == 1
        for i = 1:N_PF    
            delta_PF = Z(k) - x_PF(i)^2/20;
            w_PF(i) = w_PF(i) * (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta_PF)^2/(2*m_noise_sd^2)));
        end
        if sum(w_PF) == 0
            w_PF = ones(1,N_PF) * 1/N_PF;
        else
            w_PF = w_PF/sum(w_PF);
        end
        N_eff = 1/sum(w_PF.^2);

        if N_eff < 0.8*N_PF    
            w_PF = w_PF/sum(w_PF);
            CDF(1) = 0;
            CDF(2:N_PF+1) = cumsum(w_PF);
            u1 = rand/N_PF;
            i = 0;
            for j = 1:N_PF
                uj = u1 + (j-1)/N_PF;
                while uj > CDF(i+1)
                    i = i+1;
                end
                x_resampled(j) = x_PF(i);
            end   
            x_PF = x_resampled;
            w_PF = ones(1,N_PF) * 1/N_PF;
        end
    end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if SHOW_FIG == 1
         figure
         fnplt(prior);
         title(num2str(k));
         hold on;
         plot(knots,zeros(1,length(knots)),'o');
         plot_PF(x_PF,w_PF);
    end
end
end

function sp = opt_sp(sp,knots,N,NN_num)
    x = min(sp.knots):((max(sp.knots) - min(sp.knots)))/(NN_num - 1):max(sp.knots);
    vals = fnval(sp,x);
    sp = spap2(augknt(unique(knots),N),N,x,vals);
end

function val = integrationOfTwoSingleSpline(knot1, knot2)
val = 0;
a = knot1(1);
b = knot1(2);
c = knot1(3);
d = knot1(4);

if knot1(1) == knot2(1) && knot1(2) == knot2(2) && knot1(3) == knot2(3) && knot1(4) == knot2(4)
    for i = 1:3
        switch i
            case 1
                intg = 0.2 * (b-a)^3/(c-a)^2;
            case 2
                A = a+b-c-d;
                B = 2*(c*d-a*b);
                C = a*b*c - d*a*c + a*b*d - c*b*d;
                intg = 0.2 * A^2*c^5 + 0.5*A*B*c^4 + (1/3)*(2*A*C + B^2)*c^3 + B*C*c^2 + C^2*c...
                    -...
                    (0.2 * A^2*b^5 + 0.5*A*B*b^4 + (1/3)*(2*A*C + B^2)*b^3 + B*C*b^2 + C^2*b);
                intg = intg/(c-b)^2/(d-b)^2/(c-a)^2;
            case 3
                intg = -0.2 * (c-d)^3/(d-b)^2;
        end
        if isnan(intg) == 1
            intg = 0;
        end
        val = val + intg;
    end
elseif knot2(1) == b && knot2(2) == c && knot2(3) == d
    e = knot2(4);
    for i = 1:2
        switch i
            case 1
                A = a+b-c-d;
                B = 2*(c*d-a*b);
                C = a*b*c - d*a*c + a*b*d - c*b*d;
                intg = 0.2 * A*c^5 + 0.25*(B-2*b*A)*c^4 + (1/3)*(b^2*A + C - 2*b*B)*c^3 + 0.5*(B*b^2 - 2*b*C)*c^2 + C*b^2*c...
                    -...
                    (0.2 * A*b^5 + 0.25*(B-2*b*A)*b^4 + (1/3)*(b^2*A + C - 2*b*B)*b^3 + 0.5*(B*b^2 - 2*b*C)*b^2 + C*b^2*b);
                intg = intg/(c-b)/(d-b)/(c-a)/(c-b)/(d-b);
            case 2
                A = b+c-d-e;
                B = 2*(d*e-b*c);
                C = b*c*d - e*b*d + b*c*e - d*c*e;
                intg = 0.2 * A*d^5 + 0.25*(B-2*d*A)*d^4 + (1/3)*(d^2*A + C - 2*d*B)*d^3 + 0.5*(B*d^2 - 2*d*C)*d^2 + C*d^2*d...
                    -...
                    (0.2 * A*c^5 + 0.25*(B-2*d*A)*c^4 + (1/3)*(d^2*A + C - 2*d*B)*c^3 + 0.5*(B*d^2 - 2*d*C)*c^2 + C*d^2*c);
                intg = intg/(d-c)/(e-c)/(d-b)/(d-c)/(d-b);
        end
        if isnan(intg) == 1
            intg = 0;
        end
        val = val + intg;
    end
elseif knot2(1) == c
    e = knot2(3);
    f = knot2(4);
    val = 0.2 * d^5 - 0.5*(c+d)*d^4 + (1/3)*(c^2 + d^2 + 4*c*d)*d^3 - (c*d^2+d*c^2)*d^2 + c^2*d^2*d...
        -...
        (0.2 * c^5 - 0.5*(c+d)*c^4 + (1/3)*(c^2 + d^2 + 4*c*d)*c^3 - (c*d^2+d*c^2)*c^2 + c^2*d^2*c);
    if isnan(val) == 1
            intg = 0;
        end
    val = val/(d-c)/(d-b)/(d-c)/(e-c);
end
end