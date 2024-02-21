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
function [updatedPDF] = SF_movingKnots_2D(prior_sp,prior_knots,K,spOrder,m_noise_sd,Z,variance,nOfKnots,SCALE,DO_PLOT)
N = spOrder;
M1 = nOfKnots(1);%numOfSpline
M2 = nOfKnots(2);
knots{1} = prior_knots{1};
knots{2} = prior_knots{2};
prior = prior_sp;
NN_num1 = 100;
NN_num2 = 10;
F = [1 1;0 1];
Tao = [.5;1];
for k = 1:K
    %>>> predict >>>
    %>>> construct transition function >>>
    x_range{1} = min(knots{1}):((max(knots{1}) - min(knots{1})))/(NN_num1 - 1):max(knots{1});
    x_range{2} = min(knots{2}):((max(knots{2}) - min(knots{2})))/(NN_num2 - 1):max(knots{2});

    cur_x_range_min(1) = min(x_range{1}) + min(x_range{2}) - 3*sqrt(variance)*Tao(1);
    cur_x_range_max(1) = max(x_range{1}) + max(x_range{2}) + 3*sqrt(variance)*Tao(1);
    cur_x_range{1} = cur_x_range_min(1):((cur_x_range_max(1) - cur_x_range_min(1)))/(NN_num1 - 1):cur_x_range_max(1);
    
    cur_x_range_min(2) = min(x_range{2}) - 3*sqrt(variance)*Tao(2);
    cur_x_range_max(2) = max(x_range{2}) + 3*sqrt(variance)*Tao(2);
    cur_x_range{2} = cur_x_range_min(2):((cur_x_range_max(2) - cur_x_range_min(2)))/(NN_num2 - 1):cur_x_range_max(2);

    P=zeros(length(cur_x_range{1}),length(cur_x_range{2}),length(x_range{1}),length(x_range{2}));
    for m = 1:length(cur_x_range{1})
        for n = 1:length(cur_x_range{2})
            for p = 1:length(x_range{1})
                for q = 1:length(x_range{2})
                    delta(1) = cur_x_range{1}(m) - (x_range{1}(p) + x_range{2}(q));
                    delta(2) = cur_x_range{2}(n) - x_range{2}(q);
                    P(p,q,m,n) = (1/(sqrt(2*pi*variance/4))) * exp(-(delta(1))^2/(2*variance/4)) *...
                        (1/(sqrt(2*pi*variance))) * exp(-(delta(2))^2/(2*variance));
                end
            end
        end
    end
    spTran = spap2({M1-N+1,M2-N+1,M1-N+1,M2-N+1},{N,N,N,N},{x_range{1},x_range{2},cur_x_range{1},cur_x_range{2}},P);
    knots{1} = spTran.knots{1};
    knots{2} = spTran.knots{2};
    cur_knots{1} = spTran.knots{3};
    cur_knots{2} = spTran.knots{4};
    prior = opt_sp(prior,knots,N,NN_num1,NN_num2);
    tran_coefs = spTran.coefs;
    
%     sp_Basic = spap2({augknt(unique(knots(1,:)),N),augknt(unique(knots(2,:)),N)},[N N],{x_range(1,:),x_range(2,:)}, ...
%         ones(length(x_range(1,:)),length(x_range(2,:))));
    
    sp_Basic1 = spap2(augknt(unique(knots{1}),N),N,x_range{1}, ones(1,length(x_range{1})));    
    r=M1;
    IntOfBasicSP1 = zeros(r,r);
    coefs = sp_Basic1.coefs;
    L = length(coefs);
    
    for i = 1:r
        coef_temp_i = zeros(1,L);
        coef_temp_i(i) = coefs(i); 
        sp_temp_i = sp_Basic1;
        sp_temp_i.coefs = coef_temp_i;
        for j = i:min(i+N,r)
            coef_temp_j = zeros(1,L);
            coef_temp_j(j) = coefs(j);
            sp_temp_j = sp_Basic1;
            sp_temp_j.coefs = coef_temp_j;
            product = fncmb(sp_temp_i ,'*', sp_temp_j);
            IntOfBasicSP1(i,j) = fnval(fnint(product), product.breaks(end));
            IntOfBasicSP1(j,i) = IntOfBasicSP1(i,j);
        end
    end
    
    sp_Basic2 = spap2(augknt(unique(knots{2}),N),N,x_range{2}, ones(1,length(x_range{2})));
    r=M2;
    IntOfBasicSP2 = zeros(r,r);
    coefs = sp_Basic2.coefs;
    L = length(coefs);
    
    for i = 1:r
        coef_temp_i = zeros(1,L);
        coef_temp_i(i) = coefs(i); 
        sp_temp_i = sp_Basic2;
        sp_temp_i.coefs = coef_temp_i;
        for j = i:min(i+N,r)
            coef_temp_j = zeros(1,L);
            coef_temp_j(j) = coefs(j);
            sp_temp_j = sp_Basic2;
            sp_temp_j.coefs = coef_temp_j;
            product = fncmb(sp_temp_i ,'*', sp_temp_j);
            IntOfBasicSP2(i,j) = fnval(fnint(product), product.breaks(end));
            IntOfBasicSP2(j,i) = IntOfBasicSP2(i,j);
        end
    end
    
    IntOfBasicSP = zeros(M1,M2,M1,M2);
    
    for i = 1:M1
        for j = 1:M2
            for p = 1:M1
                for q = 1:M2
                    IntOfBasicSP(i,j,p,q) = IntOfBasicSP1(i,p) * IntOfBasicSP2(j,q);
                end
            end
        end
    end
    
    pre_coefs = zeros(M1,M2);
    D_ij = zeros(M1,M2);
    pri_coefs = zeros(M1,M2);
    for i = 1:M1
        for j = 1:M2
            pri_coefs(i,j) = prior.coefs(1,i,j);
        end
    end
    
    for m = 1:M1
        for n = 1:M2
            for p = 1:M1
                for q = 1:M2
                    B_pq(p,q) = tran_coefs(1,p,q,m,n);
                end
            end
            for i = 1:M1
                for j = 1:M2
                    for p = 1:M1
                        for q = 1:M2
                            A_pq(p,q) = IntOfBasicSP(i,j,p,q);
                        end
                    end
                    D_ij(i,j) = trace(B_pq'*A_pq);
                end
            end
            pre_coefs(m,n) = trace(pri_coefs'*D_ij);
        end
    end
    
    x_range{1} = min(cur_knots{1}):((max(cur_knots{1}) - min(cur_knots{1})))/(NN_num1 - 1):max(cur_knots{1}); 
    x_range{2} = min(cur_knots{2}):((max(cur_knots{2}) - min(cur_knots{2})))/(NN_num2 - 1):max(cur_knots{2});
    temp = spap2({augknt(unique(cur_knots{1}),N),augknt(unique(cur_knots{2}),N)},[N N],{x_range{1},x_range{2}},...
        ones(length(x_range{1}),length(x_range{2})));
    
    predict = temp;
    predict.coefs = pre_coefs;
    %cump = fnval(fnint(predict),{predict.knots{1}(end),predict.knots{2}(end)});
    %predict.coefs = predict.coefs/cump;
    
    %>>>>> update >>>>
    x{1} = min(cur_knots{1}):((max(cur_knots{1}) - min(cur_knots{1})))/(NN_num1 - 1):max(cur_knots{1});
    x{2} = min(cur_knots{2}):((max(cur_knots{2}) - min(cur_knots{2})))/(NN_num2 - 1):max(cur_knots{2});
    for t = 1:length(x{1})
        for r = 1:length(x{2})
            delta = Z(k) - atan2(20,x{1}(t) - 4*k);
            Corr(t,r) = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
        end
    end
    %likelihoodDistribution = spaps({x{1},x{2}},Corr,1e-6);
    
    vals = zeros(length(x{1}),length(x{2}));
    for i = 1:length(x{1})
        for j = 1:length(x{2})
            vals(i,j) = fnval(predict,{x{1}(i),x{2}(j)}) * Corr(i,j);
        end
    end
    
    updatedDistribution = spap2({M1-N+1,M2-N+1},[N N],{x{1},x{2}},vals);
    
    knots{1} = updatedDistribution.knots{1};
    knots{2} = updatedDistribution.knots{2};
    sp_Basic1 = spap2(augknt(unique(knots{1}),N),N,x{1}, ones(1,length(x{1})));
    sp_Basic2 = spap2(augknt(unique(knots{2}),N),N,x{2}, ones(1,length(x{2})));
    cump = 0;
    for i = 1:M1
        coef_temp_i = zeros(1,M1);
        coef_temp_i(i) = sp_Basic1.coefs(i); 
        sp_temp_i = sp_Basic1;
        sp_temp_i.coefs = coef_temp_i;
        for j = 1:M2
            coef_temp_j = zeros(1,M2);
            coef_temp_j(j) = sp_Basic2.coefs(j);
            sp_temp_j = sp_Basic2;
            sp_temp_j.coefs = coef_temp_j;
            cump = cump + updatedDistribution.coefs(1,i,j) * fnval(fnint(sp_temp_i),max(sp_temp_i.knots)) * fnval(fnint(sp_temp_j),max(sp_temp_j.knots)); 
        end
    end
    updatedDistribution.coefs = updatedDistribution.coefs/cump;
    
    %>>> cut the zero probability region >>>
    ind = find(vals > 1e-3);
    AREA = zeros(length(x{1}),length(x{2}));
    AREA(ind) = 1;
    x1_range = sum(AREA,2);
    x1_min_ind = min(find(x1_range > 0));
    x1_min = x{1}(x1_min_ind);
    x1_max_ind = max(find(x1_range > 0));
    x1_max = x{1}(x1_max_ind);
    
    x2_range = sum(AREA,1);
    x2_min_ind = min(find(x2_range > 0));
    x2_min = x{2}(x2_min_ind);
    x2_max_ind = max(find(x2_range > 0));
    x2_max = x{2}(x2_max_ind);
    
    x{1} = x1_min:((x1_max - x1_min))/((x1_max_ind - x1_min_ind + 1) - 1):x1_max;
    x{2} = x2_min:((x2_max - x2_min))/((x2_max_ind - x2_min_ind + 1) - 1):x2_max;
    vals = vals(x1_min_ind:x1_max_ind,x2_min_ind:x2_max_ind);
    updatedDistribution = spap2({M1-N+1,M2-N+1},[N N],{x{1},x{2}},vals);
    
    knots{1} = updatedDistribution.knots{1};
    knots{2} = updatedDistribution.knots{2};
    sp_Basic1 = spap2(augknt(unique(knots{1}),N),N,x{1}, ones(1,length(x{1})));
    sp_Basic2 = spap2(augknt(unique(knots{2}),N),N,x{2}, ones(1,length(x{2})));
    cump = 0;
    for i = 1:M1
        coef_temp_i = zeros(1,M1);
        coef_temp_i(i) = sp_Basic1.coefs(i); 
        sp_temp_i = sp_Basic1;
        sp_temp_i.coefs = coef_temp_i;
        for j = 1:M2
            coef_temp_j = zeros(1,M2);
            coef_temp_j(j) = sp_Basic2.coefs(j);
            sp_temp_j = sp_Basic2;
            sp_temp_j.coefs = coef_temp_j;
            cump = cump + updatedDistribution.coefs(1,i,j) * fnval(fnint(sp_temp_i),max(sp_temp_i.knots)) * fnval(fnint(sp_temp_j),max(sp_temp_j.knots)); 
        end
    end
    updatedDistribution.coefs = updatedDistribution.coefs/cump;  
    prior = updatedDistribution;
    
    if DO_PLOT == 1
        figure;
        fnplt(prior);
    end
    
    k
    updatedPDF(k) = prior;
    knots{1} = unique(prior.knots{1});
    knots{2} = unique(prior.knots{2});
end
end

function sp = opt_sp(sp,knots,N,NN_num1,NN_num2)
    x{1} = min(sp.knots{1}):((max(sp.knots{1}) - min(sp.knots{1})))/(NN_num1 - 1):max(sp.knots{1});
    x{2} = min(sp.knots{2}):((max(sp.knots{2}) - min(sp.knots{2})))/(NN_num2 - 1):max(sp.knots{2});
    vals = zeros(length(x{1}),length(x{2}));
    for i = 1:length(x{1})
        for j = 1:length(x{2})
            vals(i,j) = fnval(sp,{x{1}(i),x{2}(j)});
        end
    end
    sp = spap2({augknt(unique(knots{1}),N),augknt(unique(knots{2}),N)},[N N],{x{1} x{2}},vals);
end