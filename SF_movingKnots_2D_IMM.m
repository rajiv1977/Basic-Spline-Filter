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
function [updatedPDF modeProb] = SF_movingKnots_2D_IMM(prior_sp1,prior_knots1,prior_sp2,prior_knots2,K,spOrder,m_noise_sd,Z,variance,nOfKnots,SCALE,DO_PLOT,F1,Tao1,G1,F2,Tao2,G2)
N = spOrder;
M1 = nOfKnots(1);%numOfSpline
M2 = nOfKnots(2);
knots1{1} = prior_knots1{1};
knots1{2} = prior_knots1{2};

knots2{1} = prior_knots2{1};
knots2{2} = prior_knots2{2};

prior1 = prior_sp1;
prior2 = prior_sp2;
NN_num1 = 100;
NN_num2 = 20;

Markov = [0.9 0.1;0.2 0.8];
u_prior1 = 0.5;
u_prior2 = 0.5;
for k = 1:K
    %>>> predict >>>
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
    
    prior1_0 = merge_sp(prior1,prior2,NN_num1,NN_num2,u11_predict,u21_predict,N,M1,M2);
    prior2_0 = merge_sp(prior1,prior2,NN_num1,NN_num2,u12_predict,u22_predict,N,M1,M2);
    
    knots1{1} = prior1_0.knots{1};
    knots1{2} = prior1_0.knots{2};
    knots2{1} = prior2_0.knots{1};
    knots2{2} = prior2_0.knots{2};
    
    [x_range1 cur_x_range1] = getXandCurX(knots1,NN_num1,NN_num2,variance,F1,Tao1,G1,k);
    [x_range2 cur_x_range2] = getXandCurX(knots2,NN_num1,NN_num2,variance,F2,Tao2,G2,k);

    [predict1 knots1 cur_knots1] = prediction(knots1,N,M1,M2,prior1_0,x_range1,cur_x_range1,F1,variance,NN_num1,NN_num2,Tao1,G1,k);
    [predict2 knots2 cur_knots2] = prediction(knots2,N,M1,M2,prior2_0,x_range2,cur_x_range2,F2,variance,NN_num1,NN_num2,Tao2,G2,k);
    
    %>>>>> update >>>>
    [prior1 likelihood1 knots1] = update(cur_knots1,Z,predict1,knots1,NN_num1,NN_num2,k,m_noise_sd,M1,M2,N);
    [prior2 likelihood2 knots2] = update(cur_knots2,Z,predict2,knots2,NN_num1,NN_num2,k,m_noise_sd,M1,M2,N);
    
    knots1{1} = prior1.knots{1};
    knots1{2} = prior1.knots{2};
    knots2{1} = prior2.knots{1};
    knots2{2} = prior2.knots{2};
    
    u_prior1 = (u11_predict + u21_predict) * likelihood1;
    u_prior2 = (u12_predict + u22_predict) * likelihood2;
    sum_u = u_prior1 + u_prior2;
    u_prior1 = u_prior1/sum_u;
    u_prior2 = u_prior2/sum_u;
    updatedPDF(k) = merge_sp(prior1,prior2,NN_num1,NN_num2,u_prior1,u_prior2,N,M1,M2);
    modeProb(1,k) = u_prior1;
    modeProb(2,k) = u_prior2;
    figure
    fnplt(updatedPDF(k));
    modeProb
    k
end
end






function [x_range cur_x_range] = getXandCurX(knots,NN_num1,NN_num2,variance,F,Tao,G,k)
    x_range{1} = min(knots{1}):((max(knots{1}) - min(knots{1})))/(NN_num1 - 1):max(knots{1});
    x_range{2} = min(knots{2}):((max(knots{2}) - min(knots{2})))/(NN_num2 - 1):max(knots{2});

    cur_x_range_min(1) = F(1,:) * [min(x_range{1});min(x_range{2})] + G(1) * (k-1) - 3*sqrt(variance)*Tao(1);
    cur_x_range_max(1) = F(1,:) * [max(x_range{1});max(x_range{2})] + G(1) * (k-1) + 3*sqrt(variance)*Tao(1);
    cur_x_range{1} = cur_x_range_min(1):((cur_x_range_max(1) - cur_x_range_min(1)))/(NN_num1 - 1):cur_x_range_max(1);
    
    cur_x_range_min(2) = F(2,:) * [min(x_range{1});min(x_range{2})] + G(2) * (k-1)- 3*sqrt(variance)*Tao(2);
    cur_x_range_max(2) = F(2,:) * [max(x_range{1});max(x_range{2})] + G(2) * (k-1) + 3*sqrt(variance)*Tao(2);
    cur_x_range{2} = cur_x_range_min(2):((cur_x_range_max(2) - cur_x_range_min(2)))/(NN_num2 - 1):cur_x_range_max(2);
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

function [predict knots cur_knots] = prediction(knots,N,M1,M2,prior,x_range,cur_x_range,F,variance,NN_num1,NN_num2,Tao,G,k)
    variance1 = variance * Tao(1)^2;
    variance2 = variance * Tao(2)^2;
    %S = variance * Tao * Tao';
    %inv_S = inv(S);
    P=zeros(length(cur_x_range{1}),length(cur_x_range{2}),length(x_range{1}),length(x_range{2}));
    for m = 1:length(cur_x_range{1})
        for n = 1:length(cur_x_range{2})
            for p = 1:length(x_range{1})
                for q = 1:length(x_range{2})
                    delta = [cur_x_range{1}(m);cur_x_range{2}(n)] - (F*[x_range{1}(p); x_range{2}(q)] + G*(k-1));
                    P(p,q,m,n) = (1/(sqrt(2*pi*variance1))) * exp(-(delta(1))^2/(2*variance1)) *...
                        (1/(sqrt(2*pi*variance2))) * exp(-(delta(2))^2/(2*variance2));
                    %P(p,q,m,n) = 1/sqrt(det(2*pi*S))*exp(-0.5*(delta)'*inv_S*(delta)); 
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
end


function [prior likelihood knots] = update(cur_knots,Z,predict,knots,NN_num1,NN_num2,k,m_noise_sd,M1,M2,N)
    x{1} = min(cur_knots{1}):((max(cur_knots{1}) - min(cur_knots{1})))/(NN_num1 - 1):max(cur_knots{1});
    x{2} = min(cur_knots{2}):((max(cur_knots{2}) - min(cur_knots{2})))/(NN_num2 - 1):max(cur_knots{2});
    for t = 1:length(x{1})
        for r = 1:length(x{2})
            delta = Z(k) - atan(20/(x{1}(t) - 4*k));
            Corr(t,r) = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
        end
    end
    
    vals = zeros(length(x{1}),length(x{2}));
    for i = 1:length(x{1})
        for j = 1:length(x{2})
            vals(i,j) = fnval(predict,{x{1}(i),x{2}(j)}) * Corr(i,j);
        end
    end
    
    ind = find(vals < 0);
    vals(ind) = 0;
    
    
    
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
    vals_backup = vals;
    x_backup = x;
    magic = 5e-3;
    ind = find(vals > magic);
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
    
    counter = 1;
    while ((length(find(x1_range > 0)) < 2) || (length(find(x2_range > 0)) < 2)) && (counter < 10)
        counter = counter + 1;
        vals = vals_backup;
        x = x_backup;
        magic = magic/counter;
        ind = find(vals > magic);
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
    end
    if ((length(find(x1_range > 0)) < 2) || (length(find(x2_range > 0)) < 2))
        vals = vals_backup;
        x = x_backup;
    else
        x{1} = x1_min:((x1_max - x1_min))/((x1_max_ind - x1_min_ind + 1) - 1):x1_max;
        x{2} = x2_min:((x2_max - x2_min))/((x2_max_ind - x2_min_ind + 1) - 1):x2_max;
        vals = vals(x1_min_ind:x1_max_ind,x2_min_ind:x2_max_ind);
    end
    
    counter = 1;
    while ((length(x{1}) < N - 1 + M1-N+1) || (length(x{2}) < N - 1 + M2-N+1)) && (counter < 10)
        counter = counter + 1;
        vals = vals_backup;
        x = x_backup;
        magic = magic/counter;
        ind = find(vals > magic);
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
    end
    if ((length(x{1}) < N - 1 + M1-N+1) || (length(x{2}) < N - 1 + M2-N+1))
        vals = vals_backup;
        x = x_backup;
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
    
    likelihood = cump;
    prior = updatedDistribution;
    
    knots{1} = unique(prior.knots{1});
    knots{2} = unique(prior.knots{2});
end

function sp = merge_sp(sp1,sp2,NN_num1,NN_num2,w1,w2,N,M1,M2)
    sum_w = w1 + w2;
    w1 = w1/sum_w;
    w2 = w2/sum_w;
    cur_knots{1} = [sp1.knots{1}, sp2.knots{1}];
    cur_knots{2} = [sp1.knots{2}, sp2.knots{2}];
    x_range{1} = min(cur_knots{1}):((max(cur_knots{1}) - min(cur_knots{1})))/(NN_num1 - 1):max(cur_knots{1}); 
    x_range{2} = min(cur_knots{2}):((max(cur_knots{2}) - min(cur_knots{2})))/(NN_num2 - 1):max(cur_knots{2});
    for i = 1:length(x_range{1})
        for j = 1:length(x_range{2})
            val(i,j) = w1*fnval(sp1,{x_range{1}(i),x_range{2}(j)}) + w2*fnval(sp2,{x_range{1}(i),x_range{2}(j)});
        end
    end
    %sp = spap2({augknt(unique(cur_knots{1}),N),augknt(unique(cur_knots{2}),N)},[N N],{x_range{1},x_range{2}},val);
    sp = spap2({M1-N+1,M2-N+1},{N,N},{x_range{1},x_range{2}},val);
end