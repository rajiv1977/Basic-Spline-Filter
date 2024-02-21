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
function MovingKonts()
clear all;
close all;
clc;

N = 3;

X0 = 1;
K = 3;
X_true(1) = X0;
sigma_x = 0.1;
for i = 2:K
    noise = randn * sigma_x;
    X_true(i) = 1 * X_true(i-1) +5+ noise; 
end
figure(1);
plot(X_true);

knots = 0:.1:2;
t(1:length(knots)) = 0;
for i = 1:length(t)
    t(i) = (1/(sqrt(2*pi)*sigma_x)) * exp(-(knots(i) - X0)^2/2/sigma_x/sigma_x);
end
prior = spapi(N, knots, t);
fnplt(prior);

for k = 1:K
    cur_knots = 1 .* knots + 5;
    [x, y] = ndgrid(knots,cur_knots);
    P = (1/(sqrt(2*pi)*sigma_x)) * exp(-((y-1*x - 5).^2./(2*sigma_x^2)));
    spTran = spapi({N,N},{knots,cur_knots},P);
    tran_coefs = spTran.coefs;
    
    sp_Basic = spapi(N,knots,ones(1,length(knots)));

    r = length(sp_Basic.coefs);
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
    numOfBasicSpline = length(pri_coefs);
    pre_coefs = zeros(1,numOfBasicSpline);
    for i = 1:numOfBasicSpline
        for j = 1:numOfBasicSpline
            temp_t = pri_coefs * IntOfBasicSP(j,:)';
            pre_coefs(i) = pre_coefs(i) + temp_t * tran_coefs(1,i,j);
        end
        pre_coefs(i) = pri_coefs * V(i,:)';
    end
    
    temp = spapi(N, knots, ones(1,length(knots)));
    predict = temp;
    predict.coefs = pre_coefs;
    
    knots = cur_knots;
    prior = predict;
    updatedPDF(k) = prior;
end
    
for kk = 1:K
figure
fnplt(updatedPDF(kk));
hold on
axis([0 20 0 1]);
end
%Conditional Mean

end