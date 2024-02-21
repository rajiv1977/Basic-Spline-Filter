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
function [X1lims V dx1 IntOfBasicSP sp_Prior] = SP_preconstruction(N,variance,spOrder,NTsteps)
dx1=50/(N-1);
X1lims = -25:dx1:25;

lx1=length(X1lims);

for k = 1:NTsteps
    for ii=1:lx1
        x=X1lims(ii);
        for jj=1:lx1
            y=X1lims(jj);
            delta = y - (x/2 + 25*x/(1+x^2) + 8*cos(1.2*k));
            P(k,jj,ii)= (1/(sqrt(2*pi*variance))) * exp(-((delta)^2/(2*variance)));
        end
    end
    x = X1lims; 
    y = x; 
    [xx,yy] = ndgrid(x,y);
    z = P(k,:,:);
    sp_Transition = spap2({length(x) - spOrder + 1, length(x) - spOrder + 1},[spOrder spOrder],{x,y},z);
    tran_coefs(k,:,:) = sp_Transition.coefs;
end

%Here we are essentially trying to multiply the Transition probabilty
%spline with the prior.



t(1:lx1) = 1;
sp_Basic = spap2(length(x) - spOrder + 1, spOrder, X1lims, t(1,:));

r = length(sp_Basic.coefs);
IntOfBasicSP = zeros(r,r);
coefs = sp_Basic.coefs;
L = length(coefs);

for i = 1:r
    coef_temp_i = zeros(1,L);
    coef_temp_i(i) = coefs(i); 
    sp_temp_i = sp_Basic;
    sp_temp_i.coefs = coef_temp_i;
    for j = i:min(i+spOrder,r)
        coef_temp_j = zeros(1,L);
        coef_temp_j(j) = coefs(j);
        sp_temp_j = sp_Basic;
        sp_temp_j.coefs = coef_temp_j;
        product = fncmb(sp_temp_i ,'*', sp_temp_j);
        IntOfBasicSP(i,j) = fnval(fnint(product), product.breaks(end));
        IntOfBasicSP(j,i) = IntOfBasicSP(i,j);
    end
end
V = zeros(NTsteps,r,r);
for k = 1:NTsteps
    for i = 1:r
        for j = 1:r
            temp_i(j) = tran_coefs(k,i,j);
        end
        for t = 1:r
            V(k,i,t) = temp_i * IntOfBasicSP(:,t);
        end
    end
end
t(1:length(X1lims)) = 1;
sp_Prior = spap2(length(X1lims) - spOrder + 1, spOrder, X1lims, t(1,:));

