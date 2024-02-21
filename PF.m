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
function [xPF wPF] = PF(cdY1,NTsteps,variance,m_noise_sd,N,SCALE,U,Yp, x_min, x_max)
w = ones(1,N) * 1/N;
len = x_max - x_min;
x = x_min:(len/(N-1)):x_max;
for k = 1:NTsteps
    for i = 1:N
        x(i) = x(i)/2 + SCALE*25*x(i)/(1 + x(i)^2) + 8 *cos(1.2 * k) + randn * sqrt(variance) + U;
        %x(i) = SCALE*x(i) + U+ randn * sqrt(variance);
        delta = cdY1(k) - x(i)^2/20;
        %delta = cdY1(k) - atan(Yp/(x(i) - 4*k));
        w(i) = w(i) * (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta)^2/(2*m_noise_sd^2)));
    end
    if sum(w) == 0
        w = ones(1,N) * 1/N;
    else
        w = w/sum(w);
    end
    
    
%         figure
%     for kk = 1:1:NTsteps
%         hist(x,length(unique(x)));
%         title(num2str(k));
%     end
    
    
    N_eff = 1/sum(w.^2);
    
    if N_eff < 0.8*N    
        w = w/sum(w);
        CDF(1) = 0;
        CDF(2:N+1) = cumsum(w);
        u1 = rand/N;
        i = 0;
        for j = 1:N
            uj = u1 + (j-1)/N;
            while uj > CDF(i+1)
                i = i+1;
            end
            x_resampled(j) = x(i);
        end   
        x = x_resampled;
        xPF(k,:) = x;
        w = ones(1,N) * 1/N;
    else
        xPF(k,:) = x;
    end
    wPF(k,:)=w;
%     figure
%     for kk = 1:1:NTsteps
%         hist(x,length(unique(x)));
%         title(num2str(k));
%     end
    
    
    
end