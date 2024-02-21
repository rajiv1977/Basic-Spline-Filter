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
function [xPF wPF] = PF_2D(cdY1,NTsteps,variance,m_noise_sd,N,SCALE)
w = ones(1,N) * 1/N;
x1_min = 0;
x1_max = 100;
x2_min = 0;
x2_max = 2;
F = [1 1;0 1];
Tao = [.5;1];

x = zeros(2,N);
for i = 1:N
    x(1,i) = rand * (x1_max - x1_min) + x1_min;
    x(2,i) = rand * (x2_max - x2_min) + x2_min;
end


for k = 1:NTsteps
    for i = 1:N
        x(:,i) = F*x(:,i) + randn * sqrt(variance) * Tao;
        delta = cdY1(k) - atan2(20,(x(1,i) - 4*k));
        w(i) = w(i) * (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta)^2/(2*m_noise_sd^2)));
    end
    if sum(w) == 0
        w = ones(1,N) * 1/N;
    else
        w = w/sum(w);
    end
   
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
            x_resampled(:,j) = x(:,i);
        end   
        x = x_resampled;
        xPF(k,:,:) = x;
        w = ones(1,N) * 1/N;
    else
        xPF(k,:,:) = x;
    end
    wPF(k,:)=w;   
end