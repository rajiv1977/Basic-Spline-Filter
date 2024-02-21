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
function [X1lims u] = numerical(cdY1,NTsteps,variance,m_noise_sd,dx1)
X1lims = 65:dx1:160;
lx1 = length(X1lims);
for k = 1:NTsteps
    for ii=1:lx1
        x=X1lims(ii);
        for jj=1:lx1
            y=X1lims(jj);
            %delta = y - (x/2 + 25*x/(1+x^2) + 8*cos(1.2*k));
            delta = y - (2*x + 1);
            P(k,jj,ii)= (1/(sqrt(2*pi*variance))) * exp(-((delta)^2/(2*variance)));
        end
    end
end
%Original from Balaji paper
u(1:NTsteps,1:lx1)=1;
u0=u(1,:);
for kk=1:NTsteps
    for i = 1:length(u0)
        for j = 1:length(u0)
            PP(i,j) = P(kk,i,j);
        end
    end
    if kk == 1
        upred=PP*u0.'*dx1; %Prediction
        x=X1lims;
        %delta = cdY1(kk) - x.^2/20;
        delta = cdY1(kk) - atan2(20,(x - 4*kk));
        Corr = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
        ucorr=Corr.*upred.'; %Correction
        u(kk,:)=ucorr/(sum(ucorr*dx1)); %Normalized
    else
        upred=PP*u(kk-1,:).'*dx1; %Prediction
        x=X1lims;
        %delta = cdY1(kk) - x.^2/20;
        delta = cdY1(kk) - atan2(20,(x - 4*kk));
        Corr = (1/(sqrt(2*pi*m_noise_sd^2))) * exp(-((delta).^2/(2*m_noise_sd^2)));
        ucorr=Corr.*upred.'; %Correction
        u(kk,:)=ucorr/(sum(ucorr*dx1)); %Normalized
    end
end