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
function CRLB = crlb(X1,m_noise_sd,F,Tao,variance,J_0,Yp)
Q = Tao * Tao' * variance;
q = Q(1,1);
f = F(1,1);
for k = 1:length(X1(1,:))
    if k == 1
        J(k) = (f * J_0^(-1)*f + q)^(-1) +  1/(m_noise_sd^2 * ((Yp^2+(X1(1,k) - 4*k)^2) / (-Yp))^2);
    else
        J(k) = (f * J(k-1)^(-1)*f + q)^(-1) +  1/(m_noise_sd^2 * ((Yp^2+(X1(1,k) - 4*k)^2) / (-Yp))^2);
    end
    CRLB(k) = 1/J(k);
end