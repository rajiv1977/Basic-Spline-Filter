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
function x = findMean(coeff,knots)
x = 0;
for i = 1:length(coeff)
    a = knots(i);
    b = knots(i+1);
    c = knots(i+2);
    d = knots(i+3);
    A = a+b-c-d;
    B = 2*(c*d-a*b);
    C = a*b*c - d*a*c + a*b*d - c*b*d;
    
    intg1 = 0.25*b^4 - (2/3)*a*b^3 + 0.5*a^2*b^2 ...
    -...
        (0.25*a^4 - (2/3)*a*a^3 + 0.5*a^2*a^2);
    intg1 =intg1/(b-a)/(c-a);
    if isnan(intg1)
        intg1 = 0;
    end
    intg2 = 0.25 * A*c^4 + (1/3)*B*c^3 + 0.5*C*c^2 ...
        -...
        (0.25 * A*b^4 + (1/3)*B*b^3 + 0.5*C*b^2);
    intg2 =intg2/(c-b)/(d-b)/(c-a);
    if isnan(intg2)
        intg2 = 0;
    end
    intg3 = 0.25*d^4 - (2/3)*d*d^3 + 0.5*d^2*d^2 ...
    -...
        (0.25*c^4 - (2/3)*d*c^3 + 0.5*d^2*c^2);
    intg3 =intg3/(d-c)/(d-b);
    if isnan(intg3)
        intg3 = 0;
    end
    x = x + (intg1 + intg2 + intg3) * coeff(i);
end
end