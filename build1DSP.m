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
function coeff = build1DSP(knots,f,x,nOfSpline)
b = f';
n = nOfSpline;
A = zeros(length(f),n);
for i = 1:length(f)
    for j = 1:n
        A(i,j) = evaluateSingleSpline(knots(j:j+3),x(i));
    end
end
coeff = pinv(A)*b;
end

function y = evaluateSingleSpline(knot1,x)
y=0;
a = knot1(1);
b = knot1(2);
c = knot1(3);
d = knot1(4);
sum1 = (x-a)^2/(b-a)/(c-a);
if isnan(sum1) || isinf(sum1)
    sum1 = 0;
end
A = a+b-c-d;
B = 2*(c*d-a*b);
C = a*b*c - d*a*c + a*b*d - c*b*d;
sum2 = A*x^2 + B*x + C;
if isnan(sum2)|| isinf(sum2)
    sum2 = 0;
end
sum3 = (x-d)^2/(d-c)/(d-b);
if isnan(sum3)|| isinf(sum3)
    sum3 = 0;
end
y = sum1 + sum2 + sum3;
end