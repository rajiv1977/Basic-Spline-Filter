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
clear all;
tic;
% Initialize the grid and initial and final points
nx1 = 30; nx2 = 30;
x1l = 0; x1u = 100;
x2l = 0; x2u = 100;
x1 = linspace(x1l,x1u,nx1+1);
x2 = linspace(x2l,x2u,nx2+1);
limsf1 = 1:nx1+1; limsf2 = 1:nx2+1;
% Initalize other variables
T = 10;
sigma0=1;
sigma = sigma0^2*[T^3/3 T^2/2;T^2/2 T];
invSig = inv(sigma);
detSig = det(sigma);
expF = [1 T; 0 1];
n = size (expF, 1); %2
gausThresh = 10;
small = 0; subs = []; vals = [];
% Iterate through all possible initial
% and final positions and calculate
% the values of exponent and out
% if exponent > gausThresh.
for i1 = 1:nx1+1
for i2 = 1:nx2+1
for f1 = limsf1
for f2 = limsf2
% Initial and final position
xi = [x1(i1) x2(i2)]';
xf = [x1(f1) x2(f2)]';
exponent = 0.5 * (xf - expF * xi)'...
* invSig * (xf - expF * xi);
if exponent > gausThresh
small = small + 1;
else
out = 1 / (sqrt((2 * pi)^n * detSig))...
* exp(-exponent);
subs = [subs; i1 i2 f1 f2];
vals = [vals; out];
end
end
end
end
end
small = ((nx1+1)^2*(nx2+1)^2)-length(subs);
siz=[nx1+1 nx2+1 nx1+1 nx2+1];
ttp0=sptensor(subs,vals,siz);
toc;
