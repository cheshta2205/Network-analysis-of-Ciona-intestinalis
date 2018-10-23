% Calculating efficiency of neurons in the poor periphery(neurons not present in the rich club)
%Copyright 2018. Cheshta Bhatia & Lav R. Varshney

% Load the adjacency matrix
function A = Poorperiphery_efficiency(varargin)
if (nargin == 0)
    %load the chemical network
    A = datareader('chem','unweighted');
elseif (nargin == 1)
    A = varargin{1};
else
    error('TRIPCOUNT_CHEM: incorrect number of inputs');
end
% To find distances between all neurons in this matrix
for n = 1:1:231
    dist(n,:) = graphshortestpath(A,n);
end
% Choosing distances of all neurons except those belonging to the rich
% club.
dist(:,[80 222 225 142 147 172 173 110 198 143 98 153 163 175 196 197]) = [];
dist([80 222 225 142 147 172 173 110 198 143 98 153 163 175 196 197],:) = [];
%Equating zero values to inf.
dist(dist==0) = inf;
x = 1./dist;
N = sum(x,2);
%Efficiencies of neurons not in the rich club defined by E.
E = N./214;
%Average efficiency or overall efficiency of poor periphery
eff_poor = sum(E)./215;