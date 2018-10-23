% Calculating efficiency of neurons in Rich Club.

%Copyright 2018. Cheshta Bhatia & Lav R. Varshney

%Load the adjacency matrix A.

function A = RichClub_efficiency(varargin)
%adjacency matrix
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
% Choosing distances of only those neurons that belong to the rich club
y = dist(:,[80 222 225 142 147 172 173 110 198 143 98 153 163 175 196 197])
x = y([80 222 225 142 147 172 173 110 198 143 98 153 163 175 196 197],:)
%Equating zero values to inf.
x(x==0) = inf;
p = 1./x;
N = sum(p,2);
%Efficiencies of neurons of the rich club defined by E.
E = N./15;
%Average efficiency or overall efficiency of the Rich club
eff_rich = sum(E)./16
