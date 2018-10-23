%Superfamily analysis of C.elegans and Ciona intestinalis.
%Finding Z normalized values for C.elegans and Ciona intestinalis and
%grouping them into same superfamily if they show similar triad
%significance profiles.

%Copyright 2018. Cheshta Bhatia & Lav R. Varshney

% Load matrix for C.elegans
%adjacency matrix
function A = superfamily_analysis(varargin)
if (nargin == 0)
    %load the chemical network
    A = datareaderk('chem','unweighted');
elseif (nargin == 1)
    A = varargin{1};
else
    error('TRIPCOUNT_CHEM: incorrect number of inputs');
end

D = (1-A).*(1-A)'; D = D - D.*eye(size(D));
U = A.*(1-A)';
B = A.*A';
tripletcount(1) = 1/2*sum(sum(D.*(U'*U)));
tripletcount(2) = 1/2*sum(sum(D.*(U*U')));
tripletcount(3) = sum(sum(D.*(U^2)));
tripletcount(4) = sum(sum(D.*(U*B)));
tripletcount(5) = sum(sum(D.*(U'*B)));
tripletcount(6) = 1/2*sum(sum(D.*(B^2)));
tripletcount(7) = sum(sum(U.*(U^2)));
tripletcount(8) = 1/3*sum(sum(U'.*(U^2)));
tripletcount(9) = 1/2*sum(sum(B.*(U'*U)));
tripletcount(10) = 1/2*sum(sum(B.*(U*U')));
tripletcount(11) = sum(sum(B.*(U^2)));
tripletcount(12) = sum(sum(U.*(B^2)));
tripletcount(13) = 1/6*sum(sum(B.*(B^2)));

%triplet counts for chemical and for related random graphs

%chemical network
tripletcount = tripCount_chem;
%randomized with matched degree and matched doublets

for ii = 1:1000
    tripletcountDDM(ii,:) = tripCount_chem(GenRandomMatrix3(datareaderk('chem','unweighted')));
%    tripletcountDDM(ii,:) = tripCount_chem(CreateRandomMatrix6(double(datareader('chem','unweighted'))));
end
%Finding the mean, standard deviation and Z values.
mean_Nrand_elegans = mean(tripletcountDDM);
std_Nrand_elegans = std(tripletcountDDM);
Zi_elegans = (tripletcount-mean_Nrand_elegans) ./ std_Nrand_elegans;
Znorm_elegans = Zi_elegans ./ (sumsqr(Zi_elegans))^0.5;
x = 1:13;
y1 = Znorm_elegans;
save('TSP_elegans')

