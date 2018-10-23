%Superfamily analysis of C.elegans and Ciona intestinalis.
%Finding Z normalized values for C.elegans and Ciona intestinalis and
%grouping them into same superfamily if they show similar triad
%significance profiles.

%Copyright 2018. Cheshta Bhatia & Lav R. Varshney

% Load matrix for Ciona intestinalis
%adjacency matrix
function A = superfamily_analysis(varargin)
if (nargin == 0)
    %load the chemical network
    A = datareader('chem','unweighted');
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
    tripletcountDDM(ii,:) = tripCount_chem(GenRandomMatrix3(datareader('chem','unweighted')));
%    tripletcountDDM(ii,:) = tripCount_chem(CreateRandomMatrix6(double(datareader('chem','unweighted'))));
end
%Finding mean, standard deviationa nd Z values
mean_Nrand_ciona = mean(tripletcountDDM);
std_Nrand_ciona = std(tripletcountDDM);
Zi_ciona = (tripletcount-mean_Nrand_ciona) ./ std_Nrand_ciona;
Znorm_ciona = Zi_ciona ./ (sumsqr(Zi_ciona))^0.5;
y2 = Znorm_ciona;
%Plotting triad significance profiles of Ciona intestinalis and C.elegans.
%Further, will group them into same superfamily if we get similar Triad
%significance profiles.
load('TSP_elegans')
plot(x,y1)
hold on
plot(x,y2);


