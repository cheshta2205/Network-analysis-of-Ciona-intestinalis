%2018
%MATLAB Code for the diverse club analysis of Ciona intestinalis
%Copyright Cheshta Bhatia and Lav R. Varshney

%Load the weighted matrix of Chemical synapse network of Ciona intestinalis.
function A = DiverseClub(varargin)
%adjacency matrix
if (nargin == 0)
    %load the chemical network
    A = datareader('chem','weighted')
elseif (nargin == 1)
    A = varargin{1};
else
    error('TRIPCOUNT_CHEM: incorrect number of inputs');
end
%Find participation coefficient. It denotes how evenly distributed are
%nodes edges across various communities(Ci). If all edges of a node are to a single community, 
%its participation coefficient is zero. However, Participation coefficient is maximum if a node
%has an equal sum of edge weights to each community of a network.

n=length(A);                        %nodes present in network
K1=sum(A,2) ;                        % out degree of a node
Gc=(A~=0)*diag(Ci);                 %connection with community 
K2=zeros(n,1);                      

for i=1:max(Ci)
   K2=K2+(sum(A.*(Gc==i),2).^2);
end
%PC denotes the participation coefficient. PC is denoted by the formula
%stated below.
PC=ones(n,1)-K2./(K1.^2);
PC(~K1)=0;                           %P=0 if for nodes with no (out)neighbors
%Neurons present in the diverse club should have a participation
%coefficient with value greater than 80 percentile.
Percentile = prctile(PC,80);   %Find 80th percentile value.
diverseclub_members = find(PC>Percentile); 
save('diverseclub_members');

