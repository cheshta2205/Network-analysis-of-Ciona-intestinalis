function A = datareader(network,weights)
%DATAREADER Reads Ciona intestinalis connectivity data.
%   A = DATAREADER(N,W) takes strings N and W that specify the network,
%   N in {'gap','chem'}, and whether to consider it as weighted, W in
%   {'weighted','unweighted'}, and returns the adjacency matrix.
%
%   [A,L] = DATAREADER(N,W) additionally returns the neuron labels L.
%
%   [A,L,C] = DATAREADER(N,W) additionally returns the neuron class labels C.

%   Copyright 2018. Cheshta Bhatia and Lav R. Varshney

load ConnOrdered
load NeuronType_Ordered

if isequal(network,'gap') && isequal(weights,'weighted')
    A = ConnOrdered.gap_total;
    
elseif isequal(network,'gap') && isequal(weights,'unweighted')
    %threshold
    A = ConnOrdered.gap_total > 0;
    

elseif isequal(network,'chem') && isequal(weights,'weighted')
    A = ConnOrdered.chem_total;

elseif isequal(network,'chem') && isequal(weights,'unweighted')
    %threshold
    A = ConnOrdered.chem_total > 0;
else
    error('DATAREADER: incorrect input parameters.')
end

%first output is the adjacency matrix
varargout(1) = {A};

%second output is the neuron labels
varargout(2) = {ConnOrdered.Neuron_ordered};

%third output is the neuron class labels 
varargout(3) = {NeuronType_Ordered};


