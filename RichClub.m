%Rich club analysis. 
%Copyright 2018. Cheshta Bhatia & Lav R. Varshney

%load the adjacency matrix A
function A = RichClub(varargin)
%adjacency matrix
if (nargin == 0)
    %load the chemical network
    A = datareader('chem','unweighted');
elseif (nargin == 1)
    A = varargin{1};
else
    error('TRIPCOUNT_CHEM: incorrect number of inputs');
end
%Find in degree rich club coefficient(phi) and then characterise neurons
%that belong to rich club. Phi  = E/((N)(N-1)) where N is the number of
%neurons with indeg value>k and E is the number of incoming edges in
%those. Phi was calculated for Ciona chemical synapse network, random network and further Phi(normalised).

% 1. Phi(Ciona network)
indeg_ciona = sum(A);
k = 0:1:77;
phi_ciona = zeros(1,78);
for i = 0:1:77
    p_ciona = indeg_ciona;
    p_ciona(p_ciona <= i) = 0 ;
    N_ciona = size(p_ciona(p_ciona ~= 0 ),2);
    subg_indices = find(~p_ciona) ;
    A0 = A ;
    A0(subg_indices,:) = [] ;
    A0(:,subg_indices) = [] ;
    E_ciona = sum(sum(A0));
    phi_ciona(i+1) = (E_ciona)./((N_ciona)*(N_ciona-1)) ;
end
    
% 2. Phi(random network)
%Generate 1000 random networks that preserve directed degrees of A first.Then find phi(Rich club coefficient for the random matrix)

for j = 1:1000

% Get size of graph (number of vertices)
N = size(A,1);

% Compute in-degree and out-degree distributions for A
indegdist = sum(A);
outdegdist = sum(A');

% Create pool of available in-stubs
inpool = [];
for k = 1:N
    inpool = [inpool k*ones(1,indegdist(k))];
end

% Create pool of available out-stubs
outpool = [];
for k = 1:N
    outpool = [outpool k*ones(1,outdegdist(k))];
end

% Main Loop - keep trying to generate a random matrix until RM is generated
% sucessfully - generate RM by pairing out-stub on u with in-stub on v with
% the condition that u is not equal to v and edge (u,v) does not already
% exist
fail = true;
while fail  % Keep looping until RM is generated successfully
    % Initialize pools for current trial
    ip = inpool;
    op = outpool;
    
    % Initialize RM
    RM = false(N);
    
    % Begin current attempt to create RM by pairing in-/out-stubs
    fail = false;
    while ~isempty(op)
        % Generate a random in-stub/out-stub pair
        k1 = ceil(rand*length(op));
        k2 = ceil(rand*length(ip));
        u = op(k1);
        v = ip(k2);
        
        % If current pair violates conditions make a new pair
        pf = 0;   % Used to detect pairing failures
        while (u == v) || RM(u,v) % If either condition fails,
            % check for pairing failure (defined as 300 consecutive pairs
            % not meeting conditions)
            pf = pf + 1;
            if mod(pf,300) == 0
                fail = true;
                break;
            end
            
            % choose a new pair
            k1 = ceil(rand*length(op));
            k2 = ceil(rand*length(ip));
            u = op(k1);
            v = ip(k2);
        end        
        if fail
            break
        end
        
        % Now that a successful pair has been found,
        % add the edge to RM
        RM(u,v) = 1;
        
        % and remove the chosen stubs from their respective pools
        op(k1) = [];
        ip(k2) = [];
    end
end
    % Get size of graph (number of vertices)
indeg_rand = sum(RM);
%phi_rand = zeros(j,78);
for i = 0:1:77
    p_rand = indeg_rand
    p_rand(p_rand <= i) = 0 ;
    N_rand = size(p_rand(p_rand ~= 0 ),2);
    subg_indices_rand  = find(~p_rand) ;
    AR = RM ;
    AR(subg_indices_rand,:) = [] ;
    AR(:,subg_indices_rand) = [] ;
    E_rand = sum(sum(AR));
    phi_rand(j,i+1) = (E_rand)./((N_rand)*(N_rand-1)) ;
  end
end

%phi_rand gives the rich club coefficient for the random matrix.

indeg_rand = sum(RM);
mean_phi_rand = mean(phi_rand);
phi_ciona = phi_ciona(:,1:50);
mean_phi_rand = mean_phi_rand(:,1:50);
phi_norm = phi_ciona ./mean_phi_rand;
%To characterise neurons present in rich club
s = std(mean_phi_rand)
%Very stringent criterion used to characterize neurons in rich club. Neurons
%in which phi(norm) >=  1+3s
richclub = sum(phi_norm >=  1 + 3*s)
%richclub gives you the number of elements in rich club
indegree_richclub = k(find(phi_norm >=  1 + 3*s)) 
%indegree_richclub gives you the indegree(or k value) of neurons that %belong to rich club. 
%plot graph of rich club coefficient v/s in degree values
 
k = 0:1:49;
x = k
y1 = phi_ciona
plot(x,y1)
hold on
y2 = mean_phi_rand
plot(x,y2)
y3 = phi_norm
plot(x,y3)
legend('Phi(Ciona intestinalis)','phi(Random network)','phi(norm)')




