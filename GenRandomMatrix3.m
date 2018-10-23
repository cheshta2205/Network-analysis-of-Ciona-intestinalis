%To generate random matrix such that the random matrix preserves degree
%distributions of matrix A.

function RM = GenRandomMatrix3 (A)
% GENRANDOMMATRIX3(A) produces a directed random graph preserving the degree
% distribution of A

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
