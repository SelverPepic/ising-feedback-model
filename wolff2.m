% INPUT spin configuration, grow a cluster at a random position and flip it
% OUTPUT new spin configuration
function c = wolff2(c,B,J,beta)

%global boundary_periodic;
%global connected_global_local; % fully connected if 1, near neigh if 0
nx = size(c,1);
ny = size(c,2);
cluster = zeros(size(c));
di = [1,0,-1,0]; % neighborhod directions
dj = [0,1,0,-1]; % neighborhod directions

% random starting point
i_flip = randi([1,nx]);
j_flip = randi([1,ny]);
cluster(i_flip,j_flip) = 1;
c0 = c(i_flip,j_flip); % reference spin
neigh{1} = [i_flip,j_flip]; % generate initial neighbourhood seed!
new_neigh = cell(0,0); % temporary variable
p_wolff = 1-exp(-2*beta*J);

%% build cluster until neigh is empty
while ~isempty(neigh) 
    
    for it = 1:length(neigh)
        % coordinates of the current neighbor
        i = neigh{it}(1);
        j = neigh{it}(2);
        
        for d = 1:4
            % if not already in cluster
            % with prob 1-e^-beta*dE form bond
            r = rand();
            if (i+di(d)<=nx && i+di(d)>=1 && j+dj(d)<=ny && j+dj(d)>=1) &&...
                    cluster(i+di(d),j+dj(d))==0 && c(i+di(d),j+dj(d))==c0 &&...
                    (r < p_wolff)
                % add to cluster and to neighbourhood
                cluster(i+di(d),j+dj(d)) = 1;
                new_neigh = [new_neigh [i+di(d),j+dj(d)]];
            end
        end
        
    end % end of neigh loop
    
    % list of all new neighbours (i.e. the unvisited ones)
    neigh = new_neigh;
    new_neigh = cell(0,0); % temporary variable    
end

%% cluster formed, now flip spins inside the cluster
c(cluster==1) = -1*c(cluster==1);