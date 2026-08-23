%%% Selver Pepic, PhD Student, IST Austria
%%% 06.11.2019.
%%% Rotation #1, Tkacik group
%%% Statistics of neural avalanches near the critical point

%%% Ising model with multiple averages (nsweeps), using Wollff algorithm (wolff2 function)
%%% to build and detect same-spin cluster and flip them together, which notably reduced the 
%%% critical slow-down, i.e. number of correlated samples near the critical point, and thus 
%%% leads to much faster evolution towarsd the final system (average) state.


%% Simulation and physics parameters
niter = 300;
nsweeps = 10;
nx = 64;
ny = 64;
N = nx*ny;
boundary_periodic = 0;
connected_global_local = 0; % fully connected if 1, near neigh if 0

B = 0; % external field;
J = 1; % spin-spin coupling
beta = 0.5; % 1/kT, inverse temperature
cluster_sizes = [];

%%%%%%%%%%%%%%%
%% Main loop %%
%%%%%%%%%%%%%%%
for sweep = 1:nsweeps

    disp({num2str(sweep),' out of ', num2str(nsweeps), ' iterations started.'});
    c = rand(nx,ny);
    co = 0.8;
    c = -1*(c<=co) + 1*(c>co); % random matrix, average = 0
    m_save = zeros(niter,1);
    m_save(1) = sum(sum(c))/N;

    tic
    for it = 2:niter
        %disp({num2str(it),' out of ', num2str(niter), ' iterations started.'});
        % build a cluster and flip it
        c = wolff2(c,B,J,beta);
        m_save(it) = sum(sum(c))/N;    
        %if it>20
        %    m_cut = m_save(round(0.7*it):it);    
        %    if std(abs(m_cut))/mean(abs(m_cut)) < 0.01
        %        break;
        %    end
        %end
    end
    toc
    
% calculate cluster histogram for this sweep
cluster_sizes = [cluster_sizes [cluster_detect(c) cluster_detect(-1*c)]];
%cluster_detect( ((sum(sum(c)))<0)*c -((sum(sum(c)))>0)*c )];
end

%% Plot
nbins = 100;
clsize_cut = 3;
figure(1)
%h = histogram(log(cluster_sizes),length(cluster_sizes),'Normalization','probability');
h = histogram(log(cluster_sizes),nbins,'Normalization','probability');
title('Cluster size histogram')
xlabel('cluster size')
ylabel('cluster frequency')

%%
hold on;
figure(2)
prob = h.Values;
clsize = h.BinEdges(2:end);
clsize = clsize(prob~=0);
prob = prob(prob~=0);
scatter(clsize,log(prob),'*')
hold on;
%p = polyfit(clsize,log(prob),1)
%exp_raw = p(1);
%plot(clsize,polyval(p,clsize));

% cleaner version
prob = prob(clsize<clsize_cut);
clsize = clsize(clsize<clsize_cut);
p = polyfit(clsize,log(prob),1)
exp_clean = p(1);
hold on;
plot(clsize,polyval(p,clsize));

title('Cluster histogram log-log')
xlabel('log(cluster size)')
ylabel('log(cluster frequency)')
legend(['Histogram data points'],...
       ['p = ', num2str(exp_clean), ', fit to all clusters size < ', num2str(clsize_cut)])
      %['p = ', num2str(exp_raw), ', fit to all cluster sizes'],...