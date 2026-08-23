%%% Selver Pepic, PhD Student, IST Austria
%%% 06.11.2019. v1
%%% 29.11.2019. v2/3
%%% Rotation #1, Tkacik group
%%% Statistics of neural avalanches near the critical point

%%% implements a negative (temporal) feedback Ising model, which leads to interesting dynamics
%%% especailly around the critical point, where random activations (avalanches) and orderly oscillations 
%%% can be detected by examining the time-series of total spin (magnetization).
%%% Based on the model developed by: De Martino et al. 2019 "Feedback-induced self-oscillations in large interacting systems subjected to phase transitions. Journal of Physics A"
%%% - see project report and publication for more details

%% Simulation and physics parameters
niter = 5000000;
nsweeps = 10;
nx = 32;
ny = 32;
N = nx*ny;
boundary_periodic = 0;
connected_global = 0; % fully connected if 1, near neigh if 0

B = 0; %magnetic field, initial value
tau = 0.1;
J = 1; % spin-spin coupling
beta = 0.45; % 1/kT, inverse temperature
theta = 0.0; % threshold

%% variables
di = [1,0,-1,0]; % neighborhod directions
dj = [0,1,0,-1]; % neighborhod directions
co = 0.9;
c = rand(nx,ny);
c = -1*(c<=co) + 1*(c>co); % random matrix, average = 0
%c_save = zeros(nx,ny,niter,'single');
%c_save(:,:,1) = c;
m_save = zeros(niter,1);
m_save(1) = sum(sum(c))/N;
B_save = zeros(niter,1);
B_save(1) = B;
E_save = zeros(niter,1);
if connected_global == 1
   E_save(1) = -B*N*m_save(1) - 0.5*J*4*N*m_save(1)^2;
else
   E_save(1) = energy(c,B,J);
end

%%%%%%%%%%%%%%%
%% Main loop %%
%%%%%%%%%%%%%%%
av_duration = [];
av_size = [];
inter_duration = [];

for sweep=1:nsweeps    
    
    tic
    for it = 2:niter    
        % flip one spin at i/j_flip
        c_new = c;
        i_flip = randi([1,nx]);
        j_flip = randi([1,ny]);
        c_new(i_flip,j_flip) = -1*c(i_flip,j_flip);

        % calculate dE
        % fully connected
        if connected_global == 1
            dE = 2*c(i_flip,j_flip) * (B + J*4*(m_save(it-1)-2*c(i_flip,j_flip)/N));
        else
            c_neigh_sum = 0;
            for d = 1:4
                if boundary_periodic == 1
                    c_neigh_sum = c_neigh_sum + ...
                        c(mod(i_flip+di(d)-1,nx)+1,mod(j_flip+dj(d)-1,ny)+1 );
                    %(0,nx-1) mod nx = (0,nx-1), add 1 => (1,nx)
                elseif (i_flip+di(d)>=1 && i_flip+di(d)<=nx && ...
                        j_flip+dj(d)>=1 && j_flip+dj(d)<=ny)            
                    c_neigh_sum = c_neigh_sum + c(i_flip+di(d),j_flip+dj(d));
                end
            end
            % nearest neighbour interaction only    
            dE = 2 * c(i_flip,j_flip) * (B + J*c_neigh_sum);
        end

        % accept or reject proposed flip
        if dE<0
            c = c_new;
            E_save(it) = E_save(it-1) + dE;
        elseif rand()<exp(-beta*dE)
            c = c_new;
            E_save(it) = E_save(it-1) + dE;
        else
            E_save(it) = E_save(it-1);
        end

        %save new configuration
        %c_save(:,:,it) = c;
        m_save(it) = sum(sum(c))/N;

        % update magnetic field
        B = B - m_save(it)/tau/N; % update external field
        B_save(it) = B;
    end
    toc
    [av_duration_temp, av_size_temp,inter_duration_temp] = avalanche_detection(m_save,theta);
    av_duration = [av_duration av_duration_temp];
    av_size = [av_size av_size_temp];
    inter_duration = [inter_duration inter_duration_temp];    
end

%% Plot
figure(1)
plot(B_save)%/max(B_save))
title({["B (norm.) and m vs. ""time"" "],...
    [num2str(nx) ' x ', num2str(ny) ' spins'],[...
    'beta = ', num2str(beta)]})
hold on;
plot(m_save)
xlabel('# iteration')
ylabel('B, m')
legend('B','m')

figure(2)
plot(m_save,B_save)
title({["B vs m"],...
    [num2str(nx) ' x ', num2str(ny) ' spins'],[...
    'beta = ', num2str(beta)]})
xlabel('m')
ylabel('B')

%% Histograms 1
figure(3)
hold on;
% average of av_size for fixed av_duration
av_size_mean = [];
av_size_num = [];
for i=1:length(av_duration)
    av_size_mean(i) = mean(av_size(av_duration == av_duration(i)));
    av_size_num(i) = length(av_size(av_duration == av_duration(i)));
end
scatter(log10(av_duration),log10(av_size_mean),'.')

% cut data, fit
cut_dur_low = 1.5;
cut_dur_up = 2.5;
av_duration_cut = av_duration( log10(av_duration)<cut_dur_up & log10(av_duration)>cut_dur_low );
av_size_mean_cut = av_size_mean( log10(av_duration)<cut_dur_up & log10(av_duration)>cut_dur_low );
p = polyfit(log10(av_duration_cut),log10(av_size_mean_cut),1)
%cut_size_low = -2;
%cut_size_up = 4;
cut_size_low = p(1)*cut_dur_low + p(2);
cut_size_up = p(1)*cut_dur_up + p(2);
hold on;
plot(log10(av_duration_cut),polyval(p,log10(av_duration_cut)),'LineWidth',2);
title({['Avalanche size vs. duration (log-log)'],...
    [num2str(nx) 'x', num2str(ny) ' spins', ', \beta = ', num2str(beta), ', \theta = ', num2str(theta)],...
     ['slope = ', num2str(p(1))]})
xlabel('log(avalanche duration)')
ylabel('log(avalanche size)')
legend('simulation data',['y = ' num2str(p(1)),'*x + ', num2str(p(2))])

%% Histograms 2
nbins = 30;

figure(4)
%hold on;
edges_duration_min = min(av_duration);
edges_duration_max = max(av_duration);
edges_duration = exp(linspace(log(edges_duration_min),log(edges_duration_max),nbins));
h = histogram(av_duration,edges_duration,'Normalization','pdf');
fig = gca;
fig.XScale = 'log';
fig.YScale = 'log';

% extract histogram data
av_dur_prob = h.Values;
av_dur = [];
for i=1:nbins-1
    av_dur(i) = 0.5*(edges_duration(i)+edges_duration(i+1));
    %av_dur(i) = edges_duration(i);
end
av_dur = av_dur(av_dur_prob~=0);
av_dur_prob = av_dur_prob(av_dur_prob~=0);
av_dur_log = log10(av_dur);
av_dur_prob_log = log10(av_dur_prob);

% plot data
figure(5)
plot(av_dur_log,av_dur_prob_log,'--o')
p = polyfit(av_dur_log(av_dur_log<cut_dur_up & av_dur_log>cut_dur_low),...
    av_dur_prob_log(av_dur_log<cut_dur_up & av_dur_log>cut_dur_low),1)
hold on;
plot(av_dur_log(av_dur_log<cut_dur_up & av_dur_log>cut_dur_low),...
    polyval(p,av_dur_log(av_dur_log<cut_dur_up & av_dur_log>cut_dur_low)));
title({['Avalanche duration histogram (log-log)'],...
     [num2str(nx) 'x', num2str(ny) ' spins', ', \beta = ', num2str(beta), ', \theta = ', num2str(theta)],...
     ['slope = ', num2str(p(1))]})
xlabel('log(avalanche duration)')
ylabel('log(probability)')
legend('simulation data',['y = ' num2str(p(1)),'*x + ', num2str(p(2))])

figure(6)
%hold on;
edges_size_min = min(av_size);
edges_size_max = max(av_size);
edges_size = exp(linspace(log(edges_size_min),log(edges_size_max),nbins));
h = histogram(av_size,edges_size,'Normalization','pdf');
fig = gca;
fig.XScale = 'log';
fig.YScale = 'log';

% get histogram data
av_siz_prob = h.Values;
av_siz = [];
for i=1:nbins-1
    av_siz(i) = 0.5*(edges_size(i)+edges_size(i+1));
end
av_siz = av_siz(av_siz_prob~=0);
av_siz_prob = av_siz_prob(av_siz_prob~=0);
av_siz_log = log10(av_siz);
av_siz_prob_log = log10(av_siz_prob);

% plot data
figure(7)
plot(av_siz_log,av_siz_prob_log,'--o')
p = polyfit(av_siz_log(av_siz_log<cut_size_up & av_siz_log>cut_size_low),...
    av_siz_prob_log(av_siz_log<cut_size_up & av_siz_log>cut_size_low),1)
hold on;
plot(av_siz_log(av_siz_log<cut_size_up & av_siz_log>cut_size_low),...
    polyval(p,av_siz_log(av_siz_log<cut_size_up & av_siz_log>cut_size_low)));
title({['Avalanche size histogram (log-log)'],...
     [num2str(nx) 'x', num2str(ny) ' spins', ', \beta = ', num2str(beta), ', \theta = ', num2str(theta)],...
     ['slope = ', num2str(p(1))]})
xlabel('log(avalanche size)')
ylabel('log(probability)')
legend('simulation data',['y = ' num2str(p(1)),'*x + ', num2str(p(2))])

%% interduration (ISI)
nbins = 30;
cut_inter_dur_low = cut_dur_low;
cut_inter_dur_up = cut_dur_up;

figure(6)
hold on;
inter_duration = inter_duration(inter_duration~=0);
edges_duration_min = min(inter_duration);
edges_duration_max = max(inter_duration);
edges_duration = exp(linspace(log(edges_duration_min),log(edges_duration_max),nbins));
h = histogram(inter_duration,edges_duration,'Normalization','pdf');
%fig = gca;
%fig.XScale = 'log';
%fig.YScale = 'log';
% extract histogram data
inter_dur_prob = h.Values;
inter_dur = [];
for i=1:nbins-1
    inter_dur(i) = 0.5*(edges_duration(i)+edges_duration(i+1));
end
inter_dur = inter_dur(inter_dur_prob~=0);
inter_dur_prob = inter_dur_prob(inter_dur_prob~=0);
inter_dur_log = log(inter_dur);
inter_dur_prob_log = log(inter_dur_prob);
% plot data
plot(inter_dur_log,inter_dur_prob_log,'--o')
p = polyfit(inter_dur_log(inter_dur_log<cut_inter_dur_up & inter_dur_log>cut_inter_dur_low),...
    inter_dur_prob_log(inter_dur_log<cut_inter_dur_up & inter_dur_log>cut_inter_dur_low),1)
hold on;
plot(inter_dur_log(inter_dur_log<cut_inter_dur_up & inter_dur_log>cut_inter_dur_low),...
    polyval(p,inter_dur_log(inter_dur_log<cut_inter_dur_up & inter_dur_log>cut_inter_dur_low)));
title({['ISI histogram (log-log)'],...
     [num2str(nx) 'x', num2str(ny) ' spins', ', \beta = ', num2str(beta), ', \theta = ', num2str(theta)],...
     ['slope = ', num2str(p(1))]})
xlabel('log(ISI)')
ylabel('log(probability)')
legend('simulation data',['y = ' num2str(p(1)),'*x + ', num2str(p(2))])