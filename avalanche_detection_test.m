%%% avalanche detection
%%% 27.11.2019.

%% function [av_duration, av_size] = avalanche_detection(m,theta)

% Input
%load('m_save.mat');
%x = smooth(m_save); % timeseries to be analized
t = 0:0.1:5.5*2*3.14;
x = sin(t);

% Parameters
theta = 0.0; % threshold
to = 0; % single avalanche start
summ = 0; % avalanche size variable
count = 0; % avalanche count
av_duration = [];
av_size = [];

% Loop
for it = 1:length(x)
    
    if x(it)>theta
        if to==0
            count = count + 1; % avalanche count
            to = it;
        end
        summ = summ + (x(it)-theta)*1;
    end
    
    % if currently inside an avalanche and it ends:
    % evaluated only if avalanche was detected but hasn?t yet finished (t~=0)
    if to~=0 && x(it)<theta
        t0(count) = to;
        tfin(count) = it;
        av_duration(count) = it-to;
        av_size(count) = summ;
        to = 0; % reset to zero
        summ = 0;
    end
    
end

% if last avalanche was not finished, reject/delete it
if to~=0
        count = count - 1;
end
 
%% Plot
figure(1)
plot(x)
hold on;
plot(theta*ones(size(x)),'-')
for it=1:count
    xline(t0(it),'b--',{num2str(it)});
    xline(tfin(it),'r--',{num2str(it)});
end

% Histograms
figure(3)
histogram(av_duration,100)
title('Avalanche duration histogram')
xlabel('avalanche duration')
figure(4)
histogram(av_size,100)
title('Avalanche size histogram')
xlabel('avalanche size')
figure(5)
scatter(av_duration,av_size)
title('Avalanche size vs. duration')
xlabel('avalanche duration')
ylabel('avalanche size')