function wdist(mode)

if ~exist('mode','var'), mode='old'; end;

S_MAX = 10;
N = 1000;
Ne = 800;
Ni = 200;
M  = 100;
D  = 50;

switch mode
    case 'old', files = dir(sprintf('%d/synapses*.dat',D)); % csv of rows: n_src, n_dest, d, w
    case 'all', files = dir(sprintf('%d/synapses-all*.dat',D)); % csv of rows: n_src, n_dest, d, w
    case 'big', files = dir(sprintf('%d/synapses-big*.dat',D)); % csv of rows: n_src, n_dest, d, w
end;

delays = 1:D';
delays_distn = zeros(1+length(files), length(delays));

for fi=1:length(files)
    synapses = load(sprintf('%d/%s', D, files(fi).name));
    
    exc_synapses = synapses(ismember(synapses(:,1),0:(Ne-1)),:);

    switch mode
        case {'old','all'}, big_exc_idx   = (exc_synapses(:,4)>=S_MAX*0.95);
        case 'big',         big_exc_idx   = true(size(exc_synapses,1),1);
    end;

    big_exc_delays = exc_synapses(big_exc_idx,3);
    all_exc_delays = exc_synapses(:,3);

    % pre-training distribution
    if fi==1
        switch mode
            case {'old','all'}, delays_distn(1,:) = histc(all_exc_delays,delays) / length(all_exc_delays); 
            case 'big',         delays_distn(1,:) = ones(size(delays))/length(delays);
        end;
    end;
    delays_distn(1+fi,:)  = histc(big_exc_delays,delays) / length(all_exc_delays);
    
end;



figure;
surf(delays, 1:size(delays_distn,1),delays_distn);
xlabel('delay');
ylabel('timestep in training');
zlabel('proportion of synapses');
view(33.5, -40);

%figure;
%subplot(1,2,1);
%bar(delays,h);
%set(gca, 'xlim', [delays(1)-0.5, delays(end)+0.5]);
%set(gca, 'ylim', [0 max(h(:))+0.01]);

%subplot(1,2,2);
%bar(delays,h2);
%set(gca, 'xlim', [delays(1)-0.5, delays(end)+0.5]);
%set(gca, 'ylim', [0 max(h(:))+0.01]);


