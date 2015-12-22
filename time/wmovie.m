function M = wmovie()

S_MAX = 10;
N = 1000;
Ne = 800;
Ni = 200;
M  = 100;
D  = 50;

files = dir(sprintf('%d/synapses-all-*.dat',D)); % csv of rows: n_src, n_dest, d, w
files = files(1:10);

delays = 1:D';
weight_bins = 1:S_MAX;%+0.5);


for fi=1:length(files)
    delays_distn = zeros(length(weight_bins), length(delays));
    
    synapses = load(sprintf('%d/%s', D, files(fi).name));
    exc_synapses = synapses(ismember(synapses(:,1),0:(Ne-1)),:);
    n_exc_synapses = size(exc_synapses,1);
    
    
    for di=1:length(delays)
        d_synapses = exc_synapses(:,3)==delays(di);

        delays_distn(:,di) = histc(exc_synapses(d_synapses,4), weight_bins)/n_exc_synapses;
    end;
    

    f = figure;
    surf(delays, weight_bins, delays_distn);
    set(gca,'ylim',[0 S_MAX]);
    xlabel('delay');
    ylabel('timestep in training');
    zlabel('proportion of synapses');

    if fi==1, M = getframe(gcf);
    else, M(fi) = getframe(gcf);
    end
    
    close(f);
end;

%subplot(1,2,1);
%bar(delays,h);
%set(gca, 'xlim', [delays(1)-0.5, delays(end)+0.5]);
%set(gca, 'ylim', [0 max(h(:))+0.01]);

%subplot(1,2,2);
%bar(delays,h2);
%set(gca, 'xlim', [delays(1)-0.5, delays(end)+0.5]);
%set(gca, 'ylim', [0 max(h(:))+0.01]);

if nargout==0,
    movie(M,1,2/length(M)); % make the movie 2s long
end;
