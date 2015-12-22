function hp_figures(figs, summary)
%
%
%

    close all
    
    if (~exist('figs','var'))
        figs = 1:10;
    end;

    if (~exist('summary','var'))
        load('summary.mat', 'summary');
    end;
    
    % Get out labels
    lbls = cell(length(summary),1);
    for i=1:length(summary)
%        lbls{i} = sprintf('M=%d, M_{IH}=%d, D=%d', summary{i}.args.M, summary{i}.args.M_INTER, summary{i}.args.D);
        if (summary{i}.args.M_INTER == 0)
            lbls{i} = sprintf('SPLIT, N=%d', summary{i}.args.N);
        else
            lbls{i} = sprintf('D=%d, N=%d', summary{i}.args.D, summary{i}.args.N);
        end;
        
        if (summary{i}.args.pct_poly ~= 1.0)
            lbls{i} = [lbls{i} ' **'];
        end;
    end;

    fig_firing_rate  (intersect(1,figs),     summary, lbls);
    fig_synapses     (intersect([2 3],figs), summary, lbls);
    fig_ih_targets   (intersect(4,figs),     summary, lbls);
    fig_poly_maxtime (intersect(5,figs),     summary, lbls);
    fig_poly_maxlen  (intersect(6,figs),     summary, lbls);
    fig_poly_qty     (intersect([7 8],figs), summary, lbls);
    %fig_poly_distn   (intersect(9,figs),     summary, lbls);

    
    
%function fig_poly_len_distn(fn, summary, lbls)
% Show the 
    %figure(fn);
    %



%%%%%%%%%%%%%%%%%%%%
function fig_firing_rate(fn, summary, lbls)
    if (isempty(fn)), return; end;
    
    figure(fn);
    b = zeros(length(summary), 2);
    for i=1:length(summary)
        b(i,:) = [summary{i}.fr.exc summary{i}.fr.inh]; 
    end;
    bar(b');
    set(gca, 'xtick', 1:2, 'xticklabel', {'Excitatory', 'Inhibitory'});
    set(gca, 'FontSize', 18);
    xlabel('Neuron type'); ylabel('Hz');
    legend(lbls, 'Location', 'NorthWest');
    title('Mean Firing Rates', 'FontSize', 24);

    
    
%%%%%%%%%%%%%%%%%%%%
function fig_synapses(fns, summary, lbls)
    if (isempty(fns)), return; end;

    b = zeros(length(summary), 3);
    for i=1:length(summary)
                b(i,:) = [length(summary{i}.cxns.sm_tn)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                              length(summary{i}.cxns.sm_en)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                              length(summary{i}.cxns.sm_in)       / (summary{i}.args.Ne*summary{i}.args.M) ];
    end;
    bar(b');
    set(gca, 'xtick', 1:3, 'xticklabel', {'All', 'Excitatory', 'Inhibitory'});
    set(gca, 'FontSize', 18);
    xlabel('Neuron type'); ylabel('%');
    legend(lbls, 'Location', 'NorthEast');
    title('% Strong Synapses', 'FontSize', 24);


%%%%%%%%%%%%%%%%%%%%
function fig_synapses_old(fns, summary, lbls)
    if (isempty(fns)), return; end;

    figure(fns(1)); 
    %Plots for total, intra, and inter-hemsipheric excitatory connections
    sp = 1;
    for ih=1:4
        for p=1:3
            subplot(4,3,sp);

            b        = zeros(0, 3);
            loc_lbls = {};

            for i=1:length(summary)
                switch (ih)
                    case 1, if (~ismember(summary{i}.args.M_INTER, [0])),  continue; end;
                    case 2, if (~ismember(summary{i}.args.M_INTER, [2])),  continue; end;
                    case 3, if (~ismember(summary{i}.args.M_INTER, [20])), continue; end;
                    case 4, if (~ismember(summary{i}.args.M_INTER, [50])), continue; end; 
                end;

                switch (p)
                    case 1
                        b(end+1,:) = [length(summary{i}.cxns.sm_tn)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                                      length(summary{i}.cxns.intra.sm_tn) / (summary{i}.args.Ne*summary{i}.args.M_INTRA) ...
                                      length(summary{i}.cxns.inter.sm_tn) / (summary{i}.args.Ne*summary{i}.args.M_INTER)];
                    case 2
                        b(end+1,:) = [length(summary{i}.cxns.sm_en)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                                      length(summary{i}.cxns.intra.sm_en) / (summary{i}.args.Ne*summary{i}.args.M_INTRA) ...
                                      length(summary{i}.cxns.inter.sm_en) / (summary{i}.args.Ne*summary{i}.args.M_INTER)]; 
                    case 3
                        b(end+1,:) = [length(summary{i}.cxns.sm_in)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                                      length(summary{i}.cxns.intra.sm_in) / (summary{i}.args.Ne*summary{i}.args.M_INTRA) ...
                                      length(summary{i}.cxns.inter.sm_in) / (summary{i}.args.Ne*summary{i}.args.M_INTER)]; 
                end;
                loc_lbls{end+1} = lbls{i};
            end;

            if (~isempty(b))
                bar(b');
                set(gca, 'xtick', 1:3, 'xticklabel', {'Total', 'Intra', 'Inter'});
                 ylabel('% Synapses');
                legend(loc_lbls, 'Location', 'NorthEast');
                set(gca, 'ylim', [0 1]);
            end;

            switch (p)
                case 1, xlabel('% Strong Exc Synapses');
                case 2, xlabel('% Strong Exc->Exc Synapses');
                case 3, xlabel('% Strong Exc->Inh Synapses');
            end;        
            sp = sp + 1;
        end;
    end;
    
    
    figure(fns(2)); 
    %Plots for total, intra, and inter-hemsipheric excitatory connections
    sp = 1;
    for ih=1:3
        for p=1:3
            subplot(3,3,sp);

            b        = zeros(0, 3);
            loc_lbls = {};

            for i=1:length(summary)
                switch (ih)
                    %case 1, if (~ismember(summary{i}.args.M_INTER, [0])), continue; end;
                    case 1, if (~ismember(summary{i}.args.D, [20 21])),   continue; end;
                    case 2, if (~ismember(summary{i}.args.D, [50])),      continue; end;
                    case 3, if (~ismember(summary{i}.args.D, [100])),     continue; end; 
                end;

                switch (p)
                    case 1
                        b(end+1,:) = [length(summary{i}.cxns.sm_tn)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                                      length(summary{i}.cxns.intra.sm_tn) / (summary{i}.args.Ne*summary{i}.args.M_INTRA) ...
                                      length(summary{i}.cxns.inter.sm_tn) / (summary{i}.args.Ne*summary{i}.args.M_INTER)];
                    case 2
                        b(end+1,:) = [length(summary{i}.cxns.sm_en)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                                      length(summary{i}.cxns.intra.sm_en) / (summary{i}.args.Ne*summary{i}.args.M_INTRA) ...
                                      length(summary{i}.cxns.inter.sm_en) / (summary{i}.args.Ne*summary{i}.args.M_INTER)]; 
                    case 3
                        b(end+1,:) = [length(summary{i}.cxns.sm_in)       / (summary{i}.args.Ne*summary{i}.args.M) ...
                                      length(summary{i}.cxns.intra.sm_in) / (summary{i}.args.Ne*summary{i}.args.M_INTRA) ...
                                      length(summary{i}.cxns.inter.sm_in) / (summary{i}.args.Ne*summary{i}.args.M_INTER)]; 
                end;
                loc_lbls{end+1} = lbls{i};
            end;

            bar(b');
            set(gca, 'xtick', 1:3, 'xticklabel', {'Total', 'Intra', 'Inter'});
             ylabel('% Total Synapses');
            legend(loc_lbls, 'Location', 'NorthEast');
            set(gca, 'ylim', [0 1]);


            switch (p)
                case 1, xlabel('% Strong Exc Synapses');
                case 2, xlabel('% Strong Exc->Exc Synapses');
                case 3, xlabel('% Strong Exc->Inh Synapses');
            end;        
            sp = sp + 1;
        end;
    end;

    
%%%%%%%%%%%%%%%%%%%%
function fig_ih_targets(fn, summary, lbls)
    if (isempty(fn)), return; end;

    figure(fn);
    b = zeros(length(summary),2);
    for i=1:length(summary)
        b(i,:) = [length(summary{i}.cxns.inter.exc)   / (summary{i}.args.Ne*summary{i}.args.M_INTER) ...
                  length(summary{i}.cxns.inter.inh)   / (summary{i}.args.Ne*summary{i}.args.M_INTER)]; 
    end;
    bar(b');
    set(gca, 'xtick', 1:2, 'xticklabel', {'Excitatory', 'Inhibitory'});
    legend(lbls, 'Location', 'NorthEast');
    title('Interhemispheric Target Neuron Types');
    xlabel('Neuron Type'); ylabel('% Total Synapses');
    



%%%%%%%%%%%%%%%%%%%%
function fig_poly_maxtime(fn, summary, lbls)
    if (isempty(fn)), return; end;

    %% Polychron summary
    % subplot 1: average max time
    % subplot 2: average max time, interhemispheric excitatory 
    % subplot 3: average max time, interhemispheric inhibitory
    figure(fn);
    b = zeros(length(summary), 3);
    for i=1:length(summary)
        if (isfield(summary{i}, 'poly'))
            b(i,:) = [mean(summary{i}.poly.max_time) ...
                      mean(summary{i}.poly.max_time(summary{i}.poly.ih_exc)) ...
                      mean(summary{i}.poly.max_time(summary{i}.poly.ih_inh))]; 
        end;
    end;
    bar(b');
    set(gca, 'xtick', 1:3, 'xticklabel', {'All', 'Ih-Excitatory', 'Ih-Inhibitory'});
    set(gca, 'FontSize', 18);
    legend(lbls, 'Location', 'NorthWest');
    title('Mean maximum time', 'FontSize', 24);
    xlabel('Post-synaptic neuron type');
    ylabel('max time (ms)');



%%%%%%%%%%%%%%%%%%%%
function fig_poly_maxlen(fn, summary, lbls)
    if (isempty(fn)), return; end;
    % subplot 1: average max len
    % subplot 2: average max len, interhemispheric excitatory 
    % subplot 3: average max len, interhemispheric inhibitory
    %mean(out.poly.max_time(out.poly.ih_exc))
    %mean(out.poly.max_time(out.poly.ih_inh))
    figure(fn);
    b = zeros(length(summary), 3);
    for i=1:length(summary)
        if (isfield(summary{i}, 'poly'))
            b(i,:) = [mean(summary{i}.poly.max_len) ...
                      mean(summary{i}.poly.max_len(summary{i}.poly.ih_exc)) ...
                      mean(summary{i}.poly.max_len(summary{i}.poly.ih_inh))]; 
        end;
    end;
    bar(b');
    set(gca, 'xtick', 1:3, 'xticklabel', {'All', 'Ih-Excitatory', 'Ih-Inhibitory'});
    set(gca, 'FontSize', 18);
    legend(lbls, 'Location', 'NorthWest');
    title('Mean maximum length', 'FontSize', 24);
    xlabel('Post-synaptic neuron type');
    ylabel('max length (neurons)');


 


%%%%%%%%%%%%%%%%%%%%
function fig_poly_qty(fns, summary, lbls)
    if (isempty(fns)), return; end;
    % Polychron numbers
    % subplot 1: total polychron
    % subplot 2: total interhemispheric polychron
    % subplot 3: % interhemisperhic polychron
    figure(fns(1));
    b = zeros(length(summary),1);
    for i=1:length(summary)
        if (isfield(summary{i}, 'poly'))
            b(i,:) = [summary{i}.poly.tot] / summary{i}.args.pct_poly; 
        end;
    end;
    bar(b);
    set(gca, 'xticklabel', lbls);
    set(gca, 'FontSize', 18);
    title('Total # polychronous groups', 'FontSize', 24);
    
    %legend(lbls, 'Location', 'NorthEast');

    % Polychron numbers
    % subplot 1: total polychron
    % subplot 2: total interhemispheric polychron
    % subplot 3: % interhemisperhic polychron
    figure(fns(2));
    b = zeros(length(summary), 2);
    for i=1:length(summary)
        if (isfield(summary{i}, 'poly'))
            b(i,:) = [length(summary{i}.poly.ih_exc) / summary{i}.poly.tot ...
                      length(summary{i}.poly.ih_inh) / summary{i}.poly.tot]; 
        end;
    end;
    bar(b');% (find(b'));
    set(gca, 'xtick', 1:2, 'xticklabel', {'Excitatory', 'Inhibitory'});
    legend(lbls, 'Location', 'NorthEast');
    set(gca, 'FontSize', 18);
    title('% Interhemispheric polychronous groups', 'FontSize', 24);

