function summary = hp_summary(files, sort_field)
% Returns a cell array of struct arrays
%
% Each struct array represents a simulation run
%   with particular parameter settings
%
    if (~exist('sort_field','var')), sort_field = 'D_INTER'; end;
    
    if (~exist('files','var'))
        if (exist('summary.mat', 'var'))
            load('summary');
        else
            files = dir('*.mat');
            files = setdiff({files.name}, 'summary.mat');
            summary   = hp_summary(files);
            save('summary', 'summary');
        end;
        
        return;
    end;
    

    summary  = {};
    sort_val = [];
    
    for i=1:length(files)
        % Process into usable form
        s = hp_process_to_summary(files{i});
        if (isempty(s)), continue; end;
        % Check whether to append to an existing cell,
        %   or to put in a new cell
        for j=1:length(summary)
            if (   s.p.M       == summary{j}(end).p.M ...
                && s.p.M_INTER == summary{j}(end).p.M_INTER ...
                && s.p.Ne      == summary{j}(end).p.Ne ...
                && s.p.D_INTRA == summary{j}(end).p.D_INTRA ...
                && s.p.D_INTER == summary{j}(end).p.D_INTER )
                
                summary{j}(end+1) = s;
                clear s;
                break;
            end;
        end;

        if (exist('s','var'))
            summary{end+1} = s;
            sort_val(end+1) = s.p.(sort_field);
            clear s;
        end;
    end;
    
    % Reorder based on M
    [a,idx] = sort(sort_val);
    summary = summary(idx);
    
    
    function out = hp_process_to_summary(f);
    %
    

        % Load the data
        fprintf('Processing data for %s ...', f);
        load(f);
        
        keyboard
        if (~exist('poly','var')), out = []; return; end;
		
		if (0.05 <= 1-(poly.summary(end,1)+1)/p.Ne)
			warning('Possible incomplete run: last poly group found with mother neuron #%d/%d', poly.summary(end,1)+1,p.Ne);
		end;
		
		
        % fix up c++ coding issue
        if (p.M_INTER == 0), p.D_INTER = inf; end;        
        
        if (isempty(poly.spikes) || isempty(poly.summary))
            warning('No polychronous groups detected!');
            poly.spikes = zeros(0,4); 
            poly.summary = zeros(0,4);
            
        elseif (max(poly.spikes(:,1)) ~= max(poly.summary(:,1)))
            error('Poly spikes vs summary files incompatible; %d vs %d.  Adjusting to lower number...', max(poly.spikes(:,1)), max(poly.summary(:,1)));
        end;


        %%

        % Save settings and make some local aliases
        out.p = p;
        out.p.N = out.p.Ne+out.p.Ni;
        
        M  = out.p.M; Ne = out.p.Ne; Ni = out.p.Ni;
       
        % Calculate summary
        out.fr.inh = mean(spnet.summary(:,4));
        out.fr.exc = mean(spnet.summary(:,2));

        % Calculate spnet.conns
        post = reshape(spnet.conns(:,1), [M size(spnet.conns,1)/M])';
        s    = reshape(spnet.conns(:,2), [M size(spnet.conns,1)/M])';
        lh_ih_exc= (Ne/2<= post(1:Ne/2,:) & post(1:Ne/2,:) < Ne);
        lh_ih_inh= ((Ne+Ni/2)<=post(1:Ne/2,:));
        rh_ih_exc= (Ne/2>  post(Ne/2+1:Ne,:));
        rh_ih_inh= (Ne <= post(Ne/2+1:Ne,:) & post(Ne/2+1:Ne,:) < (Ne+Ni/2));
        ih   = [lh_ih_exc | lh_ih_inh; ...
                rh_ih_exc | rh_ih_inh]; % logical indices of interhemispheric connections

        out.conns.sn    = numel(s);
        out.conns.sm_tn = find(s>=9.5); % # strong excitatory synapses
        out.conns.sm_en = find(s>=9.5 & post < Ne);  % # exc-exc strong synapses
        out.conns.sm_in = find(s>=9.5 & post >= Ne); % # exc-inh strong synapses

        out.conns.inter.sm_tn = find(ih & s>=9.5);
        out.conns.inter.sm_en = find(ih & s>=9.5 & post < Ne);  % # exc-exc strong synapses
        out.conns.inter.sm_in = find(ih & s>=9.5 & post >= Ne); % # exc-inh strong synapses

        out.conns.intra.sm_tn = find(~ih & s>=9.5);
        out.conns.intra.sm_en = find(~ih & s>=9.5 & post < Ne);  % # exc-exc strong synapses
        out.conns.intra.sm_in = find(~ih & s>=9.5 & post >= Ne); % # exc-inh strong synapses

        % Calculate dist'n of # interhemispheric connections
        out.conns.inter.exc = find([lh_ih_exc;rh_ih_exc]);
        out.conns.inter.inh = find([lh_ih_inh;rh_ih_inh]);

        %%

        % Calculate polychronous groups
    %    if (exist([file_stem '.poly.summarymary.dat'], 'file'))

            out.poly.lh_ih_exc= (poly.spikes(:,1)<Ne/2 & Ne/2<=poly.spikes(:,3) & poly.spikes(:,3) < Ne);
            out.poly.lh_ih_inh= (poly.spikes(:,1)<Ne/2 & (Ne+Ni/2)<=poly.spikes(:,3));
            out.poly.rh_ih_exc= (Ne/2<=poly.spikes(:,1) & poly.spikes(:,1)<Ne & poly.spikes(:,3)<Ne/2);
            out.poly.rh_ih_inh= (Ne/2<=poly.spikes(:,1) & poly.spikes(:,1)<Ne & Ne<=poly.spikes(:,3) & poly.spikes(:,3)<(Ne+Ni/2));

            out.poly.ih_exc   = unique(poly.spikes(out.poly.lh_ih_exc | out.poly.rh_ih_exc,2));
            out.poly.ih_inh   = unique(poly.spikes(out.poly.lh_ih_inh | out.poly.rh_ih_inh,2));
            out.poly.ih       = unique([out.poly.ih_exc; out.poly.ih_inh]);
            out.poly.tot      = size(poly.summary,1);

            out.poly.max_len = poly.summary(:,3);
            out.poly.max_time= poly.summary(:,4);

             out.poly.ih_spikes = poly.spikes(ismember(poly.spikes(:,2),out.poly.ih), :);
   %    end;



        %%
        fprintf('done.\n');