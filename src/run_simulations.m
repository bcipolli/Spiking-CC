function run_simulations(p,m)
%
% p : parameter structure
% m : string "mode"


    p = fixup_params(p);
    if (~exist('m','var')), m='run'; end;

    if (~exist(fullfile(pwd,p.dirname),'dir')), mkdir(p.dirname); end;
    cd(p.dirname);

    set_rseed = ~(isfield(p, 'rseed'));
    
    %% Run simulations
    for i=1:p.nSimulations
    	try
    		% Set some variables
	        if (set_rseed), p.rseed        = (i+1); end; % > 1 is best
			p.prefix       = sprintf('%s.%d', p.dirname, p.rseed);
			p.log_file     = sprintf('out/%s.spnet.txt', p.prefix);
			p.summary_file = sprintf('mat/%s.summary.mat', p.prefix);

			p.sp_full_fn     = sprintf('dat/%s.spnet.dat',    p.prefix);
			p.sp_spikes_fn   = sprintf('dat/%s.spikes.dat',   p.prefix);
			p.sp_summary_fn  = sprintf('dat/%s.summary.dat',  p.prefix);
			p.sp_conns_fn    = sprintf('dat/%s.conns.dat',    p.prefix);

			p.poly_links_fn  = sprintf('dat/%s.poly.links.dat',   p.prefix);
			p.poly_results_fn= sprintf('dat/%s.poly.results.dat', p.prefix);
			p.poly_spikes_fn = sprintf('dat/%s.poly.spikes.dat',  p.prefix);
			p.poly_summary_fn= sprintf('dat/%s.poly.summary.dat', p.prefix);
		
			p.spnet     = sprintf('spnet.%s', p.prefix);
			p.polychron = sprintf('polychron.%s', p.prefix);
		
			% Rewrite files
			p.start_idx = 0;
			p.end_idx = p.Ne;
		
       		% Check what we need to do
	        if (exist(fullfile(pwd, p.summary_file), 'file') && ~p.force),
    	    	fprintf('\tSkipping rseed=%d: existing mat file.\n', p.rseed);
        		continue; 
       		end;
       	
%        % Do the simulation
			switch (m)
				case {'refresh'}
					make_output(p, 0);
				
				otherwise
		        	make_binaries(p);
        			run_programs(p);
		        	make_output(p);
		        
			        % Email
 	       			my_unix('mail -s "Ran on $HOSTNAME: %s/%s" %s < %s', p.spnet, p.polychron, p.eAddr, p.log_file);
			end;
			


	   	catch 
	   		le  = lasterror;
   			str = sprintf('%s:%d: %s', le.stack(1).file, le.stack(1).line, le.message);

   			fprintf('%s\n', str);
	        my_unix('echo "%s" | mail -s "ERROR on $HOSTNAME: spnet/polychron (%d)" %s', str, p.rseed, p.eAddr);
   			continue;
	   	end;
    end;
    
    % Signal completion
    cd('..');
    if (~strcmp(m, 'refresh')),
    	my_unix('echo "" | mail -s "Completed on $HOSTNAME: run_simulations" %s', p.eAddr);
    end;

    %%%%%%%%%%%%%%%%%%%%%%%
    function make_binaries(p)

		if (~exist(fullfile(pwd,'src'),'dir')), mkdir('src'); end;
        if (~exist(fullfile(pwd,'dat'),'dir')), mkdir('dat'); end;

		% Resume
		[p.start_idx, p.poly_cnt] = find_start_idx(p);

        % Write const.cpp
        fh = fopen('src/const.cpp','w');
        fprintf(fh, 'const	int		M_INTRA=%d;\n', p.M_INTRA);
        fprintf(fh, 'const	int		M_INTER=%d;\n', p.M_INTER);
        fprintf(fh, 'const	int		M      = M_INTRA+M_INTER;\n'); 
        fprintf(fh, 'const  int		D_INTRA=%d;\n', p.D_INTRA);
        fprintf(fh, 'const	int		D_INTER=%d;\n', p.D_INTER);
        fprintf(fh, 'const	int		D=D_INTER;\n');		
        fprintf(fh, 'const	int		Ne =%d;\n', p.Ne);		
        fprintf(fh, 'const	int		Ni =%d;\n', p.Ni);
        fprintf(fh, 'const	int		N = Ne+Ni;\n');	
		fprintf(fh, '\n');
        fprintf(fh, 'const	float	NHOURS = %f;\n', p.nHours);
        fprintf(fh, 'const	int		RSEED  = %d;\n', p.rseed);
        fprintf(fh, 'const	int		POLY_START_NEURON  = %d;\n', p.start_idx);
        fprintf(fh, 'const	int		POLY_END_NEURON    = %d;\n', p.end_idx);
        fprintf(fh, 'const	int		POLY_START_NUMBER  = %d;\n', p.poly_cnt+1);
		fprintf(fh, '\n');
        fprintf(fh, 'const	double	C_max=10;\n');		
        fprintf(fh, 'const	int		W=3;	// initial width of polychronous groups\n');
        fprintf(fh, 'const	int		min_group_path = 7;		// minimal length of a group\n');
        fprintf(fh, 'const	int		min_group_time = 40;	// minimal duration of a group (ms)\n');
        fprintf(fh, 'const	int		N_firings_max	=200000; //so that this simulation can run efficiently\n');
		fprintf(fh, '\n');
		fprintf(fh, 'const  char*   sp_full_fn     = "%s"; //\n', p.sp_full_fn);
		fprintf(fh, 'const  char*   sp_spikes_fn   = "%s"; //\n', p.sp_spikes_fn);
		fprintf(fh, 'const  char*   sp_summary_fn  = "%s"; //\n', p.sp_summary_fn);
		fprintf(fh, 'const  char*   sp_conns_fn    = "%s"; //\n', p.sp_conns_fn);
		fprintf(fh, '\n');
		fprintf(fh, 'const  char*   poly_links_fn  = "%s"; //\n', p.poly_links_fn);
		fprintf(fh, 'const  char*   poly_results_fn= "%s"; //\n', p.poly_results_fn);
		fprintf(fh, 'const  char*   poly_spikes_fn = "%s"; //\n', p.poly_spikes_fn);
		fprintf(fh, 'const  char*   poly_summary_fn= "%s"; //\n', p.poly_summary_fn);
        fclose(fh);

        % Link to shared source
        src_files = {'defs.h', 'polychron.cpp', 'spnet.cpp'};
        for i=1:length(src_files)
        	full_path  = fullfile(pwd,'src',src_files{i});
        	if (exist(full_path, 'file')), continue; end;
        	
        	my_unix('rm -f %s', full_path); % silently delete any broken linked file
        	my_unix('ln -s %s %s', fullfile(p.srcPath, src_files{i}), full_path);
        end;

        % Compile binaries
        if (exist(p.spnet, 'file')), my_unix('rm %s', p.spnet); end;
        my_unix('g++ -Wno-deprecated -o %s src/spnet.cpp', p.spnet);

        if (exist(p.polychron, 'file')), my_unix('rm %s', p.polychron); end;
        my_unix('g++ -Wno-deprecated -o %s src/polychron.cpp', p.polychron);

		% Save off the source
		my_unix('mv src/const.cpp src/const.%s.cpp', p.prefix);

        
    %%%%%%%%%%%%%%%%%%%%%%%    
    function [start_idx,poly_cnt] = find_start_idx(p)
		if (validate_spnet_output(p)~=0 || validate_polychron_output(p) ~= 2)
		    start_idx = 0;
		    poly_cnt = 0;
		    return;
		end;

		psum = load(p.poly_summary_fn);  
		if (isempty(psum)), 
			start_idx = 0;
			poly_cnt = 0;
			return;
		end;
		
		fprintf('\tSetting up resume  ... ');
										psum = psum(:,1);
		pspi = load(p.poly_spikes_fn);   pspi = pspi(:,1);
		plin = load(p.poly_links_fn);    plin = plin(:,1);
			
		start_idx = psum(end);
		poly_cnt  = length(find(psum<psum(end)));
		fprintf(' at index %d (%d polygroups found) ... ', start_idx, poly_cnt);

		% Rewrite summary file 
		fhi = fopen(p.poly_summary_fn, 'r');
		fho = fopen('tmp.txt', 'w');
		for i=1:length(find(psum<psum(end)))
			fprintf(fho, fgets(fhi));
		end;
		fclose(fhi); fclose(fho);
		my_unix('mv %s %s', 'tmp.txt', p.poly_summary_fn);
			
			
		% Rewrite spikes file
		fhi = fopen(p.poly_spikes_fn, 'r');
		fho = fopen('tmp.txt', 'w');
		for i=1:length(find(pspi<psum(end)))
			fprintf(fho, fgets(fhi));
		end;
		fclose(fhi); fclose(fho);
		my_unix('mv %s %s', 'tmp.txt', p.poly_spikes_fn);
			
		
		% Rewrite links file
		fhi = fopen(p.poly_links_fn, 'r');
		fho = fopen('tmp.txt', 'w');
		for i=1:length(find(plin<psum(end)))
			fprintf(fho, fgets(fhi));
		end;
		fclose(fhi); fclose(fho);
		my_unix('mv %s %s', 'tmp.txt', p.poly_links_fn);
			
		fprintf('done.\n');
		
		
    %%%%%%%%%%%%%%%%%%%%%%%    
    function p = fixup_params(p)
        if (~isfield(p, 'M')), p.M = p.M_INTRA + p.M_INTER; end;
        if (~isfield(p, 'D')), p.D = max(p.D_INTRA, p.D_INTER); end;
        if (~isfield(p, 'dirname')), 
            if (p.M_INTER==0), p.dirname = sprintf('split.%d', p.M);
            else,               p.dirname = sprintf('%d.%d', p.M, p.D_INTER); end;
        end;
		if (~isfield(p, 'force')), p.force=false; end;
    

    %%%%%%%%%%%%%%%%%%%%%%%    
    function run_programs(p)

        if (~exist(fullfile(pwd,'out'),'dir')), mkdir('out'); end;
        
        % Run spnet
        switch (validate_spnet_output(p))
        	case 0, fprintf('\t(skipping completed spnet result.)\n');
        	
        	case 1, 
        			[ss,ws] = my_unix('./%s > %s', p.spnet, p.log_file);
        			
        	case 2, fprintf('\t(overwriting incomplete run; stopped at t=%d)\n', s(end,1)); 
        	        [ss,ws] = my_unix('./%s > %s', p.spnet, p.log_file);
        	        
        	case 3, fprintf('\t(redoing run with missing conns)\n');
         	        [ss,ws] = my_unix('./%s > %s', p.spnet, p.log_file);
	   	end;
    	
		switch (validate_polychron_output(p))
			case 0, fprintf('\t(skipping completed polychron result.)\n'); 
			
			% Resume is now implemented
			case {1,2}, [sp,wp] = my_unix('./%s >> %s', p.polychron, p.log_file);

%			case 2, warning('Possible incomplete run: last poly group found with mother neuron #%d/%d', poly.summary(end,1)+1,p.Ne);
%					fprintf('Delete file to re-run.\n');
		end;

		
    %%%%%%%%%%%%%%%%%%%%%%%
    function b = validate_spnet_output(p)
		if (~exist(p.sp_summary_fn, 'file')), b = 1; return; end;
		
        s = load(p.sp_summary_fn);
        
        if (s(end,1)+1~=60*60*p.nHours),   b = 2; return; end;
        if (~exist(p.sp_conns_fn,'file')), 
        	if (exist(strrep(p.sp_conns_fn, 'conns', 'connections'), 'file'))
        		my_unix('mv %s %s', strrep(p.sp_conns_fn, 'conns', 'connections'), p.sp_conns_fn);
        	else
        		b = 3; return; %
			end;
		end;
		
		%            	
		b = 0;

	    
    %%%%%%%%%%%%%%%%%%%%%%%
    function b = validate_polychron_output(p)

        if (~exist(p.poly_summary_fn, 'file')), b = 1; return; end;

       	ps = load(p.poly_summary_fn);
       	if (isempty(ps)),                  b = 2; return; end;
		if (0.05 <= 1-(ps(end,1)+1)/p.Ne), b = 2; return; end;
		
		b = 0;
		
	    
    %%%%%%%%%%%%%%%%%%%%%%%
    function make_output(p,m)
    %
    %	m = make message
    
    	if (~exist('m','var')), m=1; end;
    	
    	vs = validate_spnet_output(p);
		if (~ismember(vs, 0))
			if (m), warning('Bad spnet output (%d) for rseed = %d?', vs, p.rseed); end;
			return;
		end;
			
		vp = validate_polychron_output(p);
		if (~ismember(vp, 0))
			if (m), warning('Bad polychron output (%d) for rseed = %d?', vp, p.rseed); end;
			return;
		end;

		fprintf('\tmaking summary for rseed = %d\n', p.rseed);

		
		% Make summary mat file
        if (~exist(fullfile(pwd,'mat'),'dir')), mkdir('mat'); end;
        
        spnet.conns   = load(p.sp_conns_fn);
        spnet.summary = load(p.sp_summary_fn);
        poly.links    = load(p.poly_links_fn);
        poly.spikes   = load(p.poly_spikes_fn);
        poly.summary  = load(p.poly_summary_fn);
        
        if (exist(p.sp_delays_fn, 'file'))
        	spnet.delays = load(p.sp_delays_fn);
        else
			delay = zeros(p.Ne+p.Ni, p.M);
			for i=1:p.Ne
				for j=1:p.M_INTRA
					delay(i,j) = ceil(j / (p.M_INTRA/p.D_INTRA));
				end;
			
				for j=1:p.M_INTER
					delay(i,j) = p.D_INTER;
				end;
			end;
		
			for i=(p.Ne+1):(p.Ne+p.Ni)
				delay(i,1:p.M) = 1;
			end;
			
			spnet.delays = delay; 
			clear('delay');
		end;
		
		
        save(p.summary_file, 'p', 'spnet', 'poly');
        clear('spnet','poly');
        
        
		% Move out data files
%        my_unix('mv %s dat', p.sp_conns_fn);
%        my_unix('mv %s dat', p.sp_full_fn);
%        my_unix('mv %s dat', p.sp_spikes_fn);
%        my_unix('mv %s dat', p.sp_summary_fn);
%        
%        my_unix('mv %s dat', p.poly_links_fn);
%        my_unix('mv %s dat', p.poly_results_fn);
%        my_unix('mv %s dat', p.poly_spikes_fn);
%        my_unix('mv %s dat', p.poly_summary_fn);



    %%%%%%%%%%%%%%%%%%%%%%%
	function [s,t] = my_unix(varargin)
		cmd = sprintf(varargin{:});
		
        c = clock; tic;
        fprintf('%02d:%02d:%02d : ', c(4), c(5), round(c(6)));
        fprintf('Running "%s" ...', cmd);
        [s,t] = unix(cmd);
        
	    if (s == 0), fprintf('OK (%6.2fs).\n', toc);
	    else, fprintf('\nError:\t%s\n\n', t);
	    end;


    %%%%%%%%%%%%%%%%%%%%%%%
%    function load_balance();
%        persistent servers = {'desi','dino','garcia','janis','marley'};
%        [s,cur]  = my_unix('echo $HOSTNAME');
%        [s,user] = my_unix('echo $USER');
%        
%        for i=1:length(servers)
%            if (isempty(findstr(cur,servers{i})))
%                my_unix('ssh 
%            [s,str] = my_unix('ps -ef | grep %s | grep MATLAB', user);
%            
%            if (~isempty(findstr(cur,servers{i})))
%                
%              
