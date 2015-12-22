shared_settings;


% All interhemispheric
D_INTERS = {20 37 50 83 100};
p.nSimulations = 10;

for i=1:length(D_INTERS)
	p.D_INTER = D_INTERS{i};
	
	fprintf('D_INTER = %d:\n', p.D_INTER);
	run_simulations(p, 'refresh');
end;

% one more for split
p.M_INTER = 0;
p.D_INTER = p.D_INTRA;
fprintf('SPLIT:\n', p.D_INTER);
run_simulations(p, 'refresh');