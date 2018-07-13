for j1 =1:12
	for j2=(j1+1):12
		str1 = sprintf('benchmarks/test%d.mat',j1);
		str2 = sprintf('benchmarks/test%d.mat',j2);
		if exist(str1) && exist(str2) 
			load(str1,'problem_list');
			p1 = problem_list;
			load(str2,'problem_list');
			p2 = problem_list;
			if strcmp(p1{1},p2{1});

				%%%%%%%%%%%%%%%%%%%%%%%
				% 1
				%%%%%%%%%%%%%%%%%%%%%%%
				load(str1);

				tttt = size(Q,2);
				QQ=Q(:,1:(tttt-1));
				Tsum = T+Te;

				columnLabels = problem_list;
				rowLabels = solver_list;

				vars{1} = QQ;	nvars{1} = 'QQ';	
				vars{2} = T;	nvars{2} = 'T';	
				vars{3} = Te;	nvars{3} = 'Te';	
				vars{4} = NMI;	nvars{4} = 'NMI';	
				vars{5} = PUR;	nvars{5} = 'PUR';	
				vars{6} = NMI;	nvars{6} = 'NMI';	
				vars{7} = Tsum;	nvars{7} = 'Tsum';	

				%%%%%%%%%%%%%%%%%%%%%%%
				% 2
				%%%%%%%%%%%%%%%%%%%%%%%
				load(str2);

				tttt = size(Q,2);
				QQ=Q;
				Tsum = T+Te;

				columnLabels = problem_list;
				rowLabels = [rowLabels solver_list 'reference'];

				vars{1} = [vars{1} QQ];		nvars{1} = 'QQ';	ff{1} = '%-3.2f';
				vars{2} = [vars{2} T];		nvars{2} = 'T';		ff{2} = '%-3.2f';
				vars{3} = [vars{3} Te];		nvars{3} = 'Te';	ff{3} = '%-3.2f';
				vars{4} = [vars{4} NMI];	nvars{4} = 'NMI';	ff{4} = '%-3.2f';
				vars{5} = [vars{5} PUR];	nvars{5} = 'PUR';	ff{5} = '%-3.2f';
				vars{6} = [vars{6} NMI];	nvars{6} = 'NMI';	ff{6} = '%-3.2f';
				vars{7} = [vars{7} Tsum];	nvars{7} = 'Tsum';	ff{7} = '%-3.2f';


				for i=1:length(vars)
					matrix2latex(vars{i}',sprintf('Report/matlab_tables/%s%d.tex',nvars{i},j1),'columnLabels',columnLabels,'rowLabels',rowLabels,'format',ff{i});
				end
			end
		end
	end
end
