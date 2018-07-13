function [avePurity]=purityMeas(group,tag)
	
	%cleaning
	if size(group,1)<size(group,2)
		group=group';
	end
	if size(tag,1)<size(tag,2)
		tag=tag';
	end
	[~,~,group]=unique(group);
	[~,~,tag]=unique(tag);

	Purity = 0;
	for i = 1:max(group)
		index = find(group == i);
		plurality = max(hist(tag(index), 1:max(tag)));
		Purity = Purity + plurality;
	end
	avePurity=Purity/length(group);
end
