function [A,tag,ind] = get_giant_component(A,tag)

	tagg = get_components(A);

	for i=1:max(tagg)
		siz(i) = nnz(tagg==i);
	end
	siz_max = max(siz);
	i_max = find(siz==siz_max,1);
	ind = find(tagg==i_max);
	A = A(ind,ind);
	tag = tag(ind);
	%[~,~,tag] = unique(tag);
end
