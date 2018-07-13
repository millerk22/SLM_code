% Plots each frame of the evolution of the graph clustering. 
u = g2u(g);
G = graph(logical(W));

% Treat the binary and n-ary case separately
if size(u,2) > 2

	[~,nodecolors] = max(u,[],2);
	%nodecolors = zeros(size(u(:,1)));
	%for a=1:length(u(1,:))
	%	nodecolors = nodecolors + a*u(:,a);
	%end

	% Hard coded color map to make all the colors very different from each other. Limited to 11 classes right now. There is probably a more elegant solution.
	map = [1 0 0; % 1
	0 1 0;
	0 0 1;
	1 1 0;
	1 0 1; % 5
	0 1 1;
	1 1 1;
	.5 0 0;
	0 .5 0;
	0 0 .5; % 10
	.5 .5 0
	.5 0 .5;
	0 .5 .5;
	.5 .5 .5;
	1 .5 0; %15
	1 0 .5;
	1 1 .5;
	1 .5 1;
	0 1 .5;
	0 .5 1; %20
	.5 0 1;
	.5 0 .5;
	.5 1 0;
	.5 1 .5;
	.5 1 1; %25
	.25 0 0;
	0 .25 0;
	0 0 .25;
	.25 .25 0;
	.25 0 .25; %30
	0 .25 .25;
	.25 .25 .25;
	1 .25 0;
	1 0 .25;
	1 .25 .25; %35
	1 1 .25;
	1 .25 0;
	.25 1 0;
	0 1 .25;
	.25 1 .25; %40
	.25 0 1;
	0 .25 1;
	.25 .25 1;
	1 1 .25;
	1 .25 1; %45
	.25 1 1;
	1 .5 .25;
	1 .25 .5;
	.5 1 .25;
	.25 1 .5 %50
	];
	map = map(1:size(u,2),:); % Restrict to the number of classes actually needed
	colormap(map);
	caxislen = length(u(1,:));

else

	nodecolors=u(:,1);
	colormap default;
	caxislen = 1;

end

G.Nodes.nodecolors = nodecolors;
gplot = plot(G,'layout','auto');
gplot.NodeCData = G.Nodes.nodecolors;
colorbar;
caxis manual;
caxis([0 caxislen]);
drawnow;
