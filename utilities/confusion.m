function M=confusion(g1,g2)
	u1=g2u(g1);
	u2=g2u(g2);
	M=u1'*u2;
end
