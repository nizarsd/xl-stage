function y = onesidenorm(mu,sig)	
	y=0;	
	while (y < mu)
		y = normrnd(mu,sig);
	end
end


