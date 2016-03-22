function y = limitednorm(mu,sig,mn,mx)	
	y=0;	
	while (y <= mn || y >= mx)
		y = normrnd(mu,sig);
	end
end


