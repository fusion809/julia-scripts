using Pkg;
Pkg.add("FLoops");
using FLoops;
N=100000

@time @floop begin
	s = 0;
	for x in 1:N
		s += x;
	end
end

k = 0;
@time begin
	k = 0;
	for y in 1:N
		if (y == 1)
			k = y;
		else
			k += y;
		end
	end
end
