using SpecialFunctions
tol = 1.25e-16
x = -4.5;

function nm(f, fp, x, tol)
   
    ctr, max_steps = 0, 100
     
    while (abs(f(x)) > tol) && ctr < max_steps
        x = x - f(x) / fp(x)
        ctr = ctr + 1
    end

    ctr >= max_steps ? error("Method did not converge") : return (x, ctr)
    
end

function f(x)
	airyai(x)
end

function fp(x)
	airyaiprime(x)
end

val = airyai(nm(f,fp,x,tol)[1])
