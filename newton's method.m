function p = newton(f, df,i,tol)
% Solve f(p) = 0 using Newton's method. 
while 1
	p = i - (f(i)/df(i)); 
	if abs(i - p) < tol
		break 
	end
	i = p;
	end 
end
	
	yy3 = @(x) splineeval(t, a1 -1.2, b2, c2, d2, x); 
	yyd3 = @(x) diffsplineeval(t,a1-1.2,b2,c2,d2,x); 
	p1 = newton(yy3,yyd3,4.5,10^-6)
	p1 = newton(yy3,yyd3,2.1,10^-6)


function val = compositetrap(f,a,b,n) score=0;
ap = f(a);
bp = f(b);
h = (b-a)/n;
x = [a+h:h:b-h]; 
for i=1:n-1
	score =score+f(x(i));
end
score=h*(score+0.5*(ap+bp)); 
val=score;
end


newf = @(x) ((diffsplineeval(t,a,b1,c1,d1,x)).^2 + (diffsplineeval(t,a1,b2,c2,d2,x)).^2).^(1/2) 
L = compositetrap(newf, 2.317982170678300, 4.661644158641098, 10000)
L16 = compositetrap(newf, 2.317982170678300, 4.661644158641098, 16)
L32 = compositetrap(newf, 2.317982170678300, 4.661644158641098, 32)
L64 = compositetrap(newf, 2.317982170678300, 4.661644158641098, 64) 
L128 = compositetrap(newf, 2.317982170678300, 4.661644158641098, 128)



diff = 4.661644158641098-2.317982170678300
er1 = abs(L16-L); 
er2 = abs(L32-L); 
er3 = abs(L64-L); 
er4 = abs(L128-L); 
h1 = diff/16; 
h2 = diff/32; 
h3 = diff/64; 
h4 = diff/128;
erra = [er1,er2,er3,er4]; 
ha = [h1,h2,h3,h4];
p = polyfit(log(ha),log(erra),1) % slope
loglog(erra,ha);
grid on;