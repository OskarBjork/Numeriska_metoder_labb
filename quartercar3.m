function [tv,yv]=Euler(f,tspan,y0,n)
%Solve y'=f(t,y) with inital condition y(a)=y0. 
a=tspan(1);
b=tspan(2); 
h=(b-a)/n;

tv=a+h*(0:n); 
yv=zeros(size(tv));
yv(1)=y0;

for ii=1:n
    yv(ii+1)=yv(ii)+h*f(tv(ii),yv(ii));
end;