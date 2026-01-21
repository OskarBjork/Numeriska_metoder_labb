m1 = 475
m2 = 53
k1 = 5400
k2 = 135000
c1 = 310
c2 = 1200
% 
%h0 = 0
%hp0 = 0
%h(1) = h0
%hp(1) = hp0
h = 0
hp = 0

z0 = 0
zp0 = 0
zpp0 = 0
z(1) = z0
zp(1) = zp0
zpp(1) = zpp0
x0 = 0
xp0 = 0
xpp0 = 0
x(1) = x0
xp(1) = xp0
xpp(1) = xpp0
v = [z; zp;x;xp]


size([[1,1];[1 1]])

t0 = 0
t(1) = t0
tn = 3
s = .5
n = (tn-t0)/s

H = [0 0 1 0;0 0 0 1;-k1/m1 k1/m1 -c1/m1 c1/m1;k1/m2 -(k1+k2)/m2 c1/m2 -(c1+c2)/m2]
g = [0;0;0;(k2*h+c2*hp)/m2]

size(H)
size(v(1))
size(g)
size(H*v)


dv = @(t, v) H*v + g

for i=1:n
    v(i+1) = v(i) + s*dv(t(i), v(i))
    t(i+1) = t0 + i.*s
end

% dv/dt = Hv + g(t)
