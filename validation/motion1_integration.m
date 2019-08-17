%% perform computation

T = 1
h = T/4.


xp = [ 0.5; 0.1 ]

pos = [ xp ]
times = [ 0 ]

printf "\n\n*** first step ***\n\n"

x0 = pos(:,end)
t0 = times(end)
[v,dv] = motion(x0, t0)
k0 = v

x1 = x0 + 0.5*h * k0
t1 = t0 + 0.5*h
[v,dv] = motion(x1, t0)
k1 = (v + (t1 - t0) * dv)

x2 = x0 + 0.5*h * k1
t2 = t0 + 0.5*h
[v,dv] = motion(x2, t0)
k2 = (v + (t2 - t0) * dv)

x3 = x0 + h*k2
t3 = t0 + h
[v,dv] = motion(x3, t0)
k3 = (v + (t3 - t0) * dv)

xp = xp + h*(k0 + 2*k1 + 2*k2 + k3)/6.

pos(:,end+1) = xp;
times(end+1) = t3;

printf "\n\n*** second step ***\n\n"

x0 = pos(:,end)
t0 = times(end)
[v,dv] = motion(x0, t0)
k0 = v

x1 = x0 + 0.5*h * k0
t1 = t0 + 0.5*h
[v,dv] = motion(x1, t0)
k1 = (v + (t1 - t0) * dv)

x2 = x0 + 0.5*h * k1
t2 = t0 + 0.5*h
[v,dv] = motion(x2, t0)
k2 = (v + (t2 - t0) * dv)

x3 = x0 + h*k2
t3 = t0 + h
[v,dv] = motion(x3, t0)
k3 = (v + (t3 - t0) * dv)

xp = xp + h*(k0 + 2*k1 + 2*k2 + k3)/6.

pos(:,end+1) = xp;
times(end+1) = t3;

printf "\n\n*** third step ***\n\n"

x0 = pos(:,end)
t0 = times(end)
[v,dv] = motion(x0, t0)
k0 = v

x1 = x0 + 0.5*h * k0
t1 = t0 + 0.5*h
[v,dv] = motion(x1, t0)
k1 = (v + (t1 - t0) * dv)

x2 = x0 + 0.5*h * k1
t2 = t0 + 0.5*h
[v,dv] = motion(x2, t0)
k2 = (v + (t2 - t0) * dv)

x3 = x0 + h*k2
t3 = t0 + h
[v,dv] = motion(x3, t0)
k3 = (v + (t3 - t0) * dv)

xp = xp + h*(k0 + 2*k1 + 2*k2 + k3)/6.

pos(:,end+1) = xp;
times(end+1) = t3;

printf "\n\n*** forth step ***\n\n"

x0 = pos(:,end)
t0 = times(end)
[v,dv] = motion(x0, t0)
k0 = v

x1 = x0 + 0.5*h * k0
t1 = t0 + 0.5*h
[v,dv] = motion(x1, t0)
k1 = (v + (t1 - t0) * dv)

x2 = x0 + 0.5*h * k1
t2 = t0 + 0.5*h
[v,dv] = motion(x2, t0)
k2 = (v + (t2 - t0) * dv)

x3 = x0 + h*k2
t3 = t0 + h
[v,dv] = motion(x3, t0)
k3 = (v + (t3 - t0) * dv)

xp = xp + h*(k0 + 2*k1 + 2*k2 + k3)/6.

pos(:,end+1) = xp;
times(end+1) = t3;

%% output

times
pos

truepos = [ truemotion( pos(:,1), 0*h )',
	    truemotion( pos(:,1), 1*h )',
	    truemotion( pos(:,1), 2*h )',
	    truemotion( pos(:,1), 3*h )',
	    truemotion( pos(:,1), 4*h )' ]'

dpos = pos - truepos

for i=1:5
    err(i) = sqrt( dpos(:,i)' * dpos(:,i) );
end

err
