function f_t = ricker1(t,tp,ts)
aa = pi*(t-ts)/tp;
bb=pi*ts/tp;
f_t = (aa.^2 - bb^2).* exp(-aa.^2);