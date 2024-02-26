function p=mpsnr(x,y)
x=double(x);
y=double(y);
p = mean((x(:)-y(:)).^2);
m1 = max( abs(x(:)) );
    m2 = max( abs(y(:)) );
    vmax = max(m1,m2);
p=10*log10(vmax^2/p);