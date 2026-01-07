function F = extract_features(x,y,fs,W)

x0 = x - movmean(x,W);
y0 = y - movmean(y,W);

vx = diff(x0)*fs;
vy = diff(y0)*fs;

v  = sqrt(vx.^2 + vy.^2);
E  = vx.^2 + vy.^2;

win = hamming(W);
V   = fft(v(1:W).*win);
P   = abs(V).^2;
f   = (0:length(P)-1)*(fs/length(P));
[~,idx] = max(P(2:floor(end/2)));
fdom = f(idx+1);

F = [mean(v), mean(E), fdom, std(v)]
end
