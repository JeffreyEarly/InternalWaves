j_star = 3;
H = (j_star+(1:1024)).^(-5/2);
H_norm = 1/sum(H);
H = H_norm*(j_star+(1:1024)).^(-5/2);

figure, plot(squeeze(sum(sum(abs(obj.w_plus).^2+ abs(obj.w_minus).^2,2),1))./(H'/H_norm))
