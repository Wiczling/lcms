function plot_uncertainity_chromatogram(samples_tr,expidx,metidx)

map = colormap('lines');
for i=1:1:length(metidx) 
tr = squeeze(samples_tr(:,metidx(i),expidx));
hold on
[f,xi] = ksdensity(tr,0:0.1:(max(tr)+5));
plot(xi,f,'Color',map(i,:));
end