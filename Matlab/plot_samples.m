function plot_samples(samples,mu,vars, bins, name)

n = size(samples,1);

[nu,x] = hist(samples,bins);
h = x(2) - x(1);

bar(x,nu/(n*h),1);

xs = min(x):0.01:max(x);

actual_gaussian = normpdf(xs,mu,vars);

hold on;
plot(xs,actual_gaussian,'r','linewidth',1);
title(name);
xlabel('Energy','fontsize',10);
ylabel('Density','fontsize',10);
hold off;

sample_mean = mean(samples)
sample_vars = var(samples)

end
