close all;


cont = sample(sample(:,1)==1,2);

figure();
plot_samples(cont, 0.47693, 0.45096, 200, 'Samples compared to the Gaussian of the model', 'r'); 

figure();
disc = sample(:,1);
histDisc = hist(disc,2)/length(disc)
bar([0,1], histDisc);
title('Distribution of indicator in the samples');
xlabel('Percentage','fontsize',10);
ylabel('Indicator','fontsize',10);


aa = energy_CO_mat(energy_CO_mat(:,2)==1,:);
bb = energy_NH_mat(energy_NH_mat(:,2)==1,:);

figure();
xs = plot_samples(aa(aa(:,3)>-5,3),-1.6746, 0.895861 , 200, 'Energy CO', 'r'); 
hold on;
plot(xs,normpdf(xs,-2.41637,0.213894),'g','linewidth',1);


figure();
xs = plot_samples(bb(bb(:,3)>-5,3),-0.97857, 0.47761 , 200, 'Energy NH', 'r'); 
hold on;
plot(xs,normpdf(xs,-2.1716, 0.54614 ),'g','linewidth',1);


