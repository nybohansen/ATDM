close all;


cont = sample(sample(:,1)==1,2);

figure();
plot_samples(cont, 0.47693, 0.45096, 200, 'Samples compared to the Gaussian of the model'); 

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
plot_samples(aa(aa(:,3)>-10,3),-2.0957, 0.64378, 200, 'Energy CO'); 
figure();  
plot_samples(bb(bb(:,3)>-10,3),-2.0957, 0.64378, 200, 'Energy NH');


