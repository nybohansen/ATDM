
close all;

figure();
cont = samples(samples(:,1)==1,2);
histfit(cont,1000);
[muhat,sigmahat] = normfit(cont)


figure();
disc = samples(:,1);
histDisc = hist(disc,2)/length(disc)
bar([0,1], histDisc);



