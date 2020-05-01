close all;
clear all;



M = csvread("upscaledData.csv");
M2 = csvread("resolvedData.csv");
h=figure(1);
hold on;
plot(M(2:end,8),M(2:end,1),"--b",'linewidth',2)
xlabel ("x");
ylabel ("log(c)");
axis ([0 16.5 -50 0]);

plot(M2(2:end,8),M2(2:end,1),"-r",'linewidth',2)
xlabel ("x");
ylabel ("log(c)");
axis ([0 16.5 -100 0]);

W = 4; H = 3;

box on;
legend("upscaled model","direct simulation","location","northeast")
legend boxoff;
print(h,'-depsc','-color','res_up_compare.eps')