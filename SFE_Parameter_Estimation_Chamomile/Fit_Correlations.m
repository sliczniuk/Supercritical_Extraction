startup;
delete(gcp('nocreate'));

AA = xlsread('Regression.xlsx');
DI = [0.701472275, 1.331443055, 2.239307889, 2.711813187, 1.32629228, 1.485504345, 1.73827467, 2.59502961, 0.48656241, 1.363499511, 0.72227, 0.756214019];
GG = [4.274825704, 2.189390368, 2.552240039, 1.365163176, 2.830760407, 2.573487374, 1.642279591, 1.906200052, 4.287215235, 2.723682117, 3.82240, 3.35589348];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
FF = [6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 3.33, 3.33, 3.33, 3.33];

%%
x1 = RE(1:8) ; y1 = GG(1:8);
x2 = RE(9:12); y2 = GG(9:12);

[mdl1,stuff1] = fit(x1',y1','poly1');
[mdl2,stuff2] = fit(x2',y2','poly1');

p1 = coeffvalues(mdl1);
p2 = coeffvalues(mdl2);

xtest  = linspace(0.1,0.5);
ytest1 = polyval(p1,xtest);
ytest2 = polyval(p2,xtest);

hold on
plt1 = plot(x1,y1,'ko', 'LineWidth',2);
plt2 = plot(x2,y2,'kd', 'LineWidth',2);
plot(xtest,ytest1, 'LineWidth',2);
plot(xtest,ytest2, 'LineWidth',2);
%yline(0, 'k', 'LineWidth',2);
hold off

text(0.3, 4.0, sprintf('~$D_i = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p1, stuff1.rsquare]))
text(0.15, 6.0, sprintf('$~D_i = %.3f \\cdot Re + %.3f$ \n $R^2=%.3f$',[p2, stuff2.rsquare]))

legend([plt1(1),plt2(1)],'$3.33\times 10^{-5}[kg/s]$','$6.67\times 10^{-5}[kg/s]$','Location','best');
legend box off

%ylabel('$D_i^R \times 10^{-13} [m^2/s]$')
ylabel('$\Upsilon [-]$')
xlabel('Re [-]')
set(gca,'FontSize',12)

exportgraphics(figure(1), ['Correlation_Gamma_Re_.png'], "Resolution",300);

close all