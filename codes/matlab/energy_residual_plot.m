[X,Y] = meshgrid(1:1:100,1:1:100);
F = exp(-log(2)/1000*X.*Y);
surf(X,Y,F)
xlabel('$$n_{iter}$$','interpreter','latex', 'FontSize', 18)
ylabel('$$n_{sym}$$','interpreter','latex', 'FontSize', 18)
zlabel('Residual energy ($$E_{residual}$$)','interpreter','latex','FontSize', 18)
title('$$E_{residual}= e^{-\frac{\log_e(2)}{10^3}n_{sym}n_{iter}}$$',...
    'FontSize', 18,'interpreter','latex')