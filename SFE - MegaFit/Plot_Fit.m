figure(ii)
subplot(2,2,1)
hold on
plot(Time  , xx_0(Nx,:)  ,'k'   )
plot(Time  , xx_out(Nx,:),'r'   )
plot(SAMPLE, data_org    ,'b--o')
%plot(Time, 1e3 * (sum(xx_out(0*nstages+1:1*nstages,:) .* V_fluid) + sum(xx_out(1*nstages+1:2*nstages,:) .* V_solid)) + xx_out(Nx,:));
hold off

xlabel('Time [min]')
ylabel('Yield [g]')

%colorbar

%legend('Inital curve','Optimized curve','Data points','Location','northwest'); legend('boxoff')

ax=gca;
ax.FontSize = 20;

subplot(2,2,2)
hold on
plot(SAMPLE(2:end), diff(xx_0(Nx,N_Sample)  ),'k'   )
plot(SAMPLE(2:end), diff(xx_out(Nx,N_Sample)),'r'   )
plot(SAMPLE(2:end), diff(data_org           ),'b--o')
%plot(Time, 1e3 * (sum(xx_out(0*nstages+1:1*nstages,:) .* V_fluid) + sum(xx_out(1*nstages+1:2*nstages,:) .* V_solid)) + xx_out(Nx,:));
hold off

xlabel('Time [min]')
ylabel('d Yield / dt [g/s]')

%colorbar

legend('Inital curve','Optimized curve','Data points','Location','northeast'); legend('boxoff')

ax=gca;
ax.FontSize = 20;

subplot(2,2,3)
imagesc(Time,L_nstages,xx_out(1:nstages,:)); colormap turbo;
a = colorbar;
a.Label.String = 'C_f';
a.Label.FontSize = 16;
%caxis([0, 15]);

xlabel('Time [min]')
ylabel('Length [m]')

ax=gca;
ax.FontSize = 16;

subplot(2,2,4)
imagesc(Time,L_nstages,xx_out(1*nstages+1:2*nstages,:)); colormap turbo;
a = colorbar;
a.Label.String = 'C_s';
a.Label.FontSize = 16;
%caxis([0, 70]);

xlabel('Time [min]')
ylabel('Length [m]')

ax=gca;
ax.FontSize = 16;

%name = sprintf('$k_m = %g, D_i =%g e-14, D_x = %g e-6, C_{sat}=%g, Inital~ratio = %g, \\sigma$ = %g',round(KOUT,3));
%name = sprintf('%g',round(KOUT,3));
%sgtitle(name, 'Interpreter', 'latex')

set(gcf,'PaperOrientation','landscape'); print(figure(ii),['Fitting_',DATA_set{ii},'.pdf'],'-dpdf','-bestfit'); close all