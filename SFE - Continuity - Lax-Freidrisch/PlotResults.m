function [] = PlotResults(xx_0, Time, nstagesbed, Parameters, uu, bed_mask, epsi, flag)

nstages = Parameters(1);

feedTemp  = uu(:,1);
feedPress = uu(:,2);
feedFlow  = uu(:,3);


Cs_NS  = xx_0(1*nstages+1:2*nstages,:);
T_NS   = xx_0(2*nstages+1:3*nstages,:);

if isequal(flag, 'conservative')
    Cf_NS  = xx_0(0*nstages+1:1*nstages,:) ./ ( 1 - epsi .* bed_mask );
    rho_NS = xx_0(3*nstages+1:4*nstages,:) ./ ( 1 - epsi .* bed_mask );
    u_NS   = xx_0(4*nstages+1:5*nstages,:);
    %sgtitle ('conservative')
    
else
    Cf_NS  = xx_0(0*nstages+1:1*nstages,:);
    rho_NS = xx_0(3*nstages+1:4*nstages,:);
    u_NS   = xx_0(4*nstages+1:5*nstages,:);

    %sgtitle ('nonconservative')
end

%u_NS   = Velocity(feedFlow,rho_NS, num2cell(Parameters));

P_NS   = Pressure_PR( T_NS, rho_NS, num2cell(Parameters) );

Res = {Cf_NS, Cs_NS, T_NS, rho_NS, u_NS, P_NS};

%%
ind = 1:numel(Time);

subplot(3,3,1); title('T')
hold on
plot(Time,feedTemp); 
plot(Time,xx_0(2*nstages+1,:))
plot(Time,xx_0(3*nstages  ,:))
hold off
legend('T_u','First Layer', 'Last Layer')

subplot(3,3,2)
plot(Time,feedPress); title('P_u')

subplot(3,3,3)
plot(Time,feedFlow)
ylim([0.95*min(feedFlow), 1.05* max(feedFlow)]); title('F_u')

NAME = {'C_f','C_s','T','Continuity [ (1-\epsilon) \times \rho ]','Momentum'};

for i=0:numel(NAME)-1

    subplot(3,3,i+4)
    
    imagesc(Time,1:nstages, Res{i+1}); colorbar;  colormap jet

    hold on
    yline(nstagesbed(1),'w')
    yline(nstagesbed(end),'w')
    hold off
    title(NAME{i+1})
end

subplot(3,3,numel(NAME)+4)
imagesc(Time,1:nstages,P_NS); title('P'); colorbar; colormap jet
hold on
yline(nstagesbed(1),'w')
yline(nstagesbed(end),'w')
hold off

%plotyy(Time, xx_0(Nx,:), Time, 1e3 * (sum(xx_0(0*nstages+1:1*nstages,ind) .* V_fluid) + sum(xx_0(1*nstages+1:2*nstages,ind) .* V_solid)) + xx_0(Nx,ind))

end