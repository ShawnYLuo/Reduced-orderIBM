clear;
addpath('../subaxis')
load Model_data/red_traj 
load Model_data/full_traj
vals = [1000 5000 10000];
for ax=1:3
chi1sr = chisr1{ax};
chi2sr = chisr2{ax};
% center positions of two cells at each time step (reduced model)
chi1xr = mean(chi1sr(1:end/2,:));
chi1yr = mean(chi1sr(end/2+1:end,:));
chi2xr = mean(chi2sr(1:end/2,:));
chi2yr = mean(chi2sr(end/2+1:end,:));


chi1sf = chisf1{ax};
chi2sf = chisf2{ax};
% center positions of two cells at each time step (full model)
chi1xf = mean(chi1sf(1:end/2,:));
chi1yf = mean(chi1sf(end/2+1:end,:));
chi2xf = mean(chi2sf(1:end/2,:));
chi2yf = mean(chi2sf(end/2+1:end,:));
dt = tspan(end)/(length(tspan)-1);

% periods of center positions of full and reduced models
highf = [];
highr = [];
for i=2:length(chi1yf)-1
    if chi1yf(i)>=chi1yf(i-1) && chi1yf(i)>=chi1yf(i+1)
        highf = [highf i*dt];
    end
    if chi1yr(i)>=chi1yr(i-1) && chi1yr(i)>=chi1yr(i+1)
        highr = [highr i*dt];
    end
end

figure(1)
subaxis(3,1,ax, 'sh', 0.01, 'sv', 0.01);
h(1) = plot(chi1xf,chi1yf+.2,'b','LineWidth',2);hold on
plot(chi2xf,chi2yf-.2,'b','LineWidth',2);
h(2) = plot(chi1xr,chi1yr+.2,'--r','LineWidth',2); 
plot(chi2xr,chi2yr-.2,'--r','LineWidth',2); 
axis([3.95 5.65 2 4])
set(gca,'XTick',[], 'YTick', [])
daspect([1 4 1])
ylabel(['s = -a = ',num2str(vals(ax))])
hold off
%legend(h([1,3]),{'full','reduced'})
% figure(2)
% hold on 
% plot(highr(2:end)-highr(1:end-1),'b');
% plot(highf(2:end)-highf(1:end-1),'r');
% legend('Reduced','Full');
% set(gca,'XTick',[], 'YTick', [])
% hold off
% diff = sqrt((chi1xr-chi1xf).^2 + (chi1yr-chi1yf).^2);
% figure(3)
% plot(tspan,diff)

buff_size = 0.5;
tracker = 0;
figure(2)
subaxis(3,1,ax, 'sh', 0.01, 'sv', 0.01, 'pb', 0);
axis([3.5 8.5 2 4])
set(gca,'XTick',[], 'YTick', [])
daspect([1 1 1])
ylabel(['s = -a = ',num2str(vals(ax))])
hold on
for i=1:floor(length(tspan)/5):length(tspan)
    chi1xr = chi1sr(1:end/2,i);
    chi1yr = chi1sr(end/2+1:end,i);
    chi1xf = chi1sf(1:end/2,i);
    chi1yf = chi1sf(end/2+1:end,i);
    chi2xr = chi2sr(1:end/2,i);
    chi2yr = chi2sr(end/2+1:end,i);
    chi2xf = chi2sf(1:end/2,i);
    chi2yf = chi2sf(end/2+1:end,i);
    h(3) = plot([chi1xf ; chi1xf(1)]+tracker*buff_size,[chi1yf ; chi1yf(1)],'b','LineWidth',2); 
    h(4) = plot([chi1xr ; chi1xr(1)]+tracker*buff_size,[chi1yr ; chi1yr(1)],'r--','LineWidth',2); 
    plot([chi2xf ; chi2xf(1)]+tracker*buff_size,[chi2yf ; chi2yf(1)],'b','LineWidth',2); 
    plot([chi2xr ; chi2xr(1)]+tracker*buff_size,[chi2yr ; chi2yr(1)],'r--','LineWidth',2);
    tracker = tracker + 1;
end
hold off
end
legend(h([1,2]),{'FOM','ROM'},'Orientation','horizontal')
legend(h([3,4]),{'FOM','ROM'},'Orientation','horizontal')
rmpath('../subaxis')