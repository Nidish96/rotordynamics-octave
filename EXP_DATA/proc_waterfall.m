clc
clear all
addpath('../ROUTINES/')

configs = {'ROD', 'LD'};
runs = {'CWUP', 'CWDOWN', 'CCWUP', 'CCWDOWN'};
speeds = [25 50 75 100 125 150 175 200];
fresp_xy = zeros(3, length(speeds));

config = configs{2};
run = runs{4};

## Load Data
load(sprintf('./DATS/%s_%s.mat', run, config), 'DAT')

## Waterfall Plots
figure(1)
clf()
for i=1:length(DAT)
  if DAT(i).found
    Nf = size(DAT(i).fdata,1);
    plot3(DAT(i).fdata(:,1), DAT(i).rpm*ones(Nf,1), abs(DAT(i).fdata(:,2)), 'k-', 'LineWidth', 2); hold on
    ## Interpolate for Fresp: Too Frigid
    ## fresp_xy(:, i) = interp1(DAT(i).fdata(:,1), DAT(i).fdata(:,2:3), DAT(i).rpm/60);

    ## ## "Loose" Interpolation
    ## [~, iP] = findpeaks(abs(DAT(i).fdata(:,2))+abs(DAT(i).fdata(:,3)));
    ## [mI,iM1] = min(abs(DAT(i).fdata(iP,1)-DAT(i).rpm/60)); iM1 = iP(iM1);

    ## ## Pick largest peak
    ## [~, iP] = findpeaks(abs(DAT(i).fdata(:,2))+abs(DAT(i).fdata(:,3)));
    ## [~, iM2] = max(abs(DAT(i).fdata(iP,2))+abs(DAT(i).fdata(iP,3))); iM2 = iP(iM2);

    ## if abs(DAT(i).fdata(iM2,1)/(DAT(i).rpm/60)-1)>0.1
    ##   iM = [iM1;iM2];
    ##   [~, ii] = min(abs(DAT(i).fdata(iM,1)-DAT(i).rpm/60)); iM = iM(ii);
    ## else
    ##   iM = iM2;
    ## end
    
    ## Search local peaks
    i1 = find(DAT(i).fdata(:,1)>DAT(i).rpm*0.9/60); i1 = i1(1);
    i2 = find(DAT(i).fdata(:,1)<DAT(i).rpm/60*1.1); i2 = i2(end);

    i12 = i1:i2;
    [~, iP] = findpeaks(abs(DAT(i).fdata(i12, 2)));  iP = i12(iP);
    [~, iM] = max(abs(DAT(i).fdata(iP, 2)));  iM = iP(iM);

    fresp_xy(:, i) = DAT(i).fdata(iM, :)';     

    plot3(fresp_xy(1, i), DAT(i).rpm, abs(fresp_xy(2, i)), 'ro', 'LineWidth', 2, 'MarkerFaceColor', 'r')
    ## plot3(DAT(i).fdata(i1:i2,1), DAT(i).rpm*ones(i2-i1+1,1), abs(DAT(i).fdata(i1:i2,2)), 'r-', 'LineWidth', 2); hold on
  else
    fresp_xy(:, i) = NaN;
  end
end
plot([0 speeds]./60, [0 speeds], 'k--', 'LineWidth', 2)
set(gca, 'fontsize', 20)
xlabel('Hz')
ylabel('RPM')
zlabel('Velocity Amplitude (m/s)')
title('X Whirl')
grid on

xlim([0 DAT(i).fdata(end,1)])
print(sprintf('./FIGS/%s_%s_X.eps',run,config), '-depsc');

xlim([0 200/60])
print(sprintf('./FIGS/%s_%s_X_ZOOM.eps',run,config), '-depsc');

## Waterfall Plots
figure(2)
clf()
for i=1:length(DAT)
  if DAT(i).found
    Nf = size(DAT(i).fdata,1);
    plot3(DAT(i).fdata(:,1), DAT(i).rpm*ones(Nf,1), abs(DAT(i).fdata(:,3)), 'k-', 'LineWidth', 2); hold on
    plot3(fresp_xy(1, i), DAT(i).rpm, abs(fresp_xy(3, i)), 'ro', 'LineWidth', 2, 'MarkerFaceColor', 'r')
  end
end
plot([0 speeds]./60, [0 speeds], 'k--', 'LineWidth', 2)
set(gca, 'fontsize', 20)
xlabel('Hz')
ylabel('RPM')
zlabel('Velocity Amplitude (m/s)')
title('Y Whirl')
grid on

xlim([0 DAT(i).fdata(end,1)])
print(sprintf('./FIGS/%s_%s_Y.eps',run,config), '-depsc');

xlim([0 200/60])
print(sprintf('./FIGS/%s_%s_Y_ZOOM.eps',run,config), '-depsc');


## Frequency Response
figure(3)
clf()
plot(fresp_xy(1,:)*60, abs(fresp_xy(2:3,:)), 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'w')
set(gca, 'fontsize', 20)

xlabel('Spin Speed (rpm)')
ylabel('Velocity Amplitude (m/s)')
legend('X Whirl', 'Y Whirl')

xlim([0 200])
print(sprintf('./FIGS/%s_%s_FRESP.eps',run,config), '-depsc')

save(sprintf('./DATS/FRESP_%s_%s.mat',run,config), 'speeds', 'fresp_xy', '-7')
