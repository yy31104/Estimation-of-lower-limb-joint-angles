clearvars; close all;
load('typical_gait.mat');

%calculate based on data from typical_gait.mat (without error)
z_LT = (lfeo - ltio) ./ vecnorm((lfeo - ltio), 2, 2);
z_RF = (rfep - rfeo) ./ vecnorm((rfep - rfeo), 2, 2);
y_RF = (rfeo - rkne) ./ vecnorm((rfeo - rkne), 2, 2);
y_LF = (lkne - lfeo) ./ vecnorm((lkne - lfeo), 2, 2);
x_LF = cross(y_RF, z_RF) ./ vecnorm(cross(y_RF, z_RF), 2, 2);
beta_LK = asin(dot(z_LT, y_LF, 2));
alpha_LK = asin(dot(z_LT, x_LF, 2) ./ cos(beta_LK));

%graph (without error)
figure;
plot(time, alpha_LK, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Left Knee Flexion Angle (rad)');
title('Left Knee Flexion Angle Over Time', 'FontSize', 14);

%add noise
%randn - generate random error following normal distribution
%repeat this part?
noise_std = 40;
lfeo_noisy = lfeo + noise_std * randn(size(lfeo));
ltio_noisy = ltio + noise_std * randn(size(ltio));
rfep_noisy = rfep + noise_std * randn(size(rfep));
rfeo_noisy = rfeo + noise_std * randn(size(rfeo));
rkne_noisy = rkne + noise_std * randn(size(rkne));
lkne_noisy = lkne + noise_std * randn(size(lkne));

%smooth datas with filter (using different fc)
fc_values = [0.1, 0.5, 1, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20];
MAE_values = zeros(size(fc_values));

set(0, 'DefaultFigureWindowStyle', 'normal')

for i = 1 : length(fc_values)
    fc = fc_values(i);

    lfeo_filt = lpfilt(lfeo_noisy, fc);
    ltio_filt = lpfilt(ltio_noisy, fc);
    rfep_filt = lpfilt(rfep_noisy, fc);
    rfeo_filt = lpfilt(rfeo_noisy, fc);
    rkne_filt = lpfilt(rkne_noisy, fc);
    lkne_filt = lpfilt(lkne_noisy, fc);

    z_LT_filt = (lfeo_filt - ltio_filt) ./ vecnorm((lfeo_filt - ltio_filt), 2, 2);
    z_RF_filt = (rfep_filt - rfeo_filt) ./ vecnorm((rfep_filt - rfeo_filt), 2, 2);
    y_RF_filt = (rfeo_filt - rkne_filt) ./ vecnorm((rfeo_filt - rkne_filt), 2, 2);
    y_LF_filt = (lkne_filt - lfeo_filt) ./ vecnorm((lkne_filt - lfeo_filt), 2, 2);
    x_LF_filt = cross(y_RF_filt, z_RF_filt) ./ vecnorm(cross(y_RF_filt, z_RF_filt), 2, 2);
    beta_LK_filt = asin(dot(z_LT_filt, y_LF_filt, 2));
    val = dot(z_LT_filt, x_LF_filt, 2) ./ cos(beta_LK_filt);
    val = max(min(val, 1), -1); 
    alpha_LK_filt = asin(val);

    figure(i+1);
    plot(time, alpha_LK_filt, 'b', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel(['Filtered Left Knee Flexion Angle(rad),' num2str(fc)]);
    title(['Filtered Left Knee Flexion Angle Over Time at fc = ', num2str(fc)], 'FontSize', 14);

    %calculate mean absolute error
    MAE_values(i) = mean(abs(alpha_LK - alpha_LK_filt));
    fprintf('fc = %.1f Hz, MAE = %.4f rad\n', fc, MAE_values(i));

end

%plot MAE at different frequencies on graph
figure;
plot(fc_values, MAE_values, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Cutoff Frequency fc (Hz)');
ylabel('Mean Absolute Error (rad)')
title('Mean Absolute Error at Different Cutoff Frequencies', 'FontSize', 14);

