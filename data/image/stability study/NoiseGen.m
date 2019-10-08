clc
clear
close all

fname = 'Degrade.png';

i = imread(strcat('..\simulation\',fname));

for noise_std = 0.01:0.05:0.31
    noisy_i = imnoise(i,'gaussian',0,noise_std);
    imwrite(noisy_i,strcat(fname,'_gaussian_',num2str(noise_std),'.tif'));
end
subplot(1,2,1); imshow(i);
subplot(1,2,2); imshow(noisy_i);

for noise_density = 0.01:0.05:0.31
    noisy_i = imnoise(i,'salt & pepper',noise_density);
    imwrite(noisy_i,strcat(fname,'salt & pepper',num2str(noise_density),'.tif'));
end

% subplot(1,2,1); imshow(i);
% subplot(1,2,2); imshow(noisy_i);
