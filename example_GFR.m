clear all;
close all;
pkg load image;
options = optimset('Display','off');

% 0 for regularized inverse, 1 for thresholded SVD
solutionType = 1;

% regularization for least squares
regularizationLambda = 5e3;

load('data/ureaSSFPImages.mat');
im = double(im_int2);

Nx = size(im, 1);  
Ny = size(im, 2);  
Nt = size(im, 4);
Nz = size(im, 3);

timeBetweenImages = 6;
t = [0:1:(Nt-1)] * timeBetweenImages;

cortexCenter = [153 58];
medullaCenter = [152 76];            
halfWindowSize = 2;


% make a square ROI centered about the AIF
ix = cortexCenter(1); 
iy = cortexCenter(2);
w = halfWindowSize;
cortexROI = squeeze(im((ix-w):(ix+w), (iy-w):(iy+w), 1, :));
ix = medullaCenter(1); 
iy = medullaCenter(2);
medullaROI = squeeze(im((ix-w):(ix+w), (iy-w):(iy+w), 1, :));


cortex = zeros([1, Nt]);
medulla = zeros([1, Nt]);

for ii = 1:Nt
  thisTimePointROI = cortexROI(:, :, ii);
  cortex(ii) = mean(thisTimePointROI(:));
  thisTimePointROI = medullaROI(:, :, ii);
  medulla(ii) = mean(thisTimePointROI(:));
  
end

% make a montage of the zoomed in kidney
% montage doesnt work though
##xlimits = 100:200; ylimits = 40:110;
##for ii = 1:Nt
##  thisImage = squeeze(im(xlimits, ylimits, 1, ii));
##  figure();
##  imagesc(thisImage);
##  set(gca, 'xtick', [], 'ytick', []);
##  colormap gray;  
##end


uptakeInds = 2:5; % time indices of uptake
medullaVolume = 0.5; % a rough estimate, in mL or g
P = polyfit (t(uptakeInds), medulla(uptakeInds), 1);
fit = polyval(P, t(uptakeInds));
gfr = P(1) * medullaVolume / mean(cortex(uptakeInds)) * 60;%mL/min


fontsize = 20;
 lw = 2;
figure()
plot(t, cortex, 'o-', 'linewidth', lw,...
     t, medulla, 'o-','linewidth', lw,...
     t(uptakeInds), fit, '--','linewidth', 2*lw);
xlim([0 t(end)]);
legend('cortex ROI', 'medulla ROI', 'fit during uptake')
xlabel('time [s]');
ylabel('signal');
    set(gca, 'fontsize', fontsize);

    
    