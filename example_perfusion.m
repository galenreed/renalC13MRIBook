clear all;
close all;
options = optimset('Display','off');

% 0 for regularized inverse, 1 for thresholded SVD
solutionType = 0;

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

AIFPixels = [136 120;
             136 117;
             138 116;
             139 116;
             142 116;
             106 139;
             109 135;
             120 118];
AIFHalfWindowSize = 2;



RES = zeros([Nx Ny Nz Nt]);
BF = zeros([Nx Ny Nz Nt]);
MFBF = zeros([Nx Ny Nz Nt]);

% save some representative traces 
sampleAIF = [];
sampleFit = [];
sampleTimeCourse = [];
sampleResidue = [];


%for kk=1:Nz
for kk=1:1
  
  %fill A matrix with an arterial input curve for each slice
  A = zeros([Nt Nt]);
  for mm = 1:Nt
    for nn = 1:mm
      % make a square ROI centered about the AIF
      ix = AIFPixels(kk, 1); 
      iy = AIFPixels(kk, 2);
      w = AIFHalfWindowSize;
      tempImage = squeeze(im(:, :, kk, mm + 1-nn));
      ROIPixels = tempImage((ix-w):(ix+w), (iy-w):(iy+w));
      A(mm, nn) = mean(ROIPixels(:));
    end
  end
  
  % svd followed by thresholding
  [U, S, V] = svd(A);
  diagS = diag(S);
  diagS(find( diagS < 0.1*max(max(diagS)))) = 0;
  nonZeroInds = find(diagS != 0);
  diagS(nonZeroInds) = 1./diagS(nonZeroInds);
  Sinv = diag(diagS);
  
  %compute perfusion (blood flow or BF) for each voxel
  for ii = 1:Nx 
    for jj = 1:Ny
      
      % b is a single pixel through time
      b = squeeze(im(ii,jj,kk,:));
      
      % solve for Ax=b
      if(solutionType == 0)
        x = inv( A' * A +regularizationLambda^2 * eye(Nt)' * eye(Nt)) *A' * b;
      elseif(solutionType == 1)
        x = V*Sinv*U'*b;
      end
      
      %residual function
      RES(ii,jj,kk,:) = x;
      
      %perfusion or blood flow (BF)
      BF(ii,jj,kk) = max(x);
      MFBF(ii,jj,kk) = sum(b) / sum(squeeze(A(:,1)));
      
      % grab a sample pixel from the renal cortex
      if((kk == 1) && (ii == 140) && (jj == 72))
        sampleAIF = A(:,1);
        sampleFit = A*x;
        sampleTimeCourse =  b;
        sampleResidual = x;
      end

    end
  end      
  progress = kk/Nz*100
end


figure();
imagesc(BF(:,:,1));
colormap gray

figure()
lw = 3;
fs = 20
plot(t, sampleAIF, '-', 'linewidth', lw, ...
     t, sampleTimeCourse, 'o', 'linewidth', lw,...
     t, sampleFit,'-','linewidth', lw);
legend('AIF', 'renal cortex pixel', 'fit');         
xlabel('time [s]', 'fontsize', fs);
ylabel('signal', 'fontsize', fs);

set(gca, 'fontsize', fs)


% next 2 lines just to make the legend with the right fontsize
copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 20);

figure();
plot(t, sampleResidual, '.-','linewidth', lw);
xlabel('time [s]', 'fontsize', fs);
ylabel('R(t) F', 'fontsize', fs);
set(gca, 'fontsize', fs)

figure()
subplot(2,1,1);
imagesc(BF(:,:,1,1)/timeBetweenImages*60, [0 10]);
set(gca, 'xtick',[],'ytick',[]);
colormap jet;
colorbar();


colorbar();
subplot(2,1,2);
imagesc(MFBF(:,:,1,1)/timeBetweenImages*60, [0 10]);
set(gca, 'xtick',[],'ytick',[]);
colormap jet;
colorbar();











