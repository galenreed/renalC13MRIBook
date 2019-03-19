function [T2Spectra, fittedSpectra, signalFit] = t2nnls(signal, te, doPlot, T2Range, ...
                                       NT2Samples, regularizationType, regularizationLambda)
%
% T2 non-negative least squares algorithm for NMR relaxation  
% based on the reference: 
% Whittall and Mackay, J Magn Reson 84, 134-152 (1989)
% this implementation written by Galen Reed, 01/07/2019 



T2 = logspace(log10(T2Range(1)), log10(T2Range(2)), NT2Samples);
T2DecayMatrix = exp( -te(:) * (1./T2(:).') );
lambda = regularizationLambda;

if (regularizationType == 0)% no regularization, works great 
    
    fittedSpectra = lsqnonneg(T2DecayMatrix, signal);

elseif(regularizationType == 2)% Tikhonov regularization, works great
   
    lambdaDiag = lambda * eye(NT2Samples);
    zeroPad = zeros([NT2Samples 1]);
    Acat = vertcat(T2DecayMatrix, lambdaDiag);
    ycat = vertcat(signal, zeroPad);
    fittedSpectra = lsqnonneg(Acat, ycat);    
    
elseif(regularizationType == 1) % L1 not tested in forever
    
    [fittedSpectra,status] = l1_ls_nonneg(T2DecayMatrix, signal, lambda, rel_tol,quiet);
   
end


signalFit = T2DecayMatrix*fittedSpectra;
T2Spectra = fittedSpectra;

if(doPlot == 1)

    fontsize = 20;
    lw = 2;
    figure();
    semilogy(te, signalFit, '-', 'linewidth', lw, ...
    te, signal, 'x', 'linewidth', lw);
    grid on;
    xlim([0 te(end)]);
    xlabel('TE [s]');
    ylabel('signal');
    set(gca, 'fontsize', fontsize);

    
    figure();
    semilogx(T2, T2Spectra, 'linewidth', lw)
    grid on;
    xlim([T2(1) T2(end)]);
    xlabel('T2 [s]', 'fontsize', fontsize);
    ylabel('signal', 'fontsize', fontsize);
    set(gca, 'fontsize', fontsize);
end





