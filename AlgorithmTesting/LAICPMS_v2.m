%% LA-ICP-MS Data reduction protocol
% For George's data

clear workspace;

%% 1. Input parameters for George's data

collectors.types = ['F'      'F'      'F'      'F'      'F'      'F'      'F' 'F'      'F'     'F'      'F'      'F'     'C' 'C' 'C' 'C'];
collectors.gains = [.9956997,1.041124,.9984599,1.005296,1.024909,1.017642,1,  1.034178,1.00102,.9846376,1.037016,.997009,1,   1,   1,   1];
%collectors.amplifierNoiseVariance = [1.6568*10^-9*ones(1,12) zeros(1,4)];
%collectors.countsPerVolt = [62415096.5*ones(1,12) 62415096.5*ones(1,4)];  %George's ion counter data reads in V assuming 10^11 resistor
collectors.amplifierNoiseVariance = [4.8568463e-9*ones(1,12) zeros(1,4)];
collectors.countsPerVolt = [2.080503e7*ones(1,12) 2.080503e7*ones(1,4)];
collectors.deadTimes = [zeros(1,12) [12,10.8,9.7,10]*10^-9];
collectors.deadTimeOneSigmas = [zeros(1,12) [1 1 1 1]*10^-9];
collectors.dataReportedAs = 'voltage'; %alternate is 'countsPerSecond'

n.isotopes = 7;
isotopes = struct('name',            {'Hg202' 'Pb204' 'Pb206' 'Pb207' 'Pb208' 'Th232' 'U238'}', ...
                  'collector',       { 15      13       12     11      10      2       1    }', ...
                  'integrationTime', { 1       1        1      1       1       1       1    }');

misc.r202Hg_204Hg = 4.346;

%misc.t = (0:n.p-1)';  %for now, integer times

%% 1. Input parameters for Jeff's data

% collectors.types = 'C';
% collectors.gains = 1;
% collectors.amplifierNoiseVariance = 0;
% collectors.countsPerVolt = 62415096.5;  %George's ion counter data reads in V assuming 10^11 resistor
% collectors.deadTimes = 12*10^-9;
% collectors.deadTimeOneSigmas = 0;% 1*10^-9;
% collectors.dataReportedAs = 'countsPerSecond'; %alternate is 'countsPerSecond'
% 
% n.isotopes = 8;
% isotopes = struct('name',            {'Hg202' 'Pb204' 'Pb206' 'Pb207' 'Pb208' 'Th232' 'U235' 'U238' }', ...
%                   'collector',       { 1      1       1        1       1       1       1      1     }', ...
%                   'integrationTime', { 0.012 0.012  0.012  0.012  0.012  0.012  0.012  0.012}');
% 
% misc.r202Hg_204Hg = 4.346;
% 
% %misc.t = (0:n.p-1)';  %for now, integer times

%% 1. Input parameters for Pieter's data

% collectors.types = 'C';
% collectors.gains = 1;
% collectors.amplifierNoiseVariance = 0;
% collectors.countsPerVolt = 62415096.5;  %George's ion counter data reads in V assuming 10^11 resistor
% collectors.deadTimes = 12*10^-9;
% collectors.deadTimeOneSigmas = 0;% 1*10^-9;
% collectors.dataReportedAs = 'countsPerSecond'; %alternate is 'voltage'
% 
% n.isotopes = 6;
% isotopes = struct('name',            {'Pb206' 'Pb207' 'Pb208' 'Th232' 'U235' 'U238' }', ...
%                   'collector',       { 1      1       1        1       1       1    }', ...
%                   'integrationTime', { 0.040  0.040   0.040    0.020   0.040   0.020}');

%misc.r202Hg_204Hg = 4.346;

%misc.t = (0:n.p-1)';  %for now, integer times


%%

choices.BLfit = 'mean';
choices.ablationFit = 'line';
%choices.ablationFit = 'exponential';

% figure out indices for isotopes
for i = 1:n.isotopes
    if     strcmp(isotopes(i).name, 'Pb206'), indx.Pb206 = i;
    elseif strcmp(isotopes(i).name, 'Pb207'), indx.Pb207 = i;
    elseif strcmp(isotopes(i).name, 'Pb208'), indx.Pb208 = i;
    elseif strcmp(isotopes(i).name, 'Th232'), indx.Th232 = i;
    elseif strcmp(isotopes(i).name, 'U235' ), indx.U235  = i;
    elseif strcmp(isotopes(i).name, 'U238' ), indx.U238  = i;
    end
end

%% Flags
flags.negativeBL = [];
flags.negativeOP = [];
flags.belowDetLim = [];

%% A. Read data

% filename = 'AQUIPE01-SL.txt';
% filename = 'PEIXE_1A.txt';
% filename = 'STDCZ1.txt';
filename = '112JBM13-SL1.txt';
fid = fopen(filename);

%intialize intensities structure
intensities(n.isotopes).back = 0;
intensities(n.isotopes).peak = 0;
intensities(n.isotopes).name = {' '};
intensities(n.isotopes).collectorNumber = 0;
intensities(n.isotopes).collector = ' ';

% burn first five lines
for i = 1:5, fgetl(fid); end 

% assign intensities
for i = 1:n.isotopes
    x = fgetl(fid); intensities(i).back = str2num(x(8:end-2)); %#ok<*ST2NM>
    x = fgetl(fid); intensities(i).peak = str2num(x(8:end-2));
    intensities(i).name = isotopes(i).name;
    intensities(i).collectorNumber = isotopes(i).collector;
    intensities(i).collector = collectors.types(intensities(i).collectorNumber);
    fgetl(fid); %burn last string
end

fclose(fid);

clear x fid i filename

%% 2. Building measured intensity data covariance matrices for each fraction

% n is # of background integrations, p is # of on-peak integrations
n.n = size(intensities(1).back,2); n.p = size(intensities(1).peak,2);
misc.t = (0:n.p-1)';  %for now, integer times... for use later

for i = 1:n.isotopes
    
    % If statement to separate collectors
    
    % If the collector is a 'Faraday' NOTE: Add multiplier for hidden division
    if strcmp(intensities(i).collector, 'F') 
        measuredIntensityFaraday = 3*[intensities(i).back intensities(i).peak]';
        isotopes(i).measuredIntensity = measuredIntensityFaraday * isotopes(i).integrationTime * ...
                                             collectors.countsPerVolt(intensities(i).collectorNumber);
        SiDiag = (measuredIntensityFaraday*collectors.countsPerVolt(intensities(i).collectorNumber)...
                  + collectors.amplifierNoiseVariance(i)*collectors.countsPerVolt(i)^2)/isotopes(i).integrationTime;
        
        %test for zero or negative variances
        if any(SiDiag(1:n.n)<=0)
            SiDiagBL = SiDiag(1:n.n);
            SiDiagBL(SiDiagBL<=0) = mean(SiDiagBL(SiDiagBL>0));
            SiDiag(1:n.n) = SiDiagBL;
            flags.negativeBL = [flags.negativeBL isotopes(i).name];
        end
        if any(SiDiag(n.n+1:end)<=0)
            SiDiagOP = SiDiag(n.n+1:end);
            SiDiagOP(SiDiagOP<=0) = mean(SiDiagOP(SiDiagOP>0));
            SiDiag(n.n+1:end) = SiDiagOP;
            flags.negativeOP = [flags.negativeOP isotopes(i).name];
        end
        
        isotopes(i).Si = diag(SiDiag);
        
        clear measuredIntensityFaraday measuredCountsFaraday SiDiag SiDiagBL SiDiagOP
        
    %if the collector is an 'Ion Counter'
    elseif strcmp(intensities(i).collector, 'C')
        
        %for George's data, reported as voltages
        if strcmp(collectors.dataReportedAs , 'voltage')
            measuredIntensityIonCounter = [intensities(i).back intensities(i).peak]';
            measuredIntensityIonCounter = measuredIntensityIonCounter * collectors.countsPerVolt(intensities(i).collectorNumber);
        %for Jeff's data, reported in cps
        elseif strcmp(collectors.dataReportedAs , 'countsPerSecond')  
            measuredIntensityIonCounter = [intensities(i).back intensities(i).peak]';
        else
            disp('bad collectors.dataReportedAs (voltage or cps)')
        end %if the data is reported as a 'voltage' or 'countsPerSecond'
        
        isotopes(i).measuredIntensity = measuredIntensityIonCounter;
        measuredVarianceFromIonCounts = measuredIntensityIonCounter / isotopes(i).integrationTime;
        SiDiag = measuredVarianceFromIonCounts;
        
        %test for negative variance terms
        if any(SiDiag(1:n.n)<=0)
            SiDiagBL = SiDiag(1:n.n);
            SiDiagBL(SiDiagBL<=0) = mean(SiDiagBL); %mean(SiDiagBL(SiDiagBL>0));
            SiDiag(1:n.n) = SiDiagBL;
            flags.negativeBL = [flags.negativeBL isotopes(i).name];
        end
        if any(SiDiag(n.n+1:end)<=0)
            SiDiagOP = SiDiag(n.n+1:end);
            SiDiagOP(SiDiagOP<=0) = mean(SiDiagOP); %mean(SiDiagOP(SiDiagOP>0));
            SiDiag(n.n+1:end) = SiDiagOP;
            flags.negativeOP = [flags.negativeOP isotopes(i).name];
        end
        
        Si = diag(SiDiag);
        
        measuredCountsIntensityIonCounterSquared = measuredIntensityIonCounter.^2;
        isotopes(i).Si = Si + collectors.deadTimeOneSigmas(intensities(i).collectorNumber)^2 * ...
            (measuredCountsIntensityIonCounterSquared*measuredCountsIntensityIonCounterSquared');
        
        clear measuredVarianceFromIonCounts SiDiag SiDiagBL SiDiagOP 
        clear measuredCountsIntensityIonCounterSquared i Si measuredIntensityIonCounter
        
    else
        
        disp('Collector is not identified as a Faraday or Ion Counter')
        
    end %end collector if statements
    
end %if data reported as voltages

clear i

%% 3. Isobaric Interference Correction
%Starting here, intensities are in cps



isotopes(2).measuredIntensity = isotopes(2).measuredIntensity - isotopes(1).measuredIntensity / misc.r202Hg_204Hg;

isotopes(2).Si = isotopes(2).Si + isotopes(1).Si / misc.r202Hg_204Hg^2;


%% 4. Model the baseline and extrapolate beneath the on-peak measurements


misc.onesVectorn = ones(n.n,1);
misc.onesVectorp = ones(n.p,1);

for i = 2:n.isotopes %skip Hg now
    
    if strcmp(choices.BLfit, 'mean') %For George's data, forced to be a mean
        
        % create Sib.  Recall: n.n is number of baseline integrations
        isotopes(i).Sib = isotopes(i).Si(1:n.n,1:n.n);
        
        %evaluate weighted mean
        isotopes(i).BL.meanIntensityVariance = 1 / (misc.onesVectorn' * (isotopes(i).Sib \ misc.onesVectorn) );
        isotopes(i).BL.meanIntensity = isotopes(i).BL.meanIntensityVariance * (misc.onesVectorn' * ...
            (isotopes(i).Sib \ isotopes(i).measuredIntensity(1:n.n) ));
        
        %Calculate chi-squared
        residualsBL = isotopes(i).measuredIntensity(1:n.n) - isotopes(i).BL.meanIntensity;
        isotopes(i).BL.X2 = residualsBL' * (isotopes(i).Sib \ residualsBL);
        
        isotopes(i).BL.extrapolatedOPvalues = isotopes(i).BL.meanIntensity * ones(n.p,1);
        
    %elseif strcmp(choices.BLfit, 'meanod')
        
        
        
    else disp('BL fitting function not coded yet')
    end %if choices.BLfit
    
    isotopes(i).measuredIntensityBLc = isotopes(i).measuredIntensity(n.n+1:end) - isotopes(i).BL.extrapolatedOPvalues;
    
    
end %for i = 2:n.isotopes

clear residualsBL i


%% 5. Propagate uncertainties in baseline-corrected on-peak measurements

for i = 2:n.isotopes
    
    if strcmp(choices.BLfit, 'mean')
        %J11 = (isotopes(i).BL.meanIntensityVariance * ( isotopes(i).Sib \ misc.onesVectorn ))';
        J11 = ones(1,n.n) / n.n;
        J21 = -ones(n.p,1);
        J22 = eye(n.p);
        
    else disp('BL fitting function not coded yet')
    end %if choices.BLfit
    
    JOnPeak1 = J21 * J11;
    JOnPeak = [JOnPeak1 J22];
    
    isotopes(i).Sopbc = JOnPeak * isotopes(i).Si * JOnPeak';
    
end %for i = 2:n.isotopes

clear J11 J21 J22 JOnPeak1 JOnPeak i

%% 6. Calculate raw ratios (not yet corrected for fractionation)
%% and their uncertainties (covariance matrices)

% A. Do some checks

% First, check that 204Pb and 207Pb are above the 'detection limit'
for i = [2 4]
    
    
    mean_opbc = mean(isotopes(i).measuredIntensityBLc);
    stdv_opbc = std(isotopes(i).measuredIntensityBLc);
    z_opbc = abs( (isotopes(i).measuredIntensityBLc - mean_opbc) ./ stdv_opbc);
    keep_opbc = isotopes(i).measuredIntensityBLc(z_opbc<=2);
    mean_opbck = mean(keep_opbc);
    stdv_opbck = std(keep_opbc);
    
    if mean_opbck <= 2*stdv_opbck 
        if     i == 2, flags.belowDetLim = [flags.belowDetLim '204Pb'];
        elseif i == 4, flags.belowDetLim = [flags.belowDetLim '207Pb'];
        end
    end
    
end

%Next, make sure there are no negative intensities
% and 'correct' the ones that are negative

for i = 2:n.isotopes
    
    isotopes(i).measuredIntensityBLc(isotopes(i).measuredIntensityBLc<0) = eps;
    
end

clear mean_opbc stdv_opbc z_opbc keep_opbc mean_opbck stdv_opbck

% B. Calculate the natural logarithm of the opbc intensities, convert cov matrix

for i = 2:n.isotopes
    
    isotopes(i).logIntensity = log( isotopes(i).measuredIntensityBLc );
    
    Jlogr = 1 ./ isotopes(i).measuredIntensityBLc;
    Jmat = Jlogr * Jlogr';
    isotopes(i).Sopbclr = Jmat .* isotopes(i).Sopbc;
    
    % temp:
    isotopes(i).Jlogr = Jlogr;
    isotopes(i).Jmat  = Jmat;
    
end

clear Jlogr Jmat

% C. Calculate log-ratios

n.ratios = 3;

%names
comp(1,1).name = 'r206Pb_207Pb'; indx.lr206207 = 1;
comp(2,1).name = 'r206Pb_238U';  indx.lr206238 = 2;
comp(3,1).name = 'r208Pb_232Th'; indx.lr208232 = 3;
%values
comp(1).logratios = isotopes(indx.Pb206).logIntensity - isotopes(indx.Pb207).logIntensity;
comp(2).logratios = isotopes(indx.Pb206).logIntensity - isotopes(indx.U238 ).logIntensity;
comp(3).logratios = isotopes(indx.Pb208).logIntensity - isotopes(indx.Th232).logIntensity;
% and covariance matrices
comp(1).Slr = isotopes(indx.Pb206).Sopbclr + isotopes(indx.Pb207).Sopbclr;
comp(2).Slr = isotopes(indx.Pb206).Sopbclr + isotopes(indx.U238 ).Sopbclr;
comp(3).Slr = isotopes(indx.Pb208).Sopbclr + isotopes(indx.Th232).Sopbclr;

for i = 1:n.ratios
    comp(i).ratios = exp(comp(i).logratios); %#ok<*SAGROW>
end

clear i 


%% 7. Calculate intercepts and uncertainties


if strcmp(choices.ablationFit, 'line')
    
    %variables to conform with Jim's old names
    a = [ones(n.p,1) (0:n.p-1)'];
    
    for i = 1:n.ratios
        
        b = comp(i).logratios;
        C = a' * (comp(i).Slr \ a);
        D = a' * (comp(i).Slr \ b);
        comp(i).fitM = C \ D;
        comp(i).Vp = inv(C);
        
        yhat = comp(i).fitM(1) + comp(i).fitM(2).*a(:,2);
        r = b - yhat;
        comp(i).X2 = r'*(comp(i).Slr \ r);
        comp(i).BIC = comp(i).X2 + 2*log(n.p); %BIC = chi^2 + p*log(n)
        
% plotting

figure
hold on
set(gca, 'XLim', [min(misc.t) max(misc.t)])
plot(misc.t, comp(i).logratios, '.')
for j = 1:n.p
    twosd = 2*sqrt(comp(i).Slr(j,j));
    line([misc.t(j) misc.t(j)], [comp(i).logratios(j)-twosd comp(i).logratios(j)+twosd], 'Color', 'b')
end
xlim_orig = get(gca, 'XLim'); ylim_orig = get(gca, 'YLim');

nInterp = 250;
tInterp = linspace(min(misc.t), max(misc.t), nInterp)';
env2s = zeros(nInterp,1);   %two-sigma uncertainties in y for unct envelope

Jp = [ones(nInterp,1) tInterp];
modelInterp = comp(i).fitM(1) + comp(i).fitM(2)*tInterp;

for j = 1:nInterp
    env2s(j) = 2*Jp(j,:)*comp(i).Vp*Jp(j,:)';
end

plot(tInterp, modelInterp, '-k')
plot(tInterp, modelInterp + env2s, '-g')
plot(tInterp, modelInterp - env2s, '-g')
set(gca, 'XLim', xlim_orig); set(gca, 'YLim', ylim_orig);
hold off
        
        
    end

    clear a b C D yhat r nInterp tInterp env2s modelInterp

%for exponential fit
elseif strcmp(choices.ablationFit, 'exponential')
    
%     maxiter = 1000;
%     chiTolerance = 0.000001;
%     
%     p0 = [1.234339605*10^-6 8.242848029151*10^-6 2.8234294350;
%          -0.138720389681    0.1009312690789     -2.314472542236;
%           0.306126290941    0.14663390466832    -2.5051004239176];
    
    for i = 1:n.ratios
        disp(comp(i).name)
%         [comp(i).expParams, comp(i).Vp, comp(i).X2ablfit] = ...
%             LevenbergMarquardt_Matrix_v1(misc.t, comp(i).logratios, comp(i).Slr, ...
%                 p0(i,:), maxiter, chiTolerance);
        [param0, ~, ~, ~, ~] = LevenbergMarquardt_MatrixOD_v9...
                            ('expfast', misc.t, comp(i).logratios, comp(i).Slr, [1 -1 0], 1000, 1e-10, 100);
        [comp(i).expParams, comp(i).Vp, comp(i).L, comp(i).MSWD, comp(i).BIC] = ...
            LevenbergMarquardt_MatrixOD_v9('expmat', misc.t, comp(i).logratios, comp(i).Slr, ...
                                                                         param0', 100, 1e-10, 100);

    comp(i).yint = comp(i).expParams(1)+comp(i).expParams(3);  %a+c
    comp(i).yintOneSigma = sqrt( [1 0 1]*comp(i).Vp*[1 0 1]' );
    %comp(i).BIC = comp(i).X2 + 3*log(n.p);
    
    hdl = figure(i);
    hold on
    plot(misc.t, comp(i).logratios)
    plot(misc.t, comp(i).expParams(1) .* exp(comp(i).expParams(2)*misc.t) + comp(i).expParams(3),'r')
    hold off
    
    end

else disp('Fit not yet supported')

end

clear i p0 maxiter chiTolerance

