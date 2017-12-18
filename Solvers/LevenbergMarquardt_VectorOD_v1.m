%% Vector version of Levenberg-Marquardt

function [params, Vp, L, MSWD, BIC] = LevenbergMarquardt_VectorOD_v1...
    (fitfun, t, y, v, p0, maxiter, chiTolerance, lambda0)

%fileID = fopen('outputSTD1.txt');

lambda = lambda0;
pod = p0';
n = length(t);
onesv = ones(n,1);

yhat = calc_yhat(fitfun, pod, t, onesv);

switch fitfun
    case 'expfast', m = 3; vod = 1; useod = false;
    case 'expmat',  m = 3; vod = v; useod = false;
    otherwise %for the overdispersion cases
        switch fitfun
            case 'meanod', m = 2;
            case 'lineod', m = 3;
            case 'expod',  m = 4;
        end
        useod = true; % inside 'otherwise', for overdispersion cases
        vod = v + pod(m); %inside 'otherwise' for od cases
end

r = y - yhat;
L = calcL(r,vod,fitfun);
disp(['  L0 = ' num2str(L)])

G = zeros(m,1);
H = zeros(m,m);
[G, H] = calc_GH(fitfun, t, r, vod, pod, onesv, G, H);

iters = 0;

% formatspec1 = '%s';
% formatspec2 = '%4.4f';
% strG = ['G' num2str(iters)];
% strH = ['H' num2str(iters)];
% fprintf(fileID, formatspec1, strG);
% fprintf(fileID, formatspec2, G);
% fprintf(fileID, formatspec1, strH);
% fprintf(fileID, formatspec2, H);

while iters < maxiter
    
    h = -(H + lambda*diag(diag(H))) \ G; %calculate update
    
%     strh = ['h' num2str(iter)];
%     fprintf(fileID, formatspec1, strh);
%     fprintf(fileID, formatspec2, h)
    
    podnew = pod + h;
    if strcmp(fitfun,'lineod'), podnew(3) = abs(podnew(3)); end %stabilize line fit
    
    yhat = calc_yhat(fitfun, podnew, t, onesv);
    %Sod_new = calcSod(S,podnew,useod,m,n); %update Sod for fits that use od
    vod_new = calcVod(v, podnew, useod, m); % vector analog of calcSod
    rnew = y - yhat; %test r
    Lnew = calcL(rnew,vod_new,fitfun);
    
    if Lnew > L  %if things got worse, try again
        
        lambda = lambda * 10;
        iters = iters + 1;
        
        disp(['Lnew = ', num2str(Lnew) ' rejected, ' num2str([pod' lambda])])
        
        %next up is top of while
        
    else % if things got better
        
        %disp(['Lnew = ' num2str(Lnew) ' accepted, ' num2str(pod')])
        
        if abs(1 - Lnew/L) < chiTolerance  %if things got better and solved
            
            params = podnew;
            
            if podnew(m)/v(1) > 10^-5; % if there is significant overdispersion
                Vp = inv(H);
            else           % if no overdispersion,
                % %REVERT TO VERSION OF THE SAME FIT FUNCTION WITHOUT OVERDISPERSION, QUIT HERE
                Vp = inv(H(1:m-1,1:m-1)); % avoid bad scaling from zero overdispersion NOTE TO JIM: IGNORE THIS LINE
            end
            L = Lnew;
           %MSWD = rnew'*(S\rnew) / (n-m);
            MSWD = sum(rnew.^2 ./ v) / (n-m);
            BIC = -2*L + m*log(n);
            disp([num2str(iters) ' iterations'])
            return
            
        else % if things got better but not solved, accept values, lambda /= 10
            
            lambda = lambda/10;
            pod = podnew;
            L = Lnew;
            r = rnew;
            vod = vod_new;
            iters = iters + 1;
            [G, H] = calc_GH(fitfun, t, r, vod, pod, onesv, G, H);
            
            disp(['Lnew = ' num2str(Lnew) ' accepted, ' num2str(pod')])
            
            %next up is top of while
            
        end %if solved or else just a step in the right direction
        
    end %if things are better/worse
    
end % while

%if you never get there in maxiters
if strcmp(fitfun, 'expfast') %if initializing for another fit
    params = pod; Vp = -999; L = -999; MSWD = -999; BIC = -999;
else
    disp('LM failure')
    params = -999; Vp = -999; L = -999; MSWD = -999; BIC = -999; % crash Noah's code
    % REVERT TO OVERDISPERSION CALCULATION FOR FIT FUNCTION WITH m-1 PARAMETERS
end %test for initialization usage

end % main function


%% Local Functions

function yhat = calc_yhat(fitfun, pod, t, onesv)
switch fitfun
    case 'meanod'
        yhat = pod(1) * onesv;
    case 'lineod'
        yhat = pod(1) + pod(2)*t;
    case {'expfast', 'expmat', 'expod'}
        yhat = pod(1) * exp(pod(2)*t) + pod(3);
end
end

function [G, H] = calc_GH(fitfun, t, r, vod, pod, onesv, G, H)

switch fitfun
    
    case {'expfast', 'expmat'} %for either non-od exp fit
        
        Ja =  exp(pod(2)*t);
        Jba = t.*Ja;  %arrayTimes
        Jb  = pod(1)*Jba;  %scalar multiplication
        Jbb = t.*Jb; %array times
        Jc =  onesv;
        Jabc = [Ja Jb Jc];
        
        if strcmp(fitfun, 'expfast') %if fitfun is 'expfast'
            
            G(1:3) = -Jabc'*r;
            
            H(1:3,1:3) = Jabc'*Jabc;
            H(2,2) = H(2,2) - r'*Jbb;
            H(1,2) = H(1,2) - r'*Jba;
            H(2,1) = H(1,2);
            
        elseif strcmp(fitfun,'expmat') %if fitfun is 'expmat'
            
            vodinvr = r./vod; % Sod \ r;  % note: Sod = S for 'expmat'
            vodinvJabc = Jabc./vod; % Sod \ Jabc;
            G(1:3) = -Jabc'*Sodinvr;
            H(1:3,1:3) = Jabc'*vodinvJabc;
            H(2,2) = H(2,2) - vodinvr'*Jbb;
            H(1,2) = H(1,2) - Sodinvr'*Jba;
            H(2,1) = H(1,2);
            
        end %the if statement that distiguishes expmat from expfast
        
    otherwise %use this branch of 'switch fitfun' for overdispersion fits
        
       %Sodinv = inv(Sod); %calculate Sodinv and Sodinvr for all od fits.
       %Sodinvr = Sod \ r;
        vodinv  = 1./vod;
        vodinvr = r./vod;
        
        switch fitfun %sort meanod, lineod, and expod
            
            case 'meanod'
                
                %Sodinv1 = Sod \ onesv;
                vodinv1 = 1 ./ vod;
                
                G(1) = -onesv' * vodinvr; % -onesv'*Sodinvr;
               %G(2) = -(1/2)*(trace(Sodinvr*Sodinvr') - trace(Sodinv));
                G(2) = -(1/2)*(vodinvr'*vodinvr - sum(vodinv));
                
                H(1,1) = onesv' * vodinv1;       % onesv'*Sodinv1;
                H(1,2) = vodinv1' * vodinvr; % trace(Sodinv1*Sodinvr');
                H(2,1) = H(1,2);
               %H(2,2) = trace(Sod \ (Sodinvr*Sodinvr')) - (1/2)*trace(Sodinv*Sodinv);
                H(2,2) = sum(vodinvr.^2 ./ vod) - (1/2)*sum(vodinv.^2);
                
            case 'lineod'
                
                %Sodinv1 = Sod \ onesv;
                vodinv1 = 1 ./ vod;
                
                G(1:2) = -[onesv t]'*vodinvr; % -[onesv t]'*Sodinvr;
               %G(3) = -(1/2)*(trace(Sodinvr*Sodinvr') - trace(Sodinv));
                G(3) = -(1/2)*(vodinvr'*vodinvr - sum(vodinv));
                
               %H(1:2,1:2) = [onesv t]'*(Sod \ [onesv t]);
                H(1:2,1:2) = [onesv./vod t./vod] * [onesv t];
               %H(1,3) = trace(Sodinv1*(Sodinvr'));   H(3,1) = H(1,3);
                H(1,3) = vodinv1'*vodinvr; H(3,1) = H(1,3);
               %H(2,3) = trace((Sod\t)*(Sodinvr'));  H(3,2) = H(2,3);
                H(2,3) = (t./vod)'*vodinvr; H(3,2) = H(2,3);
               %H(3,3) = trace(Sod \ (Sodinvr*Sodinvr'))- (1/2)*trace(Sodinv*Sodinv);
                H(3,3) = sum(vodinvr.^2 ./ vod) - (1/2)*sum(vodinv.^2);
                
            case 'expod'
                
                Ja =  exp(pod(2)*t);
                Jba = t.*Ja;  %arrayTimes
                Jb  = pod(1)*Jba;  %scalar multiplication
                Jbb = t.*Jb; %array times
                Jc =  onesv;
                Jabc = [Ja Jb Jc];
               %SodinvJabc = Sod \ Jabc;
                vodinvJabc = [Ja./vod Jb./vod Jc./vod];
               
                G(1:3) = -Jabc'*vodinvr; % -Jabc'*Sodinvr;
               %G(4) = -(1/2)*(Sodinvr'*Sodinvr - trace(Sodinv));
                G(4) = -(1/2)*(vodinvr'*vodinvr - sum(vodinv));
                
                H(1:3,1:3) = Jabc'*vodinvJabc;  % H(1:3,1:3) = Jabc'*SodinvJabc;
                H(2,2) = H(2,2) - vodinvr'*Jbb; % H(2,2) - Sodinvr'*Jbb;
                H(1,2) = H(1,2) - vodinvr'*Jba; % H(1,2) - Sodinvr'*Jba;
                H(2,1) = H(1,2);
                H(1:3,4) = vodinvJabc'*vodinvr; % SodinvJabc'*Sodinvr;
                H(4,1:3) = H(1:3,4)';
               %H(4,4) = trace(Sod \ (Sodinvr*Sodinvr')) - (1/2)*trace(Sodinv*Sodinv);
                H(4,4) = sum(sodinvr.^2 ./ vod) - (1/2)*(vodinv'*vodinv);
               
        end %switch fitfun for overdispersion cases
        
end %switch fitfun for all cases

end

function L = calcL(r,vod,fitfun)
switch fitfun
    case 'expfast'
        L = r'*r;     % note: L = SSE
    case 'expmat'
        L = sum(r.^2 ./ vod); % r'*(Sod\r); note: L = X2
    otherwise %for overdispersion cases
        L = 1/2*(sum(r.^2 ./ vod) + log(prod(vod))); % 1/2*(r'*(Sod\r) + log(det(Sod)));
end
end

function vod = calcVod(v,pod,useod,m)% vector analog of calcSod, no n required

if useod, vod = v + pod(m); % Sod = S + pod(m)*eye(n);
else vod = v; % Sod = S;
end
%note: for expfast, S is not used, and for expmat, Sod = S and is not updated.
end

