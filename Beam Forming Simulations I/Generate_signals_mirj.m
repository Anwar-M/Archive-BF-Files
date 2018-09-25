% Program to generate acoustic measurements on an array
% Check coherence and correlation

clear all
close all

% Sound speed (m/s)
c = 343;

% Microphone locations (nr, x, y)
load Array.txt
figure(1)
plot3(Array(:,1),Array(:,2),zeros(size(Array(:,1))),'k.')
axis square
h = get(gcf,'Children');
set(h,'FontSize',14)
xlabel('x [m]')
ylabel('y [m]')

%%

% Acoustic sources: number (N), locations (Loc), frequencies (F), source levels (S), phases (Ph)
N = 2;
Loc = [-3 -4  5;
        3  3  5];
F = [2*1000 2000];
S = [100 98];
Ph = [30*(pi/180) 0];
figure(1)
hold on
plot3(Loc(:,1),Loc(:,2),Loc(:,3),'g*')

%%

% signal length (sec)
T = 10;
% Sample frequency (Hz)
fs = 100e3;
% Number of samples used per time block
Nsam = 1e5;
% Construct time vector
t = 0:1/fs:T;

% Preallocate space for pm: array containing the pressure signals at the
% microphones per block and source
pm = zeros(N,size(Array,1),length(t));
% Generate signals in time domain
for ks = 1:N % loop over sources
    % Distances
    d = sqrt( (Array(:,1)-Loc(ks,1)).^2 + (Array(:,2)-Loc(ks,2)).^2 +Loc(ks,3).^2 );
    % Frequency of sound wave
    f = F(ks);
    % RMS pressure
    pe = 2e-5*10^(S(ks)/20);
    % Pressure amplitude of sound wave
    pmax = pe*sqrt(2);
    for km = 1:size(Array,1) % loop over microphones
        % Pressure signal at microphones
        pm(ks,km,:) = pmax/(4*pi*d(km))*cos(2*pi*f*(t - d(km)/c) + Ph(ks));
    end
end
% Sum over the pressure signals as originating from the different acoustic sources: Pm, containg the pressures at the
% microphones per microphone and time block
Pm = squeeze(sum(pm,1));
clear pm
% Plot over Nsam
figure(2)
plot(t(1:Nsam),Pm(:,1:Nsam))
h = get(gcf,'Children');
set(h,'FontSize',14)
xlabel('time [sec]')
ylabel('pressure amplitude [Pa]')

%% Covariances from time signal
for km = 1:size(Array,1)
    dummy = cov(Pm(1,:),Pm(km,:));
    C(km) = dummy(1,2);
end
figure(3)
plot(C,'.')
h = get(gcf,'Children');
set(h,'FontSize',14)
xlabel('Microphone')
ylabel('Cross covariance')

%% Fourier transforming the data over Nsam samples
Ni = floor(size(Pm,2)/Nsam);
for kt = 1:Ni % loop over time blocks
    p = Pm(:,1+(kt-1)*Nsam:kt*Nsam);
    p_four = fft(p');
    p_four = 2*p_four(1:Nsam/2,:);
    if kt == 1
        A  = 0.8*repmat([1:size(Array,1)]*max(max(abs(p_four))),size(p_four,1),1);
    end
    f_fft = (0:Nsam/2-1)*fs/Nsam;
    figure(4)
    plot(f_fft,abs(p_four)+A,'.-k')
    Fu = unique(F);
    for ks = 1:length(Fu) % loop over frequencies
        [dummy,ind] = min(abs(f_fft-Fu(ks)));
        display(['Evaluating frequency ' num2str(f_fft(ind)),' of time block ',num2str(kt)])
        P = transpose(p_four(ind,:)/Nsam);
        figure(4)
        hold on
        plot(f_fft(ind),abs(P*Nsam)+A(ind,:)','g*')
        Pf(kt,ks,:) = P;
        Ff(kt,ks) = f_fft(ind);
    end
end

%% Now matrix Pf is created containing Noblocks x Nofrequencies x No microphones Fourier Transformed Data points
% The corresponding frequencies are contained in matrix Ff
% Consider single block only
Ff = Ff(1,:);
Pf = squeeze(Pf(1,:,:));
if length(Ff) > 1
    Pf = Pf.';
end

% Beamforming: Define scanning grid
Ngrid = 1000; 
xmax  = 6;
ymax  = 6;
x     = linspace(-xmax,xmax,Ngrid+1);
y     = linspace(-ymax,ymax,Ngrid);
z     = 5;

for ks = 1:length(Ff) % loop over frequencies
    % Covariance matrix 
    R = Pf(:,ks)*Pf(:,ks)';
    for k1 = 1:numel(x)    % Loop over x
        for k2 = 1:numel(y) % Loop over y
            % Calculate distance from microphones to scan point
            r     = sqrt( (x(k1)-Array(:,1)).^2 + (y(k2)-Array(:,2)).^2 + z^2 );
            % Calculate time delays
            delay = r/c;
            % Create steering vector
            g = -exp(-1i*2*pi*Ff(ks)*delay)./(4*pi*r);
            % Beamsteering
            % S(k1,k2,ks) = 0.5*g'*R*g/(norm(g)^4);
            S(k1,k2,ks) = 0.5*g'*R*g/(norm(g)^2);
        end
    end
end

figure(5)
pcolor(x,y,squeeze(real(sum(S,3))'))
shading flat
hold on
plot(Loc(:,1),Loc(:,2),'w*')

return


% % Generate signals in frequency domain
% for ks = 1:N
%     for km = 1:size(Array,1)
%
%
%
%     end
% end


% End of program