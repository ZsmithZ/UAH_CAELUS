%% Phased Array Test Plan – Example Data & Plots (Subarray → Full Array → Radome)
% Author: <Your Name> | Date: <today>
% This script synthesizes representative datasets and generates 9 figures
% (3 per test level) to illustrate what will be collected/assessed during:
%   1) Subarray tests
%   2) Full array tests
%   3) Integrated tests with radome
% All numbers are synthetic but physically plausible for demonstration.

clear; close all; clc; rng(7);

%% -------------------- Global RF/Array Parameters ------------------------
c0   = 299792458;            % speed of light [m/s]
f0   = 10e9;                 % center frequency [Hz] (X-band example)
BW   = 1.0e9;                % analysis bandwidth [Hz]
fvec = linspace(f0-0.5*BW, f0+0.5*BW, 121);
lam0 = c0/f0;

% Arrays
% Subarray (tile) - square 4x4
Nx_tile = 4; Ny_tile = 4; N_tile = Nx_tile*Ny_tile;
dx = 0.5*lam0; dy = 0.5*lam0;

% Full array - square 16x16 (scale as desired)
Nx = 16; Ny = 16; N = Nx*Ny;

% Scan/angle grids
theta_deg = linspace(-80,80,1601);       % cut at phi=0
theta_rad = deg2rad(theta_deg);

% Utility for plots
set(0,'defaultTextInterpreter','none');
set(0,'defaultAxesFontSize',11);
set(0,'defaultLineLineWidth',1.4);

%% ========================= 1) SUBARRAY TESTS ============================
% GOALS:
%  - Characterize S-parameter coupling magnitudes |S_ij| (mutual coupling map)
%  - Validate intra-tile beamforming accuracy & pointing
%  - Quantify calibration residuals (amplitude/phase)

% ---- 1.1 Synthetic Coupling Matrix (|S| in dB) at f0 --------------------
% Model: nearest-neighbor ~ -18 dB, dropping with inter-element distance
% Diagonal uses |S11| ~ -12 dB (return loss)
positions_tile = elementXY(Nx_tile, Ny_tile, dx, dy); % [x,y] per element
D_ij = squareform(pdist(positions_tile));             % inter-element distances [m]
Dlam = D_ij/lam0;                                     % in wavelengths

S_mag = zeros(N_tile); % linear magnitude
S11   = db2mag(-12);   % |S11|
for i=1:N_tile
    for j=1:N_tile
        if i==j
            S_mag(i,j) = S11;
        else
            % exponential decay vs distance; cap near -45 dB for far neighbors
            S_lin = db2mag(-18)*exp(-1.25*Dlam(i,j));
            S_mag(i,j) = max(S_lin, db2mag(-45));
        end
    end
end
S_db = mag2db(S_mag + eps);

figure('Name','(1) Subarray |S| Coupling Matrix','Color','w');
imagesc(S_db); axis image; colorbar;
title('(1) Subarray Coupling Magnitude |S| (dB) at f_0');
xlabel('Element Index'); ylabel('Element Index');
caxis([-45 -8]); colormap(parula);
text(1,1,sprintf('|S_{11}|≈%.1f dB',mag2db(S11)),'Color','w','FontWeight','bold');

% ---- 1.2 Intra-Tile Beamforming: Commanded vs Achieved Pointing ----------
theta_cmd_deg = 30;       % commanded scan (phi=0 cut)
theta_cmd = deg2rad(theta_cmd_deg);
k0 = 2*pi*f0/c0;

% Ideal weights for steering to theta_cmd at f0
w_ideal = tileWeights(Nx_tile, Ny_tile, dx, dy, k0, theta_cmd);

% Apply realistic residual errors (pre-calibration -> calibration -> residual)
amp_err_dB   = 0.6*randn(N_tile,1);      % initial amplitude error
ph_err_deg   = 6*randn(N_tile,1);        % initial phase error
% Assume calibration removes 80% of error:
amp_resid_dB = 0.2*amp_err_dB;
ph_resid_deg = 0.2*ph_err_deg;

w_meas = w_ideal .* db2mag(amp_resid_dB) .* exp(1j*deg2rad(ph_resid_deg));

% Compute tile patterns (normalized) over theta cut at f0 for both cases
AF_ideal = planarAF_cut(Nx_tile, Ny_tile, dx, dy, k0, theta_rad, w_ideal);
AF_meas  = planarAF_cut(Nx_tile, Ny_tile, dx, dy, k0, theta_rad, w_meas);

% Estimate pointing (theta of max)
[~,iIdeal] = max(AF_ideal);
[~,iMeas ] = max(AF_meas);
theta_pk_ideal = theta_deg(iIdeal);
theta_pk_meas  = theta_deg(iMeas);
pointing_err_deg = theta_pk_meas - theta_pk_ideal;

figure('Name','(2) Tile Beam: Commanded vs Achieved','Color','w');
plot(theta_deg,20*log10(AF_ideal+eps)); hold on; grid on;
plot(theta_deg,20*log10(AF_meas +eps),'--');
yline(-3,'k:','-3 dB');
xlabel('\theta (deg)'); ylabel('Normalized |AF| (dB)');
title(sprintf(['(2) Tile Beam Cut (\\phi=0): Commanded %d^\\circ vs Achieved ' ...
               '(Error: %.2f^\\circ)'],theta_cmd_deg,pointing_err_deg));
legend('Ideal (commanded)','Measured (residual errors)','Location','SouthWest');
ylim([-50 0]);

% ---- 1.3 Calibration Residuals: Amplitude & Phase Histograms ------------
figure('Name','(3) Calibration Residuals','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
histogram(amp_resid_dB, 'BinWidth',0.1,'FaceAlpha',0.8); grid on;
xlabel('Residual Amplitude Error (dB)'); ylabel('Count');
title('(3a) Amplitude Residuals');
nexttile;
histogram(ph_resid_deg, 'BinWidth',0.5,'FaceAlpha',0.8); grid on;
xlabel('Residual Phase Error (deg)'); ylabel('Count');
title('(3b) Phase Residuals');

%% ========================= 2) FULL ARRAY TESTS ==========================
% GOALS:
%  - Validate array-level beam quality across scan (SLL, BW, scan loss)
%  - Quantify EIRP variation vs scan & frequency (inc. beam squint)
%  - Assess spectral regrowth (PA nonlinearity) and compute ACPR

% ---- 2.1 Far-field Cuts vs Scan Angle with SLL Annotation ---------------
scan_list = [0 20 40];             % degrees
k0 = 2*pi*f0/c0;
w_full = cell(numel(scan_list),1);
AFs    = zeros(numel(scan_list), numel(theta_rad));
SLLs   = zeros(numel(scan_list),1);

for ii=1:numel(scan_list)
    th0 = deg2rad(scan_list(ii));
    w_full{ii} = tileWeights(Nx, Ny, dx, dy, k0, th0); % ideal steering (phi=0)
    AF = planarAF_cut(Nx, Ny, dx, dy, k0, theta_rad, w_full{ii});
    AFs(ii,:) = AF / max(AF); % normalize
    SLLs(ii)  = computeSLL_dB(theta_deg, 20*log10(AFs(ii,:)+eps), scan_list(ii));
end

figure('Name','(4) Full Array Far-field Cuts & SLL','Color','w');
plot(theta_deg, 20*log10(AFs(1,:)+eps)); hold on; grid on;
plot(theta_deg, 20*log10(AFs(2,:)+eps));
plot(theta_deg, 20*log10(AFs(3,:)+eps));
xlabel('\theta (deg)'); ylabel('Normalized |AF| (dB)');
title('(4) Full Array Normalized Cuts ( \phi=0 ) with SLL');
legend(arrayfun(@(a)sprintf('\\theta_0 = %d^\\circ, SLL≈%.1f dB',a,SLLs(scan_list==a)), ...
       scan_list,'UniformOutput',false),'Location','SouthWest');
ylim([-60 0]); yline(-3,'k:','-3 dB');

% ---- 2.2 EIRP vs Scan Angle at Multiple Frequencies (Scan Loss + Squint)-
% Model: EIRP_dB = Ptx_dBW + G_array(th,f) - feed/rad losses
Ptx_dBW = 40;           % total TX power [dBW] (e.g., 10 kW peak ~ 40 dBW)
feedLoss_dB = 1.5;      % fixed loss
freqs = f0*[0.95, 1.00, 1.05];

EIRP = zeros(numel(freqs), numel(theta_rad));
for k=1:numel(freqs)
    fk  = freqs(k);
    kk  = 2*pi*fk/c0;
    % Steering designed at f0 -> evaluate at fk (beam squint appears off-design)
    th0 = deg2rad(0:2:60); th0 = th0(:)';  % evaluate EIRP vs scan command
    EIRP_k = zeros(size(th0));
    for j=1:numel(th0)
        w0 = tileWeights(Nx, Ny, dx, dy, 2*pi*f0/c0, th0(j)); % design at f0
        % Gain proxy: |AF|^2 at intended look angle, but evaluated at fk
        AF_on = planarAF_cut(Nx, Ny, dx, dy, kk, th0(j), w0);
        Gnorm = (AF_on / max(AF_on)).^2; % normalized power gain at that angle
        % Projection factor (scan loss) ~ cos(theta)
        scanLoss = cos(th0(j));
        % Absolute gain scaling: Gmax ~ N^2 for coherent (normalized out); keep relative
        EIRP_k(j) = Ptx_dBW + mag2db(Gnorm) + 10*log10(max(scanLoss,1e-3)) - feedLoss_dB;
    end
    EIRP(k,1:numel(th0)) = EIRP_k;
    if k==1, th_grid = rad2deg(th0); end
end

figure('Name','(5) EIRP vs Scan & Frequency','Color','w');
plot(th_grid, EIRP(1,1:numel(th_grid))); hold on; grid on;
plot(th_grid, EIRP(2,1:numel(th_grid)));
plot(th_grid, EIRP(3,1:numel(th_grid)));
xlabel('Commanded Scan \theta_0 (deg)'); ylabel('EIRP (dBW) (relative)');
title('(5) EIRP vs Scan Angle with Beam Squint Across Frequency');
legend(arrayfun(@(ff)sprintf('f = %.2f GHz',ff/1e9), freqs,'UniformOutput',false), ...
       'Location','SouthWest');

% ---- 2.3 Spectral Regrowth & ACPR (PA Nonlinearity on OFDM) -------------
% Simple baseband OFDM, pass through memoryless cubic AM/AM+AM/PM nonlinearity.
Fs   = 100e6;    % sample rate
Nfft = 1024; Ncp = 128; Nused = 512; Ns = 800; % OFDM symbols
subIdx = (-Nused/2:Nused/2-1) + Nfft/2 + 1;     % centered used carriers

% Generate random QPSK on used subcarriers
Syms = (randi([0 1],Nused,Ns)*2-1) + 1j*(randi([0 1],Nused,Ns)*2-1);
X    = zeros(Nfft,Ns);
X(subIdx,:) = Syms;

% IFFT + CP
x = ifft(ifftshift(X,1),Nfft,1)*sqrt(Nfft);
xcp = [x(end-Ncp+1:end,:); x];
tx = xcp(:);
tx = tx/ rms(tx);   % normalize

% PA nonlinearity: y = a1*x + a3*|x|^2*x with AM/PM (phase term)
a1 = 1.0;
a3 = -0.25;           % cubic compression
phi = 10*pi/180;      % AM/PM (radians) proportional to |x|^2
y  = (a1*tx + a3*(abs(tx).^2).*tx) .* exp(1j*phi*(abs(tx).^2));

% PSD (periodogram)
NfftPSD = 4*2^nextpow2(length(tx));
TX = fftshift(fft(tx, NfftPSD));
TY = fftshift(fft(y , NfftPSD));
faxis = linspace(-Fs/2, Fs/2, NfftPSD);

PSD_tx = 20*log10(abs(TX)/sqrt(NfftPSD) + eps);
PSD_y  = 20*log10(abs(TY)/sqrt(NfftPSD) + eps);

% ACPR calculation: integrate power in main vs adjacent bands
Bmain = (Fs*Nused/Nfft);           % occupied BW
edge  = Bmain/2;
% Define masks for main and adjacent bands
mainMask = abs(faxis) <= edge;
adjMaskL = (faxis > -3*edge) & (faxis < -edge);
adjMaskR = (faxis >  edge)   & (faxis <  3*edge);

P_main = mean(abs(TY(mainMask)).^2);
P_adj  = 0.5*( mean(abs(TY(adjMaskL)).^2) + mean(abs(TY(adjMaskR)).^2) );
ACPR_dB = 10*log10(P_adj/P_main);

figure('Name','(6) Spectral Regrowth & ACPR','Color','w');
plot(faxis/1e6, PSD_tx); hold on; grid on;
plot(faxis/1e6, PSD_y, '--');
xlabel('Frequency (MHz)'); ylabel('PSD (dBFS/Hz) (relative)');
title(sprintf('(6) Spectrum: Linear vs PA-Nonlinear | ACPR ≈ %.1f dB', ACPR_dB));
legend('Before PA','After PA Nonlinearity','Location','SouthWest');
xlim([-2.5*edge/1e6, 2.5*edge/1e6]);

%% =================== 3) INTEGRATED TESTS WITH RADOME ====================
% GOALS:
%  - Quantify radome insertion loss (IL) vs angle & frequency
%  - Assess beam depointing (boresight error) and pattern distortion
%  - Evaluate polarization purity degradation (XPD) vs scan/freq

thetas = linspace(0,70,181);                    % incidence angle (deg)
fr     = linspace(0.9*f0, 1.1*f0, 61);         % ±10% band
[TH, FF] = meshgrid(deg2rad(thetas), fr);

% Simple parametric radome model (illustrative, not material-specific)
IL0_dB   = 0.6;    % on-boresight IL at center freq
k_theta  = 0.018;  % angle sensitivity
k_freq   = 6.0;    % frequency sensitivity
IL_dB = IL0_dB + k_theta*(rad2deg(TH)/60).^2 + k_freq*((FF-f0)/f0).^2;   % dB

figure('Name','(7) Radome Insertion Loss Map','Color','w');
imagesc(thetas, fr/1e9, IL_dB); axis xy; colorbar;
xlabel('Incidence Angle \theta (deg)'); ylabel('Frequency (GHz)');
title('(7) Radome Insertion Loss (dB)');
colormap(parula);

% Boresight error (wedge/phase aberration proxy): Δθ ≈ α * tan(θ) * freq term
alpha_bse = 0.45; % deg scaling
BSE = alpha_bse * tand(thetas) .* (1 + 0.2*((fr'-f0)/f0)); % deg, outer product

% Compare pattern bare vs radome at a challenging scan (e.g., 40 deg)
scan_eval_deg = 40; scan_eval = deg2rad(scan_eval_deg);
kk = 2*pi*f0/c0;

w_bare = tileWeights(Nx, Ny, dx, dy, kk, scan_eval);
AF_bare = planarAF_cut(Nx, Ny, dx, dy, kk, theta_rad, w_bare);
AF_bare = AF_bare / max(AF_bare);

% Apply radome effects at f0, approximate as:
%  - additional IL(θ) envelope
%  - depointing by BSE(θ0)
%  - mild phase distortion via broadened mainlobe
% Choose IL vs observation angle for plotting overlay
IL_line = interp1(thetas, IL_dB(fr==f0, :), theta_deg, 'linear', 'extrap');
IL_lin  = db2mag(-IL_line);  % attenuation (linear)

% De-pointed look angle:
bse_at_40 = interp1(thetas, BSE(:,thetas==thetas(round((scan_eval_deg/70)*numel(thetas)))) , 40, 'linear', 'extrap');
scan_eff = deg2rad(scan_eval_deg + bse_at_40); % shifted pointing

w_rado = tileWeights(Nx, Ny, dx, dy, kk, scan_eff);
AF_rado = planarAF_cut(Nx, Ny, dx, dy, kk, theta_rad, w_rado);

% Apply IL envelope to the pattern
AF_rado = AF_rado .* IL_lin(:).';
AF_rado = AF_rado / max(AF_rado);

figure('Name','(8) Bare vs Radome Pattern (De-pointing)','Color','w');
plot(theta_deg, 20*log10(AF_bare+eps)); hold on; grid on;
plot(theta_deg, 20*log10(AF_rado+eps), '--');
xlabel('\theta (deg)'); ylabel('Normalized |AF| (dB)');
title(sprintf('(8) Bare vs Radome at Commanded %d^\\circ (BSE applied)', scan_eval_deg));
legend('Bare Array','With Radome (IL + depointing)','Location','SouthWest');
ylim([-60 0]); yline(-3,'k:','-3 dB');

% XPD (Cross-Polar Discrimination) degradation model:
% Base XPD ~ 30 dB; degradations vs theta & frequency
XPD0 = 30; beta_th = 0.06; beta_f = 8.0;
XPD = XPD0 - beta_th*(thetas/60).^2 - beta_f*((fr'-f0)/f0).^2; % dB

% Plot XPD vs scan angle for three frequencies
fpick = [0.95, 1.00, 1.05]*f0;
figure('Name','(9) XPD vs Scan Angle','Color','w');
hold on; grid on;
for k=1:numel(fpick)
    [~, idxf] = min(abs(fr - fpick(k)));
    plot(thetas, XPD(idxf,:), 'DisplayName', sprintf('f = %.2f GHz', fpick(k)/1e9));
end
xlabel('Scan Angle \theta (deg)'); ylabel('XPD (dB)');
title('(9) Radome-Induced Polarization Purity Degradation');
legend('Location','SouthWest');

%% --------------------------- Helper Functions ---------------------------
function pos = elementXY(Nx, Ny, dx, dy)
% Return [x,y] positions (m) of an Nx-by-Ny planar array centered at (0,0)
[xi, yi] = meshgrid( (0:Nx-1) - (Nx-1)/2, (0:Ny-1) - (Ny-1)/2 );
pos = [xi(:)*dx, yi(:)*dy];
end

function w = tileWeights(Nx, Ny, dx, dy, k, theta0)
% Ideal complex weights to steer to (theta0, phi=0) at wavenumber k
pos = elementXY(Nx, Ny, dx, dy);
% phi=0 => only x-projection matters along the cut
phase = -k * ( pos(:,1)*sin(theta0) ); % y term zero at phi=0 cut
w = exp(1j*phase);
end

function AF = planarAF_cut(Nx, Ny, dx, dy, k, theta, w)
% Array factor along phi=0 cut
pos = elementXY(Nx, Ny, dx, dy);
% Steering/weights included via w
AF = zeros(size(theta));
for it=1:numel(theta)
    psi = k * ( pos(:,1)*sin(theta(it)) );
    AF(it) = abs(sum(w .* exp(1j*psi)));
end
AF = AF / max(AF+eps);
end

function SLLdB = computeSLL_dB(theta_deg, AFdB, scan_deg)
% Estimate SLL by excluding mainlobe region around the commanded scan angle
% and finding the highest sidelobe.
[~, idx_main] = max(AFdB);
% Define a guard around mainbeam (±8 deg)
guard = 8;
mask = true(size(theta_deg));
mask(theta_deg > theta_deg(idx_main)-guard & theta_deg < theta_deg(idx_main)+guard) = false;
SLLdB = max(AFdB(mask));
end
