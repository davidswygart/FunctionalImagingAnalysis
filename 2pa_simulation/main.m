%% setup
%universal
h = 6.62607015e-34; %Plank's constant (m^2*kg/s)
c = 299792458; %speed of light in a vacuum (m/s)
L = 6.0221408e23; %Avogadro's number (mol^-1)

%rods
sigma_rh = 50e-58; %2p cross-sectional fluorophore area (m^4*s/photon) **may be much larger for rhodopsin
lambda_max = 497e-9; % wavelength of maximal sensitivity of rhodopsin (m)
eca = 0.87e-12; % effective collecting area of a rod (m^2)
rhr = 7e7; %available rhodopsin molecules per rod, subject to bleaching (unitless)

%fluorophore
sigma_f = 35e-58; %see above. Value est. for Gcamp6f from Janelia @930nm (https://www.janelia.org/lab/harris-lab/research/photophysics/two-photon-fluorescent-probes)
conc_f = 5e-6; % (mol/l) - Yang et al., natcomm, 2018
%NB: Euler used an OGB1 conc. of 100-200uM

%objective
%Olympus LUMPlanFL N, 60x/1.00 W, IC/0/FN26.5
Ana = .6; %numerical aperture (unitless)
nrefr = 1.333; %refractive index of water (unitless)
d_fp_pos = 80e-6; %distance between focal plane and photoreceptor outer segments (m)
ps = 1e-18; %equivalent volume of activation (L). NB: 1um^3 = 1e-18 L
q = .2; %detection efficiency (unitless)


%laser
tp = 200e-15; %pulse duration (s)
fp = 80e6;%67167;%80e6; % pulse rate (Hz)
p0 = 6e-3; %laser power (W)
lambda = 928e-9; %wavelength (m)
dwell = 12.8e-6; %time per pixel (s)
ppl = 128; %pixels per line (unitless)
lpf = 128; %lines per frame (unitless)
blank = .1; %temporal blanking fraction (unitless)
flyback = 1e-3; %frame flyback time (s)
scanfield = [50 25]*1e-6; %scanning area (m)
beam_width = 1.2e-3; % the 1/e^2 value of the laser (m)
beam_mag = 3.5/60; %total magnification of the beam system (unitless)

%simulation
apron = 110e-6; %distance outside of scanfield to simulate (m)

%% stationary beam

ne = p0.^2 .* sigma_f ./ (tp.*fp.^2) .* (pi .* Ana.^2 ./ (h.*c.*lambda)).^2; %emitted photons per fluorophore per pulse
ne_t = ne .* fp .* conc_f .* ps .* L; %emitted photons per second
ne_t_det = ne_t .* asin(Ana/nrefr) ./pi .* q; %detected photons per second
ne_pix_det = ne_t .* dwell .* asin(Ana/nrefr) ./ pi .* q; %detected photons per pixel

Ar = tan(asin(Ana./nrefr)) .* d_fp_pos; %radius of outer segment activation (m)
Apos = pi .* Ar.^2; %area of outer segment activation (m^2)

fq_n = lambda_max ./ lambda; %relative frequency
S_laser = (6.565996913733051e+06.*exp(-18.*fq_n) + 0.121801754912951.*exp(-2.*fq_n)).^(-4); %relative sensitivity to laser (unitless)

equiv_1p = p0 * lambda ./ (h.*c) .* S_laser ./ Apos .* eca; % R*/rod/s

% uniform laser beam power is averaged over all outer segments and then
% passed through 2pa function

equiv_2p = sigma_rh .* (p0 ./ Apos).^2 ./ (tp .* fp) .* (lambda ./ (h.*c)) .^2 .* rhr; %R*/rod/s
% NB: there is a typo in the eyecup scope paper, xi = 1./(tp.*fp)
% See Denk W, Piston DW, Webb WW (1995), p 446
% Also note that the true inverse duty cycle for their system is 6.6e4, not
% 1e5

equiv_indirect = ne_t / (4*pi*d_fp_pos.^2) .* eca;
% NB: another typo in the equation for the surface area of a sphere

%% scanned beam
line_rate = 1./(dwell .* ceil(ppl ./ (1-blank)));
frame_rate = 1./(lpf./line_rate + flyback);
pix_delta = scanfield ./ [ppl lpf]; %m per step
dxy = scanfield./lcm(10.*ppl,10.*lpf); %resolution (m) -- at least 10x the scan resolution
sample_factor = pix_delta ./ dxy;

apron = ceil(apron./dxy) .* dxy;

% TODO: sparse?
[qx,qy] = meshgrid(-scanfield(1)/2 - apron(1): dxy(1): scanfield(1)/2 + apron(1), -scanfield(2)/2 - apron(2): dxy(2): scanfield(2)/2 + apron(2));
qa = zeros(size(qx));

%TODO: clean this up
[sx,sy] = meshgrid(apron(1)./dxy(1)+1 : sample_factor(1) : apron(1)./dxy(1) + ppl .* sample_factor(1), apron(2)./dxy(2)+1 : sample_factor(2) : apron(2)./dxy(2)+ lpf .* sample_factor(2));
qa(sub2ind(size(qa), sx, sy)) = 1; %scan locations

Arq = round(Ar./dxy) .* dxy;

[spreadx,spready] = meshgrid(-Arq(1):dxy(1):Arq(1), -Arq(2):dxy(2):Arq(2));%meshgrid((1:size(spread,1)) - size(spread,1)/2 -.5,(1:size(spread,2)) - size(spread,2)/2 - .5); 
% spread = double((spreadx.^2 + spready.^2) < Ar.^2);

% intensity = power / area
% I = 2 * p0 ./ pi ./ w.^2 * exp( -2 * r.^2 ./ w.^2) -- w is 1/e^2 of the
%   beam leaving the objective

% as the beam travels through the system, the 1/e^2 width changes according
% to the magnification (for a thin lens)
w = beam_width .* beam_mag;

% it remains to find p0. because the beam is restricted by the NA of the
% objective, and the theoretical beam is infinite, p0 should be higher
% than the measured power
p0inf = p0 ./ (1 - exp(- 2 * Ar.^2 ./ w.^2));

% gaussian
intensity = 2  .* p0inf ./ pi ./ w.^2 .* exp(-2 .* (spreadx.^2 + spready.^2) ./ w.^2) .* ((spreadx.^2 + spready.^2) < Ar.^2);

% gaussian annulus
% intensity = 2  .* p0inf ./ pi ./ w.^2 .* exp(-2 .* (Arq(1).^2 - (spreadx.^2 + spready.^2)) ./ w.^2) .* ((spreadx.^2 + spready.^2) < Ar.^2);

% TODO: cf. https://www.hindawi.com/journals/ijs/2017/7560141/

% uniform
% intensity = p0 ./ Apos .* ((spreadx.^2 + spready.^2) < Ar.^2);

spread = sigma_rh .* intensity.^2 ./ (tp .* fp) .* (lambda ./ (h.*c)) .^2 .* rhr;

qa = imfilter(qa, spread) .* dwell  .* frame_rate;

%% plot
cmax = max(spread,[],'all');

figure;clf;
subplot(221)
imagesc(qx(1,:)*1e6, qy(:,1)'*1e6, qa)
rectangle('position',[-scanfield./2,scanfield]*1e6,'edgecolor','k')
caxis([0 cmax]);
hcb = colorbar;
hcb.Label.String = 'R*/rod/s';
xlabel('Distance (\mum)')
ylabel('Distance (\mum)')
title('Scanned 2P Activation')

subplot(222)
imagesc(qx(1,:)*1e6, qy(:,1)'*1e6, spread)
caxis([0 cmax]);
hcb = colorbar;
hcb.Label.String = 'R*/rod/s';
xlabel('Distance (\mum)')
ylabel('Distance (\mum)')
title('Stationary 2P Activation')


subplot(223)
imagesc(qx(1,:)*1e6, qy(:,1)'*1e6, intensity)
hcb = colorbar;
hcb.Label.String = 'Intensity (J/s/um^2';
xlabel('Distance (\mum)')
ylabel('Distance (\mum)')
title('Stationary Beam Intensity')









