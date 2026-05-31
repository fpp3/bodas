% Drawing asymptotic bode
%
% Auralius Manurung
% manurung.auralius@gmail.com
% 
% MODIFICADO para normalizar ceros y polos a 0 dB en baja frecuencia.
%

function [G, w] = bodas(sys)
[z, p, k] = zpkdata(sys);
z = z{1};
p = p{1};

i = 1;
while(~isempty(z))
    if isreal(z(i)) == false
        z(find(z == conj(z(i)))) = [];
    end
    z(i) = -z(i);
    i = i + 1;
    if i > length(z)
        break;
    end
end

i = 1;
while(~isempty(p))
    if isreal(p(i)) == false
        p(find(p == conj(p(i)))) = [];
    end
        p(i) = -p(i);

    i = i + 1;
    if i > length(p)
        break;
    end
end

% -------------------------------------------------------------------------
% NUEVO: Calcular la Ganancia Normalizada (K_norm)
% -------------------------------------------------------------------------
K_norm = k;
for j = 1:length(z)
    if z(j) ~= 0
        if isreal(z(j))
            K_norm = K_norm * z(j);
        else
            wn = sqrt(z(j)*conj(z(j)));
            K_norm = K_norm * (wn^2);
        end
    end
end
for j = 1:length(p)
    if p(j) ~= 0
        if isreal(p(j))
            K_norm = K_norm / p(j);
        else
            wn = sqrt(p(j)*conj(p(j)));
            K_norm = K_norm / (wn^2);
        end
    end
end
% -------------------------------------------------------------------------

z = sort(z);
p = sort(p);

W = unique(nonzeros([z; p]));
for j  = 1 : length(W)
    if isreal(W(j)) == false
        W(j) = ceil(log10(sqrt(W(j) * conj(W(j)))));
    else
        W(j) = ceil(log10(W(j)));
    end
end

W = sort(W);

if isempty(W)
    W = 1;
end

wmin = W(1)-2;
wmax = W(end)+2;
omega = logspace(wmin, wmax, 5000);

s = tf('s');

% -------------------------------------------------------------------------
% Data calculation
% -------------------------------------------------------------------------

% Gain data
A_gain = zeros(length(z), length(omega));
B_gain = zeros(length(p), length(omega));
C_gain = zeros(1, length(omega));
legend_text_gain = {};

G = 1;
for i = 1:length(z)
    if z(i) == 0 % Zero at the origin
        A_gain(i,:) = 20*log10(omega);
        legend_text_gain{end+1} = "s";
        G = G * s;
    elseif isreal(z(i)) == true
        for j = 1:length(omega)
            if omega(j) >= z(i)
                A_gain(i,j) = 20*log10(omega(j)/z(i));
            else
                A_gain(i,j) = 0;
            end
        end
        legend_text_gain{end+1} = ["(1 + s/", num2str(z(i)), ")"];
        G = G * (1 + s/z(i));
    else
        peak = false;
        wn = sqrt(z(i)*conj(z(i)));
        zeta = (z(i)+conj(z(i)))/2/wn;
        for j = 1:length(omega)
            if omega(j) >= wn
                A_gain(i,j) = 40*log10(omega(j)/wn);
                if peak == false
                    peak_idx = j;
                    peak = true;
                end
            else
                A_gain(i,j) = 0;
            end
        end
        if zeta < 0.5 && peak == true
            A_gain(i,peak_idx) = A_gain(i,peak_idx) - 20*log10(1/(2*zeta));
        end
        legend_text_gain{end+1} = ["(1 + s/", num2str(z(i)), ")(1 + s/", num2str(conj(z(i))), ")"];
        G = G * (1 + s/z(i)) * (1 + s/conj(z(i)));
    end
end

for i = 1:length(p)
    if p(i) == 0  % Pole at the origin
        B_gain(i,:) = -20*log10(omega);
        legend_text_gain{end+1} = "1/s";
        G = G * 1/s;
    elseif isreal(p(i)) == true
        for j = 1:length(omega)
            if omega(j) >= p(i)
                B_gain(i,j) = -20*log10(omega(j)/p(i));
            else
                B_gain(i,j) = 0;
            end
        end
        legend_text_gain{end+1} = ["1/(1 + s/", num2str(p(i)), ")"];
        G = G * 1/(1 + s/p(i));
    else
        peak = false;
        wn = sqrt(p(i)*conj(p(i)));
        zeta = (p(i)+conj(p(i)))/2/wn;
        for j = 1:length(omega)
            if omega(j) >= wn
                B_gain(i,j) = -40*log10(omega(j)/wn);
                if peak == false
                    peak_idx = j;
                    peak = true;
                end
            else
                B_gain(i,j) = 0;
            end
        end
        if zeta < 0.5 && peak == true
            B_gain(i,peak_idx) = B_gain(i,peak_idx) + 20*log10(1/(2*zeta));
        end
        legend_text_gain{end+1} = ["1/( (1+s/", num2str(p(i)), ")(1+s/", num2str(conj(p(i))), ") )"];
        G = G * 1/( (1+s/p(i)) * (1+s/conj(p(i))) );
    end
end

if K_norm ~= 0
    kDb = 20*log10(abs(K_norm));
    C_gain(:) = kDb;
    legend_text_gain{end+1} = ['K = ', num2str(kDb), ' dB'];
    G = G * K_norm;
end

[mag, phase, wout] = bode(G, {10^wmin, 10^wmax});
phase = squeeze(phase);
mag = squeeze(mag);
M = sum([sum(A_gain,1); sum(B_gain,1); C_gain], 1);

% Phase data
A_phase = zeros(length(z), length(omega));
B_phase = zeros(length(p), length(omega));
C_phase = zeros(1, length(omega));
legend_text_phase = {};

for i = 1:length(z)
    if z(i) == 0
        A_phase(i,:) = 90;
        legend_text_phase{end+1} = "s";
    elseif isreal(z(i)) == true
        for j = 1:length(omega)
            if j > 1 && omega(j) > z(i)/10 && omega(j) <= z(i)*10
                A_phase(i,j) = min(45*log10(omega(j)/omega(j-1)) + A_phase(i,j-1), 90);
            elseif omega(j) <= z(i)/10
                A_phase(i,j) = 0;
            elseif omega(j) > z(i)*10
                A_phase(i,j) = 90;
            end
        end
        legend_text_phase{end+1} = ["(1 + s/", num2str(z(i)), ")"];
    else
        wn = sqrt(z(i)*conj(z(i))); 
        zeta = (z(i)+conj(z(i)))/2/wn;
        delta = 10^zeta;
        for j = 1:length(omega)
            if j > 1 && (omega(j) > wn/delta) && (omega(j) <= wn*delta)
                A_phase(i,j) = 90/zeta*log10(omega(j)/omega(j-1)) + A_phase(i,j-1);
            elseif omega(j) <= wn/delta
                A_phase(i,j) = 0;
            elseif omega(j) > wn*delta
                A_phase(i,j) = 180;
            end
        end
        legend_text_phase{end+1} = ["(1 + s/", num2str(z(i)), ")(1 + s/", num2str(conj(z(i))), ")"];
    end
end

for i = 1:length(p)
    if p(i) == 0
        B_phase(i,:) = -90;
        legend_text_phase{end+1} = "1/s";
    elseif isreal(p(i)) == true
        for j = 1:length(omega)
            if j > 1 && omega(j) > p(i)/10 && omega(j) <= p(i)*10
                B_phase(i,j) = max(-45*log10(omega(j)/omega(j-1)) + B_phase(i,j-1), -90);
            elseif omega(j) <= p(i)/10
                B_phase(i,j) = 0;
            elseif omega(j) > p(i)*10
                B_phase(i,j) = -90;
            end
        end
        legend_text_phase{end+1} = ["1/(1 + s/", num2str(p(i)), ")"];
    else
        wn = sqrt(p(i)*conj(p(i))); 
        zeta = (p(i)+conj(p(i)))/2/wn;
        delta = 10^zeta;
        for j = 1:length(omega)
            if j > 1 && (omega(j) > wn/delta) && (omega(j) <= wn*delta)
                B_phase(i,j) = -90/zeta*log10(omega(j)/omega(j-1)) + B_phase(i,j-1);
            elseif omega(j) <= wn/delta
                B_phase(i,j) = 0;
            elseif omega(j) > wn*delta
                B_phase(i,j) = -180;
            end
        end
        legend_text_phase{end+1} = ["1/( (1+s/", num2str(p(i)), ")(1+s/", num2str(conj(p(i))), ") )"];
    end
end

if K_norm ~= 0
    if K_norm > 0
        C_phase(:) = 0;
    else
        C_phase(:) = -180;
    end
    legend_text_phase{end+1} = ['K = ', num2str(K_norm)];
end

S = sum([sum(A_phase,1); sum(B_phase,1); C_phase], 1);
dphase = S(end) - phase(end);

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------

% Figure 1: Gain with parts
figure('Name', 'Magnitude: Asymptotic with Parts', 'NumberTitle', 'off');
set(gcf,'units','normalized','outerposition',[0 0.5 0.5 0.5])
hold on; grid on;
for i = 1:size(A_gain, 1), plot(omega, A_gain(i,:), 'LineWidth', 1.5); end
for i = 1:size(B_gain, 1), plot(omega, B_gain(i,:), 'LineWidth', 1.5); end
if K_norm ~= 0, plot(omega, C_gain, 'LineWidth', 1.5); end
plot(omega, M, 'LineWidth', 3, 'Color', [.7, .7, .7]);
plot(wout, 20*log10(mag), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
legend([legend_text_gain, "Asymptotic Sum", "Exact Bode"], 'Location', 'best');
ylabel('Magnitude (dB)'); xlabel('Frequency (rad/s)'); set(gca, 'XScale', 'log');

% Figure 2: Gain only
figure('Name', 'Magnitude: Asymptotic only', 'NumberTitle', 'off');
set(gcf,'units','normalized','outerposition',[0.5 0.5 0.5 0.5])
hold on; grid on;
plot(omega, M, 'LineWidth', 3, 'Color', [.7, .7, .7]);
plot(wout, 20*log10(mag), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
legend({"Asymptotic Sum", "Exact Bode"}, 'Location', 'best');
ylabel('Magnitude (dB)'); xlabel('Frequency (rad/s)'); set(gca, 'XScale', 'log');

% Figure 3: Phase with parts
figure('Name', 'Phase: Asymptotic with Parts', 'NumberTitle', 'off');
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5])
hold on; grid on;
for i = 1:size(A_phase, 1), plot(omega, A_phase(i,:), 'LineWidth', 1.5); end
for i = 1:size(B_phase, 1), plot(omega, B_phase(i,:), 'LineWidth', 1.5); end
if K_norm ~= 0, plot(omega, C_phase, 'LineWidth', 1.5); end
plot(omega, S, 'LineWidth', 3, 'Color', [.7, .7, .7]);
plot(wout, phase + dphase, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
legend([legend_text_phase, "Asymptotic Sum", "Exact Bode"], 'Location', 'best');
ylabel('Phase (degrees)'); xlabel('Frequency (rad/s)'); set(gca, 'XScale', 'log');

% Figure 4: Phase only
figure('Name', 'Phase: Asymptotic only', 'NumberTitle', 'off');
set(gcf,'units','normalized','outerposition',[0.5 0 0.5 0.5])
hold on; grid on;
plot(omega, S, 'LineWidth', 3, 'Color', [.7, .7, .7]);
plot(wout, phase + dphase, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k');
legend({"Asymptotic Sum", "Exact Bode"}, 'Location', 'best');
ylabel('Phase (degrees)'); xlabel('Frequency (rad/s)'); set(gca, 'XScale', 'log');

w = {10^wmin 10^wmax};
end