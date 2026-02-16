function nascarCS()
% Simulació de sonido 3D tipo carrera NASCAR usando HRTF sintéticas.

    close all; 
    clc;

    % CARGA Y PREPARACIÓN DEL AUDIO--------------------------------%

    [audio, Fs] = audioread('nascar-race.wav'); % Lee el archivo de audio y obtiene la señal de audio y la frecuencia de muestreo (Hz)
    audio = mean(audio,2); % Audio a mono 
    audio = audio / (max(abs(audio)) + eps); % Normaliza el audio para evitar saturación y eps evita división por cero

    duracion_simulacion = 15; % Duración total de la simulación en segundos

    bloque_size = 4096; % Tamaño de cada bloque de audio a procesar

    hop = floor(bloque_size/2); % Desplazamiento entre bloques (50% de solapamiento)

    n_muestras_total = round(duracion_simulacion * Fs); % Número total de muestras que tendrá la simulación

    % LOOP DE LA SEÑAL---------------------------------------------%

    n_repeticiones = ceil(n_muestras_total / length(audio)); % Calcula cuántas veces hay que repetir el audio original

    audio_source = repmat(audio, n_repeticiones, 1); % Repite el audio las veces necesarias

    audio_source = audio_source(1:n_muestras_total); % Recorta el audio a la longitud exacta de la simulación

    % FADE GLOBAL--------------------------------------------------%

    w_global = tukeywin(length(audio_source), 0.05); % Crea una ventana Tukey (Hann) para suavizar inicio y final

    audio_source = audio_source .* w_global; % Aplica el fade-in y fade-out al audio

    % HRIRs SINTÉTICAS (AZIMUT)------------------------------------%

    angulos_az = 0:10:350; % Define ángulos alrededor del oyente cada 10 grados

    n_ang = length(angulos_az); % Número total de ángulos

    hrtf_length = 128; % Longitud de cada HRIR en muestras

    HRIRs = zeros(n_ang, 2, hrtf_length); % Matriz para guardar las HRIRs: (ángulo, canal izquierdo/derecho, tiempo)

    delay_max_ms = 0.6; % Retardo máximo entre oídos en milisegundos

    delay_max = round((delay_max_ms/1000) * Fs); % Retardo máximo convertido a muestras

    % GENERACIÓN DE HRIRs------------------------------------------%

    for i = 1:n_ang % Bucle que recorre todos los ángulos

        azi = angulos_az(i); % Ángulo actual

        % Convierte el ángulo al rango [-180, 180]
        if azi > 180
            azi_c = azi - 360;
        else
            azi_c = azi;
        end
        
        itd = round(delay_max * sin(deg2rad(azi_c))); % ITD: diferencia de tiempo entre oídos

        ild = 0.8 * abs(sin(deg2rad(azi_c))); % ILD: diferencia de nivel entre oídos

        base_pos = 16; % Posición base del impulso en la HRIR

        idxL = base_pos - itd;  % Posición del impulso para el oído izquierdo

        idxL = max(1, min(hrtf_length, idxL)); % Asegura que el índice esté dentro del vector

        hL = zeros(hrtf_length,1); % Inicializa la HRIR izquierda

        hL(idxL) = 1 - ild; % Coloca el impulso con atenuación por ILD

        idxR = base_pos + itd; % Posición del impulso para el oído derecho

        idxR = max(1, min(hrtf_length, idxR)); % Asegura que el índice esté dentro del vector

        hR = zeros(hrtf_length,1); % Inicializa la HRIR derecha

        hR(idxR) = 1; % Coloca el impulso del oído derecho

        fdiff = 0.1 * fir1(20, 0.5)'; % Filtro FIR para simular diferencias espectrales

        fdiff = [fdiff; zeros(hrtf_length - numel(fdiff),1)]; % Ajusta la longitud del filtro

        % Aplica diferencias de frecuencia a ambos oídos
        hL = hL + fdiff; 
        hR = hR + fdiff;
        
        hL = hL / (max(abs(hL)) + 1e-12); % Normaliza la HRIR izquierda
        hR = hR / (max(abs(hR)) + 1e-12); % Normaliza la HRIR derecha

        HRIRs(i,1,:) = reshape(hL,1,1,[]); % Guarda HRIR izquierda
        HRIRs(i,2,:) = reshape(hR,1,1,[]);% Guarda HRIR derecha
    end


    % PROCESAMIENTO - Overlap-Add FFT----------------------------------%

    Nfft = 2^nextpow2(bloque_size + hrtf_length - 1); % Tamaño de la FFT potencia de 2

    win = hann(bloque_size, 'periodic'); % Ventana Hann para cada bloque

    % BUFFER DE SALIDA-------------------------------------------------%

    n_out = n_muestras_total + bloque_size + hrtf_length; % Tamaño del buffer con margen extra

    audio3D = zeros(n_out, 2); % Buffer estéreo de salida

    num_frames = ceil((n_muestras_total - bloque_size)/hop) + 1; % Número total de bloques a procesar

    num_frames = max(num_frames, 1); % Garantiza al menos un bloque


    frame_idx = 1; % Índice inicial del bloque

    h_wait = waitbar(0, 'Simulando carros a los lados...'); % Barra de progreso

    % ELEVACIÓN--------------------------------------------------------%

    % Rango de elevación en grados
    elev_min = -30; 
    elev_max = 60;
    
    elev_freq = 1/duracion_simulacion; % Una oscilación completa durante la simulación

    % BUCLE PRINCIPAL-------------------------------------------------%

    for f = 1:num_frames % Recorre todos los bloques de audio

        idx_inicio = frame_idx;
        idx_fin = min(frame_idx + bloque_size - 1, n_muestras_total);

        bloque_audio = audio_source(idx_inicio:idx_fin);
        % Extrae el bloque de audio

        % Rellenar con ceros si el bloque es corto
        if numel(bloque_audio) < bloque_size
            bloque_audio = [bloque_audio; zeros(bloque_size - numel(bloque_audio), 1)];
        end
        

        bloque_win = bloque_audio .* win; % Aplica la ventana

        % Movimiento azimutal de 2 vueltas
        t_center = (idx_inicio + bloque_size/2) / Fs; % Tiempo central del bloque

        progreso = t_center / duracion_simulacion; % Progreso de la simulación (0 a 1)

        azimuth_actual = mod(progreso * 360 * 2, 360); % El sonido da dos vueltas completas


        % Elevación sintética
        %elev_progress = t_center / duracion_simulacion; % Progreso temporal de elevación

        elev_actual = (sin(2*pi*elev_freq*t_center) * 0.5 + 0.5); % Oscilación suave entre 0 y 1

        elev_deg = elev_min + elev_actual * (elev_max - elev_min); % Elevación final en grados

        % Interpolar HRIR por azimut
        h_interp = interpolate_hrtf_circular(HRIRs, angulos_az, azimuth_actual); % HRIR interpolada para el ángulo actual

        % Aplicar efecto de elevación sintética
        h_interp_elev = apply_elevation_synthetic(h_interp, elev_deg, Fs); % Modifica la HRIR según la elevación


        % FFT convolution
        X = fft(bloque_win, Nfft);
        % FFT del bloque de audio

        H_L = fft(h_interp_elev(:,1), Nfft);
        % FFT de la HRIR izquierda

        H_R = fft(h_interp_elev(:,2), Nfft);
        % FFT de la HRIR derecha

        yL = real(ifft(X .* H_L));
        % Convolución izquierda

        yR = real(ifft(X .* H_R));
        % Convolución derecha


        len_conv = bloque_size + hrtf_length - 1;
        % Longitud válida de la convolución

        yL = yL(1:len_conv); 
        yR = yR(1:len_conv);
        % Recorte del resultado


        idx_out1 = idx_inicio;
        idx_out2 = idx_out1 + len_conv - 1;
        % Índices de salida


        % Protección de índices
        buf_len = size(audio3D,1);
        i1 = max(1, idx_out1);
        i2 = min(buf_len, idx_out2);

        if i2 >= i1
            start_in_conv = 1 + (i1 - idx_out1);
            end_in_conv   = start_in_conv + (i2 - i1);
            segL = yL(start_in_conv:end_in_conv);
            segR = yR(start_in_conv:end_in_conv);
            audio3D(i1:i2,1) = audio3D(i1:i2,1) + segL;
            audio3D(i1:i2,2) = audio3D(i1:i2,2) + segR;
        end
        % Overlap-Add seguro


        frame_idx = frame_idx + hop;
        % Avanza al siguiente bloque

        if mod(f,10) == 0 || f==1 || f==num_frames
            waitbar(min(1, idx_inicio/n_muestras_total), h_wait, ...
                sprintf('Az: %.1f°, El: %.1f°', azimuth_actual, elev_deg));
        end
        % Actualiza la barra de progreso
    end

    close(h_wait);
    % Cerrar la barra de progreso

    % POST-PROCESADO-----------------------------------------------%
    final_len = n_muestras_total + hrtf_length; % Longitud útil final

    final_len = min(final_len, size(audio3D,1)); % Asegura tamaño correcto

    audio3D = audio3D(1:final_len, :); % Recorta el audio final


    % normalizar
    audio3D = audio3D / (max(abs(audio3D(:))) + 1e-12) * 0.98;

    % VISUALIZACIÓN Y SALIDA------------------------------------------%

    t_out = (0:size(audio3D,1)-1)/Fs;
    % Vector de tiempo


    subplot(2,1,1);
    plot((0:length(audio_source)-1)/Fs, audio_source);
    title('Original');
    xlim([0 min(2, duracion_simulacion)]);
    % Gráfica del audio original


    subplot(2,1,2);
    plot(t_out, audio3D(:,1), t_out, audio3D(:,2));
    legend('Izq','Der');
    title('Salida 3D');
    xlim([0 min(2, duracion_simulacion)]);
    % Gráfica del audio 3D


    sound(audio3D, Fs);
    % Reproducir audio final
end

% interpolate_hrtf_circular
% Interpola HRIRs entre dos ángulos cercanos de forma circular
function h_out = interpolate_hrtf_circular(HRIRs, angulos, azi_target)

    azi_target = mod(azi_target,360);
    % Asegura que el ángulo esté entre 0 y 360 grados

    n_ang = length(angulos);
    % Número total de ángulos disponibles

    h_len = size(HRIRs,3);
    % Longitud de cada HRIR (número de muestras)

    angles_ext = [angulos, angulos + 360];
    % Duplica los ángulos para manejar el caso circular (ej. 350° → 0°)

    idx_high_ext = find(angles_ext >= azi_target, 1, 'first');
    % Encuentra el primer ángulo mayor o igual al objetivo

    if isempty(idx_high_ext)
        idx_high_ext = 1;
    end
    % Si no encuentra ninguno, vuelve al inicio

    idx_low_ext = idx_high_ext - 1;
    % Índice del ángulo inferior

    if idx_low_ext < 1
        idx_low_ext = idx_low_ext + n_ang;
    end
    % Ajuste circular si se sale del rango

    idx_low = mod(idx_low_ext-1, n_ang) + 1;
    % Índice real del ángulo inferior dentro de HRIRs

    idx_high = mod(idx_high_ext-1, n_ang) + 1;
    % Índice real del ángulo superior dentro de HRIRs

    a_low = angulos(idx_low); 
    % Ángulo inferior real

    a_high = angulos(idx_high);
    % Ángulo superior real

    if a_high < a_low
        a_high = a_high + 360;
    end
    % Corrección circular (ej. de 350° a 360°)

    azi_cmp = azi_target; 
    if azi_cmp < a_low
        azi_cmp = azi_cmp + 360;
    end
    % Ajusta el ángulo objetivo para comparación correcta

    denom = (a_high - a_low);
    % Diferencia angular total

    if denom == 0
        w = 0;
    else
        w = (azi_cmp - a_low) / denom;
    end
    % Peso de interpolación (0 → ángulo bajo, 1 → ángulo alto)

    h1 = squeeze(HRIRs(idx_low, :, :));   
    % HRIR del ángulo inferior [2 x h_len]

    h2 = squeeze(HRIRs(idx_high, :, :));  
    % HRIR del ángulo superior [2 x h_len]

    h_interp = (1-w) * h1 + w * h2;
    % Interpolación lineal entre ambas HRIRs

    h_out = h_interp.';                   
    % Transpone para obtener formato [h_len x 2]

end

% FUNCIÓN: apply_elevation_synthetic
% Aplica un efecto de elevación modificando la HRIR
function h_mod = apply_elevation_synthetic(h, elev_deg, Fs)

    % h: HRIR [h_len x 2]
    % elev_deg: elevación en grados
    % Fs: frecuencia de muestreo

    h_len = size(h,1);
    % Longitud de la HRIR


    % === MAPEO DE ELEVACIÓN A FRECUENCIA ===
    elev_min = -30; 
    elev_max = 60;
    % Rango de elevación permitido

    fc_min = 1000; 
    fc_max = min(8000, Fs/2*0.9);
    % Rango de frecuencia de corte

    frac = (elev_deg - elev_min) / (elev_max - elev_min);
    % Normaliza la elevación a rango 0–1

    frac = min(max(frac, 0), 1);
    % Limita el valor entre 0 y 1

    fc = fc_min + frac * (fc_max - fc_min);
    % Calcula la frecuencia de corte según elevación


    % === FILTRO PASA-BAJOS ===
    ord = 64;
    % Orden del filtro FIR

    Wn = fc / (Fs/2);
    % Frecuencia normalizada

    Wn = min(max(Wn, 0.001), 0.999);
    % Evita valores extremos

    b = fir1(ord, Wn);
    % Diseña el filtro FIR pasa-bajos


    % === SEPARAR GRAVES Y AGUDOS ===
    lowL = filter(b, 1, h(:,1));
    % Parte grave del canal izquierdo

    lowR = filter(b, 1, h(:,2));
    % Parte grave del canal derecho

    highL = h(:,1) - lowL;
    % Parte aguda del canal izquierdo

    highR = h(:,2) - lowR;
    % Parte aguda del canal derecho


    % === GANANCIAS SEGÚN ELEVACIÓN ===
    gain_low = 1.0;
    % Graves se mantienen iguales

    gain_high = 0.7 + 0.9 * frac;
    % Agudos aumentan cuando la fuente está arriba


    hL_mod = gain_low * lowL + gain_high * highL;
    % HRIR izquierda modificada

    hR_mod = gain_low * lowR + gain_high * highR;
    % HRIR derecha modificada


    % === ILD VERTICAL (EFECTO SUTIL) ===
    ild_offset = (-0.06) + (0.12 * frac);
    % Pequeña variación de nivel según elevación

    hL_mod = hL_mod * (1 - ild_offset/2);
    % Ajuste leve del canal izquierdo

    hR_mod = hR_mod * (1 + ild_offset/2);
    % Ajuste leve del canal derecho


    % === NORMALIZACIÓN FINAL ===
    hL_mod = hL_mod / (max(abs(hL_mod)) + 1e-12);
    % Normaliza HRIR izquierda

    hR_mod = hR_mod / (max(abs(hR_mod)) + 1e-12);
    % Normaliza HRIR derecha


    h_mod = [hL_mod, hR_mod];
    % Devuelve la HRIR modificada [h_len x 2]

end


