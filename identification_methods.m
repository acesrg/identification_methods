% Ejercicio 1
clear; close all; clc

%% Apartado A
s = tf('s');
% Planta
Gs = (20*s+500)/((3*s+1)*(5*s+1))
% Periodo de muestreo
Ts = 0.1;
% Discretizacion
Gz = c2d(Gs, Ts, 'zoh')

%% Apartado B y C
% Create an input_prbs_signal
samples = 1200;
band = [0 0.05];
range = [-1, 1];
noise_amplitude = 0.1;
input_prbs_signal = idinput(samples, 'PRBS', band, range);

% Simulate output and add noise
simulated_output = sim(input_prbs_signal, idpoly(Gz));
simulated_noised_output = add_white_noise_to_func(simulated_output, noise_amplitude);

% Armado del paquete de identificacion
ident_proportion = 0.5;  % 50 percent for identification
plot_package = true;
[data_ident, data_test] = generate_ident_package(input_prbs_signal, simulated_noised_output, Ts, ident_proportion, plot_package);

%% Identifico con la herramienta de estimacion de matlab
na = 2; nb = 2; nk = 1;
focus_mode = 'prediction';
residual_analysis = true;
Gzi = discrete_ident_arx(data_ident, Ts, focus_mode, na, nb, nk, residual_analysis);

%% Identifico por minimos cuadrados en forma recursiva
plot_ident = true;
Gzi_mc = discrete_ident_recursive_least_squares(data_ident, Ts, plot_ident)
%% Validacion de resultados
validate_identifications(data_test, Gzi, Gzi_mc)


function noisy = add_white_noise_to_func(clean_signal, noise_amplitude)
	%#ADD_WHITE_NOISE_TO_FUNC agrega ruido blanco a una señal
	%#
	%# SYNOPSIS add_white_noise_to_func(clean_signal, noise_amplitude)
	%# INPUT clean_signal: (simbólico) la señal de entrada
	%# INPUT noise_amplitude: (float) amplitud de la señal de ruido
	%# OUTPUT noisy (simbólico) señal con ruido agregado
	%#

        noiseSignal = awgn (clean_signal, snr)

end

function [data_ident, data_validation] = generate_ident_package(input_signal, output_signal, sample_time, ident_proportion, plot_package)
    % 	GENERATE_IDENT_PACKAGE arma e imprime el paquete de datos
    % 	
    % 	 SYNOPSIS generate_ident_package(input_signal, output_signal, sample_time, ident_proportion, plot_package)
    % 	 INPUT input_signal: (double-sym) la señal de entrada
    % 	 INPUT output_signal (double-sym) la señal de salida
    % 	 INPUT sample_time (double) tiempo de muestreo
    % 	 OUTPUT [data_ident(iddata), data_validation(iddata)]
    % Armo el paquete de datos
    N1 = floor(N/2);
    data = iddata(y, u, Ts);
    data_ident = data(1:N1);
    data_test = data(N1+1:N);
end

function Gzi = discrete_ident_arx(data, Ts, focus_mode, na, nb, nk, residual_analysis)
    %#DISCRETE_IDENT_ARX identifica por ARX según un paquete de datos dado y constantes
	%#definidas, residual_analysis indica si se debe plotear un análisis de residuos
    %# 
	%# SYNOPSIS discrete_ident_arx(data, Ts, focus_mode, na, nb, nk, residual_analysis)
	%# INPUT data(iddata): paquete de identificación
	%# INPUT Ts(double): tiempo de muestreo
	%# INPUT focus_mode(*char): modo para Opt.Focus
	%# INPUT na(double): 
	%# INPUT nb(double):
	%# INPUT nk(double):
	%# INPUT residual_analysis(logical): opción de analizar los residuos
	%# OUTPUT Gzi(tf): función de transferencia del sistema identificado por arx

end

function Gzi_mc = discrete_ident_recursive_least_squares(data, Ts, plot_ident)
	%#DISCRETE_IDENT_RECURSIVE_LEAST_SQUARES método de los minimos cuadrados
	%#
	%# SYNOPSIS discrete_ident_recursive_least_squares(data, Ts, plot_ident)
	%# INPUT data(iddata): paquete de identificación
	%# INPUT Ts(double): tiempo de muestreo
	%# INPUT plot_ident(logical): opción de plotear la identificación
	%# OUTPUT Gzi_mc(tf): función de transferencia identificada
    
    N1=plot_ident;
     %prueba
    n = 3;
    u = data.InputData;
    y = data.OutputData;
    Theta = zeros(4, N1);
    P = 1e12*eye(4);

    for k = n:N1-1
        Phi = [-y(k-1) -y(k-2) u(k-1) u(k-2)];
        K = P*Phi'/(1+(Phi*P*Phi')); 
        Theta(:,k+1) = Theta(:,k)+K*(y(k)-Phi*Theta(:,k));
        P = P-(K*Phi*P);
    end
    Gzi_mc = tf(Theta(3:4,N1)', [1 Theta(1:2,N1)'], Ts)
end

function analyze_residuals(data, sys_id, sampling_frequency)
    %#ANALYZE_RESIDUALS análisis de residuos y plot
	%#
	%# SYNOPSIS analyze_residuals(data, sys_id, sampling_frequency)
	%# INPUT data(iddata): paquete de identificación
	%# INPUT sys_id(idpoly): polinomio de identificación arx
	%# INPUT sampling_frequency(double): frecuencia de muestreo

end


function validate_identifications(data, Gzi, Gzi_mc)
    %#VALIDATE_IDENTIFICATIONS 
	%#
	%# SYNOPSIS validate_identifications(data, Gzi, Gzi_mc)
	%# INPUT data(iddata): paquete de testeo
	%# INPUT Gzi(tf): función de transferencia identificada por arx
	%# INPUT Gzi_mc(tf): función de transferencia identificada por mínimos
  %#                   cuadrados.
end
