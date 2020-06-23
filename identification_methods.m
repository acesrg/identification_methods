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


end

function [data_ident, data_validation] = generate_ident_package(input_signal, output_signal, sample_time, ident_proportion, plot_package)
	%#GENERATE_IDENT_PACKAGE arma e imprime el paquete de datos
	%#
	%# SYNOPSIS generate_ident_package(input_signal, output_signal, sample_time, ident_proportion, plot_package)
	%# INPUT input_signal: (simbólico) la señal de entrada
	%# INPUT output_signal (simbólico) la señal de salida
	%# INPUT sample_time (float) tiempo de muestreo
	%# OUTPUT [data_ident(paquete), data_validation(paquete)]

end

function Gzi = discrete_ident_arx(data, Ts, focus_mode, na, nb, nk, residual_analysis)
    %#DISCRETE_IDENT_ARX arma e imprime el paquete de datos
	%#
	%# SYNOPSIS discrete_ident_arx(data, Ts, focus_mode, na, nb, nk, residual_analysis)
	%# INPUT data(package): 
	%# INPUT Ts(float):
	%# INPUT focus_mode(string):
	%# INPUT na(float):
	%# INPUT nb(float):
	%# INPUT nk(float):
	%# INPUT residual_analysis(boolean):
	%# OUTPUT Gzi(tf):

end

function Gzi_mc = discrete_ident_recursive_least_squares(data, Ts, plot_ident)
	%#DISCRETE_IDENT_RECURSIVE_LEAST_SQUARES método de los minimos cuadrados
	%#
	%# SYNOPSIS discrete_ident_recursive_least_squares(data, Ts, plot_ident)
	%# INPUT data(paquete): 
	%# INPUT Ts(float):
	%# INPUT plot_ident(boolean):
	%# OUTPUT Gzi_mc:

end

function analyze_residuals(data, sys_id, sampling_frequency)
    %#ANALYZE_RESIDUALS análisis de residuos
	%#
	%# SYNOPSIS analyze_residuals(data, sys_id, sampling_frequency)
	%# INPUT data(paquete): 
	%# INPUT sys_id( - ):
	%# INPUT sampling_frequency(float):

end

function validate_identifications(data, Gzi, Gzi_mc)
    %#VALIDATE_IDENTIFICATIONS 
	%#
	%# SYNOPSIS validate_identifications(data, Gzi, Gzi_mc)
	%# INPUT data(paquete): 
	%# INPUT Gzi(tf):
	%# INPUT Gzi_mc(tf):
end
