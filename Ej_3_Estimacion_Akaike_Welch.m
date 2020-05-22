%EJ_3
clear all;

y = dlmread('Ej4.csv','\n');
y=y';
N = length(y); %cant muestras


%Modelo AR por MV
M = 20; %ordenes posibles
AIC=zeros(M,1);  
for m=1:M
    [ sigma_MV_AIC, a_MV_AIC ] = MV( y, m );  % Estimacion MV de sigma para orden m
    AIC(m) = 2*(m+1)+2*N*(1+log(2*pi*sigma_MV_AIC)); % Calculo de metrica de Akaike AIC(m)
end
[AIC_min,AIC_min_idx]=min(AIC); %el m minimo esta dado por el que tiene menor metrica
%La menor metrica es AIC_min y el orden es min_idx

[sigma_AIC, a_AIC] = MV(y, AIC_min_idx); %Obtension de los parametros del modelo
 %Recordar el signo por la definicion de la ecuacion de recurrencia del AR 
 % o sea que son todos el signo inverso para los coef a_i

delta_m = AIC - AIC_min;
%%
%Grafico Akaike vs Welch-----------------------------
fft_puntos=5000;
segmento = 250; %si el segmento es menor, es menos ruidoso el grafico.
solap = segmento/2;

%Welch
[Sxx_welch , w_welch] = pwelch(y , segmento , solap);
Sxx_welch = Sxx_welch * pi;
%Sxx_welch = 20*log(Sxx_welch);

%Akaike_AIC
[H_AIC,w_AIC] = freqz(1,[1;-a_AIC],fft_puntos); 
Sxx_AIC = sigma_AIC*abs(H_AIC).^2; %aca nose cuanto vale la varianza, quiza no sea unitaria al calcularla 
%Sxx_AIC = 20*log(Sxx_AIC);

figure()
%plot(w_welch,Sxx_welch,'DisplayName','Estimador de Welch');
semilogy(w_welch,Sxx_welch,'DisplayName','Estimador de Welch');
hold on
%plot(w_AIC,Sxx_AIC,'DisplayName','Estimador por Akaike');
semilogy(w_AIC,Sxx_AIC,'DisplayName','Estimador por Akaike');
grid on
%grid minor
xlim([0 pi]);
ylabel('Densidad espectral [dB]');
xlabel('Frecuencia [rad]');
legend('show','location','SouthEast');