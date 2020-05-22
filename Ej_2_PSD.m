%EJ_2
clear all;


N=1000;
a = [1,0.3544,0.3508,0.1736,0.2401]; % Coeficientes verdaderos
W = normrnd(0,1,N,1);                % Ruido blanco gaussiano
y = filter(1,a,W);                   % filtro que lo hace AR-4, muestras Y
m=4;                                 % Orden 4, AR-m

%2.1-----------------------------------
%Estimador  coef MV ------------------------------------------
[ sigma_MV, a_MV ] = MV( y, m ) ;
a_MV;
sigma_MV;
%Akaike -------------------------------------------------
M=10; %1<m=M<10 uso el maximo y despues veo cual es el m que minimiza
AIC=zeros(M,1);  
for m=1:M
    [ sigma_MV_AIC, a_MV_AIC ] = MV( y, m );  % Estimacion MV de sigma para orden m
    AIC(m) = 2*(m+1)+2*N*(1+log(2*pi*sigma_MV_AIC)); % Calculo de etrica de Akaike AIC(m)
end
[AIC_min,AIC_min_idx]=min(AIC); %el m minimo esta dado por el que tiene menor metrica
%La menor metrica es AIC_min y el orden es min_idx

[sigma_AIC, a_AIC] = MV(y, AIC_min_idx);

%2.2------------------DELTA_AKAIKE--------------
delta_m = AIC - AIC_min;

%%
AIC;
%PAra poder ver los resultados obtenidos.
%%

%2.3-------------------------------------
fft_puntos=5000;
segmento = 250; 
solap = segmento/2;

%PSD REAL
[H_real,w_real] = freqz(1,a,fft_puntos); 
Sxx_real = abs(H_real).^2;

%Welch
[Sxx_welch , w_welch] = pwelch(y , segmento , solap);
Sxx_welch = Sxx_welch * pi;

%los valores de los coeficientes para MV_orden4 y AIC no tienen el a_0, 
%por lo que se lo debo agregar, ademas de invertir el signo
%MV_orden4
[H_4,w_4] = freqz(1,[1;-a_MV],fft_puntos);  
Sxx_orden4 = sigma_AIC*abs(H_4).^2;

%Akaike_AIC
[H_AIC,w_AIC] = freqz(1,[1;-a_AIC],fft_puntos); 
Sxx_AIC = sigma_AIC*abs(H_AIC).^2;

figure()
semilogy(w_real,Sxx_real,'r','DisplayName','PSD real');
hold on
semilogy(w_welch,Sxx_welch,'DisplayName','Estimador de Welch');
semilogy(w_4,Sxx_orden4,'m','DisplayName','Estimador de MV de orden 4');
semilogy(w_AIC,Sxx_AIC,'k','DisplayName','Estimador por Akaike');
grid on
%grid minor
xlim([0 pi]);
ylabel('Densidad espectral [dB]');
xlabel('Frecuencia [rad]');
legend('show','location','SouthEast');
%%
%Paso a db a mano para comparar y veo luego cual voy a utilizar.
Sxx_real_db = 10*log(Sxx_real);
Sxx_welch_db = 10*log(Sxx_welch);
Sxx_orden4_db = 10*log(Sxx_orden4);
Sxx_AIC_db = 10*log(Sxx_AIC);

figure()
plot(w_real,Sxx_real_db,'DisplayName','PSD real');
hold on
plot(w_welch,Sxx_welch_db,'DisplayName','Estimador de Welch');
plot(w_4,Sxx_orden4_db,'DisplayName','Estimador de MV de orden 4');
plot(w_AIC,Sxx_AIC_db,'DisplayName','Estimador por Akaike');
grid on
%grid minor
xlim([0 pi]);
ylabel('Densidad espectral [dB]');
xlabel('Frecuencia [rad]');
legend('show','location','SouthEast');

%%
%2.4--------------------------------
k=2000;
AIC_min_idx_hist=zeros(1,k);

for i=1:k
    % iid, por lo tanto cada vez que lo haga las genero de nuevo.
    W = normrnd(0,1,N,1);      % Ruido blanco gaussiano
    y = filter(1,a,W);         % filtro que lo hace AR-4, muestras Y
    
    %Estimacion AKAIKE muchas veces
    M=10;
    AIC_hist=zeros(M,1);  
    for m=1:M
        [ sigma_MV_AIC_hist, a_MV_AIC_hist ] = MV( y, m );  % Estimacion MV de sigma para orden m
        AIC_hist(m) = 2*(m+1)+2*N*(1+log(2*pi*sigma_MV_AIC_hist)); % Calculo de metrica de Akaike AIC(m)
    end
    
    [AIC_min_hist,AIC_min_idx_hist(i)]=min(AIC_hist); %el m minimo esta dado por el que tiene menor metrica
    %La menor metrica es AIC_min_hist y el orden es AIC_min_idx_hist
    
end

figure();
hold on;
grid on;
%grid minor;
histogram(AIC_min_idx_hist,'Normalization','pdf');
xlabel('Orden de Akaike');
ylabel('Probabilidad');
xlim([0,10.8]);