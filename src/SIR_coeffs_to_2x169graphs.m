%% SIR model for coronavirus outbreak dynamics driven by social stress
% 
% MATLAB R2016b required
% 
% The script allows you to get the figures from Appendix A of the article
% regarding the COVID-19 outbreak spread model SIR_SS and the global
% analysis of epidemics' dynamics in populations affected by social stress.
% 
% Authors: Kastalskiy, I. A., Zinovyev, A. Y., Mirkes, E. M.,
%          Kazantsev, V. B., and Gorban, A. N.
% 
% A model of the COVID-19 pandemic is proposed, combining the dynamics of
% social stress described by the tools of sociophysics with classical
% epidemic models.
% The model parameters have been successfully fitted to best match the
% statistical observations of epidemics in different countries worldwide.
% 
% Please note that the fit result may vary slightly from that in the paper,
% depending on the chosen time step, 'dt' (default is 1 day).
% Additionally, the result may be influenced by potential data updates in
% the file 'owid-covid-data.xlsx', such as changes in countries' populations
% over time and retrospective data corrections.
% 
% Innokentiy Kastalskiy, kastalskiy@neuro.nnov.ru
% Copyright 2023
% 
% 
%% Loading COVID-19 data and parameters of countries of the world
% 
% Preferred use of the data file format, as in the project "Our World in Data": https://ourworldindata.org/coronavirus
% 
% Direct link to Excel file: https://covid.ourworldindata.org/data/owid-covid-data.xlsx
% 
% Data repository of Our World in Data project:
% https://github.com/owid/covid-19-data/tree/master/public/data
% 
% File Table_parameters.xlsx is also required.

if ~exist( 'table_data_par', 'var' )
    clearvars
    [ filename, pathname ] = uigetfile( '*.xls*', 'Select xlsx data file', 'owid-covid-data-10.05.2021.xlsx' );
    filename_par = uigetfile( '*.xls*', 'Select xlsx data file with parameters', 'Table_parameters.xlsx' );
    cd( pathname );
    table_data = readtable( filename );
    table_data_par = readtable( filename_par );
end
 
iso3 = table_data_par.ISO_3;
iso2 = table_data_par.ISO_2;
country = table_data_par.Country;
cont = table_data_par.Continent;
reg = table_data_par.Region;
pop = table_data_par.Population;
Ndays = table_data_par.Ndays;
thr_cases = table_data_par.Threshold;
GDP = table_data_par.GDP_per_capita_PPP;
a_ini = table_data_par.a;
K2_ini = table_data_par.K2*1e-3;
q_ini = table_data_par.q*1e3;
Io_ini = table_data_par.Io*1e-6;
R2 = table_data_par.R2;
R2_result = zeros(length(iso3),3);
 
for ii = 1:length(iso3)
 
raw_data = table_data( strcmp( table_data.iso_code, iso3{ii} ), : );
countryname = raw_data.location{1};
pop = raw_data.population(1);
 
%% Setting the parameters of the SIR_SS model
 
T0 = tic;

 
% -------------------------------------------------------------------------------------------------------
% dt - time step for simulations (you can set 1, 0.5, 0.25, 0.2, 0.125, 0.1, ..., 0.01 or even smaller)
%      reducing the time step improves the accuracy of calculations, but increases their duration

dt = 1;
recdt = round(1/dt);

% -------------------------------------------------------------------------------------------------------
 
 
arrsize = (Ndays(ii)-1)*recdt+1;
TCC = [ 0; raw_data.total_cases ];
ix = find( TCC >= thr_cases(ii), 1 );
data = TCC( ix:ix+Ndays(ii)-1 )/pop;
datanorm = norm( data-mean(data) );
 
 
% -----------------------------------------------------------------------------------------------------
% Main parameters of the model
% q - stress response rate (see model description [Kastalskiy et al., 2021]). For European countries,
%     q usually takes a value in the range from several 10k to several 100k: (20..300)*1000
% a - infection rate, usually takes a value in the range 0.1..0.4
% b - recovery rate, we take the value 0.1 as the reciprocal of the characteristic recovery period
%     and as a reference for a
% K2 - exhaustion rate (see model description [Kastalskiy et al., 2021]). For European countries,
%      K2 usually takes a value in the range (4..8)*10e-3
% K3 - relaxation rate to ignorant mode (slow), we take the value 0.01 as a reference for K2
% Io - the initial fraction of infected people in the population
 
 
a = a_ini(ii);
q = q_ini(ii);
K2 = K2_ini(ii);
Io = Io_ini(ii);
 
 
b  = 0.1;
K3 = 0.01;
 
 
I = zeros(arrsize,1);
Sign = I; Sres = I; Sexh = I; R = I;
 
%% Performing calculations
 
I(1) = Io;
Sign(1) = 1-Io;
for i = 2:arrsize
    v1 = Sign(i-1);
    v2 = Sres(i-1);
    v3 = Sexh(i-1);
    v4 = I(i-1);
    v5 = R(i-1);
    Sign(i) = v1 + dt*( -q*v1*v4^2 - a*v1*v4 + K3*v3 );
    Sres(i) = v2 + dt*( -K2*v2 + q*v1*v4^2 );
    Sexh(i) = v3 + dt*( -K3*v3 - a*v3*v4 + K2*v2 );
       I(i) = v4 + dt*( -b*v4 + a*v1*v4 + a*v3*v4 );
       R(i) = v5 + dt*b*v4;
end
R2_result(ii,1) = 1-(norm(data-I(1:recdt:end)-R(1:recdt:end))/datanorm)^2;
 
 
%% Plotting figures
 
figure(1)
subplot(17,10,ii+1)
colororder({'k','#0072BD'})
yyaxis left
set( gca, 'FontSize',4 )
plot( 1:dt:dt*(arrsize-1)+1, [Sign'; Sres'; Sexh'], 'LineWidth',0.7 )%, 'Color','#7E2F8E' )
axis([ 0.5 round((arrsize-1)*dt)+1.5 0 1 ])
title(countryname, 'FontWeight','normal' , 'FontSize',4 )
 
yyaxis right
plot( 1:dt:dt*(arrsize-1)+1, I, 'LineWidth',0.7 , 'Color','#0072BD' )
hold on
plot( 1:dt:dt*(arrsize-1)+1, R+I, '-', 'LineWidth',0.7 , 'Color',[0.6350 0.0780 0.1840]*1.2 )
ylim([ 0 R(end)+I(end) ])

 
figure(2)
subplot(17,10,ii+1)
colororder([0 0.4470 0.7410; [0.6350 0.0780 0.1840]*1.2])
yyaxis right
set( gca, 'FontSize',4 )
plot( 2-ix:length(TCC)-ix+1, TCC, '-', 'Color',[1 0.25 0.25] , 'LineWidth',0.7 , 'MarkerSize',4 )  
hold on
plot( 1:dt:round(dt*(arrsize-1))+1, (I+R)*pop, '-', 'Color',[0.6350 0.0780 0.1840]*1.2 , 'LineWidth',0.7 )
axis([ 0.5 round((arrsize-1)*dt)+1.5 0 max(TCC(round((arrsize-1)*dt)+ix),(I(end)+R(end))*pop) ])
title(countryname, 'FontWeight','normal' , 'FontSize',4 )
 
yyaxis left
stem( 3-ix:length(TCC)-ix+1, diff(TCC), 'Marker','none' , 'LineWidth',0.2 , 'Color','#4DBEEE' )
hold on
plot( 1:dt:round(dt*(arrsize-1))+1-dt, diff(I+R)*pop/dt, '-', 'LineWidth',0.7 , 'Color','#0072BD' )
ylim([ 0 max(max(diff(TCC(1:min(length(TCC),round((arrsize-1)*dt)+ix+1)))),max(diff(I+R))*pop/dt) ])
 
end
R2_result(:,2) = R2_result(:,1)-R2;
R2_result(:,3) = R2;
 
toc(T0)

figure(1)
subplot(17,10,2)
legend( 'S\_ign', 'S\_res', 'S\_exh', 'I', 'CC', 'Location','Best', 'EdgeColor','None' )
saveas( gcf, 'SIRss model dynamics.fig' )

figure(2)
subplot(17,10,2)
legend( 'DNC', 'CC''', 'TCC', 'CC', 'Location','Best', 'EdgeColor','None' )
saveas( gcf, 'covid-19 fittings.fig' )
save( [ 'SIR_SS_model_results_dt=' num2str(dt) '.mat' ] )
