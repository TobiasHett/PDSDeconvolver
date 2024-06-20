%GUI Deconvolution of PDS time traces
function Deconvmaster(src,event)
clc, clear, close all
fig=figure('Name','PDS-Deconvolver GUI','resize','on','NumberTitle','off');

a = uicontrol(fig,'Style','checkbox','String','Force a+b=1','ToolTip','The side condition a+b=1 must apply','Units','Normalized'); % Checkbox.
a.Position = [0,0.95,0.2,0.05]; % Position of pushbutton (xpos, ypos, xwidth, ywidth)
b = uicontrol(fig,'Style','checkbox','String','Save images','ToolTip','Save images as .png files','Units','Normalized'); % Push button to select folder.
b.Position = [0.20,0.95,0.2,0.05]; % Position of pushbutton (xpos, ypos, xwidth, ywidth)
c = uicontrol(fig,'Style','checkbox','String','Save ASCII files','ToolTip','Save data as ASCII files','Units','Normalized'); % Push button to select folder.
c.Position = [0.37,0.95,0.2,0.05]; % Position of pushbutton (xpos, ypos, xwidth, ywidth)
Run = uicontrol(fig,'Style','pushbutton','String','RUN','ToolTip','Run deconvolution','Units','Normalized');
Run.Position = [0.60,0.95,0.1,0.05];
About = uicontrol(fig,'Style','pushbutton','String','About','ToolTip','About...','Units','Normalized');
About.Position = [0.90,0.95,0.1,0.05];

Run.Callback = @deconvolve;
About.Callback = @aboutwindow;

% %Executes when clicking start.
function deconvolve(src,event)
% Get file paths and names.
 [apo_file,apo_path]=uigetfile('*.dat','Select Apo-file');
 [holo_file,holo_path]=uigetfile('*.dat','Select Holo-file');
 [mix_file,mix_path]=uigetfile('*.dat','Select Mix-file');

%Read input files.
 apo_raw=dlmread(strcat(apo_path,apo_file));
 holo_raw=dlmread(strcat(holo_path,holo_file));
 mix_raw=dlmread(strcat(mix_path,mix_file));
 
% Get length of data vector.
 datadim(1)=length(apo_raw);
 datadim(2)=length(holo_raw);
 datadim(3)=length(mix_raw);
 
% Bring all time traces to the same length.
 rownr=(1:min(datadim))';
 apo=apo_raw((1:min(datadim)),2);
 holo=holo_raw((1:min(datadim)),2);
 mix=mix_raw((1:min(datadim)),2);
 
% Determine time step of data.
 dt=apo_raw(2,1)-apo_raw(1,1);
 
% Plot raw data
subplot(2,2,1)
plot((rownr-1)*dt,apo,(rownr-1)*dt,holo,(rownr-1)*dt,mix)
xlabel('t (탎)')
ylabel('Norm. Int. (a.u.)')
legend('Apo','Holo','Mix')

% Do modulation depth scaling as introduced by G. Jeschke: See Appl.
% Magn. Reson. 30, 473-798 (2006)
% All traces will be scaled to the modulation depth of the Apo state
 
scalefac_holo_to_apo=(sum(log(apo).^2))/(sum(log(apo).*log(holo))); % Relation between Apo and Holo
holosc=exp(scalefac_holo_to_apo*log(holo)); %Holo trace scaled to modulation depth of Apo trace

scalefac_mix_to_apo=(sum(log(apo).^2))/(sum(log(apo).*log(mix))); % Relation between Apo and Mix
mixsc=exp(scalefac_mix_to_apo*log(mix)); %Mix trace scaled to modulation depth of Apo trace

aposc=apo; %Apo remains the same since we scale to the modulation depth of the Apo trace.

% Plot scaled data
subplot(2,2,2)
plot((rownr-1)*dt,aposc,(rownr-1)*dt,holosc,(rownr-1)*dt,mixsc)
xlabel('t (탎)')
ylabel('Norm. Int. (a.u.)')
legend('Apo scaled','Holo scaled','Mix scaled')

% % Do the deconvolution via linear combination method. Mix = a*Apo+b*Holo
if a.Value==1
% Apo+Holo=1 has to apply
fprintf('Order of coefficients:\nApo\nHolo\n')
coeff_holo=[(holosc-aposc)]\(mixsc-aposc);
coeff_apo=[(aposc-holosc)]\(mixsc-holosc);
coeff=[coeff_apo ; coeff_holo]

recalcsc=coeff_apo*aposc+coeff_holo*holosc; %Compute the in-silico mixed trace based on the percentages of apo & holo and scaled apo & scaled holo traces.

Korrelationsmatrix_rescaled=corrcoef(mixsc,recalcsc);
Korrelationsmatrix_rescaled=Korrelationsmatrix_rescaled.^2;

%Rescale back the final results to the original modulation depth of the "mix"-trace
fit_backscaled=exp(log(recalcsc)/scalefac_mix_to_apo);
subplot(2,2,3)
plot((rownr-1)*dt,mix,(rownr-1)*dt,fit_backscaled)
xlabel('t (탎)')
ylabel('Norm. Int. (a.u.)')
legend('Mix','recalc')
Korrelationsmatrix_backscaled=corrcoef(mix,fit_backscaled);
Korrelationsmatrix_backscaled=Korrelationsmatrix_backscaled.^2

else
% Apo+Holo=1 does not need to apply
 fprintf('Order of coefficients:\nApo\nHolo\n')
 coeffs=[aposc holosc]\mixsc

 recalcsc=coeffs(1)*aposc+coeffs(2)*holosc; %Compute the in-silico mixed trace based on the percentages of apo & holo and scaled apo & scaled holo traces.
 
 Korrelationsmatrix_rescaled=corrcoef(mixsc,recalcsc);
 Korrelationsmatrix_rescaled=Korrelationsmatrix_rescaled.^2;
 
% %Rescale back the final results to the original modulation depth of the "mix"-trace
 fit_backscaled=exp(log(recalcsc)/scalefac_mix_to_apo);

 subplot(2,2,3)
 plot((rownr-1)*dt,mix,(rownr-1)*dt,fit_backscaled)
 xlabel('t (탎)')
 ylabel('Norm. Int. (a.u.)')
 legend('Mix','recalc')
 Korrelationsmatrix_backscaled=corrcoef(mix,fit_backscaled);
 Korrelationsmatrix_backscaled=Korrelationsmatrix_backscaled.^2
end
 
% Compute residual of fit
subplot(2,2,4)
residual=mix-fit_backscaled;
plot((rownr-1)*dt,residual)
xlabel('t (탎)')
ylabel('Int. (a.u.)')
legend('Residual')
 
if b.Value==1 %Save png-files
plot1=figure('visible','off');
plot((rownr-1)*dt,apo,(rownr-1)*dt,holo,(rownr-1)*dt,mix)
xlabel('t (탎)')
ylabel('Norm. Int. (a.u.)')
legend('Apo','Holo','Mix')
print(plot1,'-dpng','-r0','01_Raw Data.png')

plot2=figure('visible','off');
plot((rownr-1)*dt,aposc,(rownr-1)*dt,holosc,(rownr-1)*dt,mixsc)
xlabel('t (탎)')
ylabel('Norm. Int. (a.u.)')
legend('Apo scaled','Holo scaled','Mix scaled')
print(plot2,'-dpng','-r0','02_Scaled Data.png')

plot3=figure('visible','off');
plot((rownr-1)*dt,mix,(rownr-1)*dt,fit_backscaled)
xlabel('t (탎)')
ylabel('Norm. Int. (a.u.)')
legend('Mix','recalc')
print(plot3,'-dpng','-r0','03_Fitted Data.png')

plot4=figure('visible','off');
plot((rownr-1)*dt,residual)
xlabel('t (탎)')
ylabel('Int. (a.u.)')
legend('Residual')
print(plot4,'-dpng','-r0','04_Residual.png')
end
 
if c.Value==1 %Save ASCII-files
dlmwrite('01_Raw Data.dat',{(rownr-1)*dt,apo,holo,mix},'delimiter','\t');
dlmwrite('02_Scaled Data.dat',{(rownr-1)*dt,aposc,holosc,mixsc},'delimiter','\t');
dlmwrite('03_Fitted Data.dat',{(rownr-1)*dt,mix,fit_backscaled},'delimiter','\t');
dlmwrite('04_Residual.dat',{(rownr-1)*dt,residual},'delimiter','\t');
end
end

function aboutwindow(src,event)
msgbox({'PDS-Deconvolver GUI';'Version 1.0';'June 2020'},'About','Help')
end

end