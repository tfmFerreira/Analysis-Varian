clear all

%% Experiments to analyse
expno=1:6;
foldername='1H_2Dnut_pulsecalib_';

%% Control parameters

%%% graphics
fontsize=12;
linewidth=1;

%%% baseline 
basel=[0.01 .35 .7  .99];
baseline = 'Y';

%%% phase correction                    
phc0 = 20;
phc1 = 0;
pivot = .45;

%%% ppm calibration
ppmcalib = 0.3964;
ppmcalibval = 0.9;

%%% line broadening
lb = 1;
lb2 = 100;

%% Processing

figure(3), clf, 

pwr=[];
t90=[];

for nexp=1:length(expno)

eval(['cd ' foldername num2str(expno(nexp)) '.fid'])

load output_real.txt
load output_imag.txt
load power.txt
load increment.txt
load Bf1.txt

pwr(nexp)=power(1)+20*log10(power(2)/4095);


Areal=output_real;
Aimag=output_imag;
S=Areal(:,2:length(Areal(1,:)))+i*Aimag(:,2:length(Areal(1,:)));
t=Areal(:,1)%+initial;

%figure(1), clf
%plot(t,real(S(:,1)),'-o')


td=length(t);
td2=length(S(1,:));

if exist('si') == 0, 
            powers=1:20;
            powers=2.^powers;
            sipos=find(powers>td/2);
            si=powers(sipos(3));end
            %si = td/2; end
		if exist('si2') == 0, 
        si2 = td2*8-1; %set to odd number 
        end
 swh = 1e6/(t(2)-t(1));
 dw = t(2)-t(1);
 nu = swh/2*linspace(-1, 1, si)';
 ppm = nu/(Bf1*1e6);
        
        
dw2 = increment*1e-6;
t2 = dw2*(0:(td2-1));
swh2 = 1/dw2;
nu2 = swh2/2*linspace(-1, 1, si2);


%% processing

%%%%%%%% First FID of the experiment                             %%%%%%%%%%        
        Stemp = S(:,1);
        
%%%%%%%% Exponential multiplication for good aesthetics          %%%%%%%%%%
        lbfun = exp(-lb*pi*t);
        Stemp = lbfun.*Stemp;
        
%%%%%%%% Gives Stemp the size (si) by adding zeros               %%%%%%%%%% 
        if si>(td/2), Stemp = [Stemp; zeros(si-td/2,1)]; end
        
        Stemp(1) = 0.5*Stemp(1);

%%%%%%%% Time domain to frequency domain, Stemp - Itemp          %%%%%%%%%%         
        Itemp = flipud(fftshift(fft(Stemp,si)));
        %figure(1), clf, plot(abs(Itemp)), return


%%%%%%%% To check if baseline limits are ok comment on the phc0  %%%%%%%%%%
%%%%%%%% or phc1 lines in the First Section                      %%%%%%%%%%
        basl=[];
        for n=1:2:length(basel)-1
        basl = [basl round(basel(n)*si):round(basel(n+1)*si)];
        end
        basl=basl';
        if exist('phc0') == 0
            figure(1), clf, plot((1:si)/si,real(Itemp),basl/si,zeros(size(basl)),'r.'), cd .., return
        end  
        
        if exist('phc0') == 0
            figure(1), clf, plot((1:si)/si,real(Itemp),basl/si,zeros(size(basl)),'r.'), cd .., return
        end

%%%%%%%% 0th order phase correction (without 1st order)          %%%%%%%%%%        
        pivotpoint = round(pivot*si);
        if exist('phc1') == 0
            Itemp = fPhCorr(Itemp',phc0,0,pivotpoint); Itemp = Itemp.';
            figure(1), clf, plot((1:si)/si,real(Itemp),pivotpoint/si,0,'ro'),cd .., return
        end

%%%%%%%% 0th order and 1st order phase correction                %%%%%%%%%%       
        Itemp = fPhCorr(Itemp',phc0,phc1,pivotpoint); Itemp = Itemp.';
        if exist('ppmcalibval') == 0
            figure(1), clf, plot((1:si)/si,real(Itemp),pivotpoint/si,0,'ro'),cd .., return
        end

%%%%%%%% Cleaning the spectrum in the baseline intervals         %%%%%%%%%%
        if baseline == 'Y'
        PolyCoeff = polyfit(basl,Itemp(basl),2);
        Itemp = Itemp - polyval(PolyCoeff,(1:si)');
        else
        end
        %figure(1), clf, plot(1:si,real(Itemp),basl,zeros(size(basl)),'r.'), return
       

        
        ppmcalibpoint = round(ppmcalib*si);
        ppm = ppm - ppm(ppmcalibpoint) + ppmcalibval;

       

        I1 = real(Itemp); I1 = I1./max(I1);

%%%%%%%% Same procedure as above for all the FID's, that is,     %%%%%%%%%%
%%%%%%%% for each FID along the indirect dimension calculates    %%%%%%%%%%
%%%%%%%% the corresponding Itemp.                                %%%%%%%%%%
%%%%%%%% I will be a matrix with the several Itemp's as columns  %%%%%%%%%%
        for m = 1:td2
            Stemp = S(:,m);
            Stemp = lbfun.*Stemp;
            if si>(td/2), Stemp = [Stemp; zeros(si-td/2,1)]; end
            Stemp = flipud(fftshift(fft(Stemp)));
            %Stemp = Stemp .* exp(-i*(grpdly*2*pi*count));
            Stemp = ifft(fftshift(flipud(Stemp)));
            Stemp(1) = 0.5*Stemp(1);
            Itemp = flipud(fftshift(fft(Stemp,si)));
            Itemp = fPhCorr(Itemp',phc0,phc1,pivotpoint); Itemp = Itemp.';
            if baseline == 'Y'
            PolyCoeff = polyfit(basl,Itemp(basl),2);
            Itemp = Itemp - polyval(PolyCoeff,(1:si)');
            end
            %figure(1), clf, plot(1:si,real(Itemp),basl,zeros(size(basl)),'r.'), title(num2str(m)), pause
            I(:,m) = real(Itemp);
        end

        

    Itemp = I(:,1);
    Itemp1=Itemp;
    INEPT=Itemp;
   % figure(1), plot(ppm,I(:,6)/max(I(:,1))/4,'k'), hold on
    figure(1), clf
    plot(ppm,I(:,1)/max(I(:,1))/4,'k'), hold on

    I2 = real(Itemp); I2 = I2/max(I2);
        
        Stemp2 = I;
        Stemp0=I;
        
        FIDbasl = mean(I(:,round(.75*td2):td2),2);
        [FIDbasl,t2array] = ndgrid(FIDbasl,t2);
        lbfun2 = exp(-lb2*pi*t2array);
        %Stemp2 = Stemp2 - FIDbasl;
        Stemp2 = lbfun2.*Stemp2;
        
        
        Stemp2(:,1) = 0.5*Stemp2(:,1);
        I2 = fftshift(fft(Stemp2,si2,2),2);
        I2 = real(I2); I2 = I2./max(max(I2));
        

        sliceposv=0.86;
        
        
figure(3),

        for nn=1:length(sliceposv)
        slicepos=sliceposv(nn);
        pos=find(ppm > slicepos);
        subplot(1,2,1),
        plot(t2*1e6,Stemp0(pos(1),:)/max(Stemp0(pos(1),:))*0.6+nexp-1,'b-o'), hold on
        plot(t2*1e6,lbfun2(1,:)*0.6+nexp-1,'g-'),
        hold on
        %xlim([t2(1) t2(length(t2))]*1000)
        set(gca,'YTick',[],'LineWidth',1.5*linewidth,'Box','off','TickDir','out','Ycolor',[1 1 1],'FontSize',fontsize*.8)
        xlabel('time / {\mu}s')

        subplot(1,2,2),
        plot(nu2/1000,abs(I2(pos(1),:))/max(abs(I2(pos(1),1:round(si2/2))))+nexp-1,'b-'),hold on
        pos2=find(abs(I2(pos(1),1:round(si2/2))) == max(abs(I2(pos(1),1:round(si2/2)))));
        plot(nu2(pos2)/1000,abs(I2(pos(1),pos2))/max(abs(I2(pos(1),1:round(si2/2))))+nexp-1,'bo')
        title(['power = ',num2str(pwr(nexp)) ' dB'])
        set(gca,'YTick',[],'LineWidth',1.5*linewidth,'Box','off','TickDir','out','Ycolor',[1 1 1],'FontSize',fontsize*.8)
        text(nu2(1)/1000,nexp-0.5,[num2str(0.25e6/abs(nu2(pos2))) ' {\mu}s'])
        xlabel('nutation frequency / kHz')
        end
    
t90(nexp)= 0.25/abs(nu2(pos2));

clear I I2 Stemp0
cd ..        
end

