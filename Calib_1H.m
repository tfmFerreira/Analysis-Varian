clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Code for processing results from calibration experiment  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Varian Version    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% First Section.                          %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Definition of the parameters to be used.%%%%%%%%%%%%%%%%%%%%%%%%%%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Experiment to analyse                                    %%%%%%%%%
expno = [1];



%%%%%%%% Cicle for the experiments set in expno                  %%%%%%%%%%
for nexp = 1:length(expno)

%%%%%%%% Goes to folder number expno(nexp)                       %%%%%%%%%%
eval(['cd ' num2str(expno(nexp))])


%%%%%%%%                   Procedure:                             %%%%%%%%%

%%%%%%%%   I- Comment on phc0 and ph1 first and define the        %%%%%%%%%
%%%%%%%%      baseline                                            %%%%%%%%%
%%%%%%%%  II- Uncomment phc0 to do the 0th order phase correction %%%%%%%%%      
%%%%%%%% III- Uncomment phc1 after doing II                       %%%%%%%%%
%%%%%%%%      (might be useful to use pivot)                      %%%%%%%%%
%%%%%%%%  IV- Comment/Uncomment ppmcalibval to calibrate ppm axis %%%%%%%%%
%%%%%%%%   V- Set si2 to odd number                               %%%%%%%%%



%%%%%%%% I - Defining the baseline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basel=[0.01 .35 .7  .99];
baseline = 'Y';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% II and III - PHASE CORRECTION                        %%%%%%%%%%%    
% 
       phc0 = 40;
      %  phc1 = 0;

%%%%%%%% Pivot (used for the first order phase correction)        %%%%%%%%%
pivot = .45;

%%%%%%%%%% IV - CALIBRATION OF THE 13C INEPT SPECTRUM %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%      uncomment these two lines bellow %%%%%%%%%%%%%%%%%%%%%%%%%%
ppmcalib = 0.3964;
ppmcalibval = 0.9;

%%%%%%%% Line broadening for the first and second dimension       %%%%%%%%%
lb = 10;
lb2 = 1;
%%%%%%%% Extra Plotting parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% ppm limits for plotting the spectrum                     %%%%%%%%%
ppmmireturnn = 0;
ppmmax = 80;

%%%%%%%% Data points for the first and second dimension           %%%%%%%%%
LoadPara = 'Y';
LoadRaw = 'Y';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Second Section.                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Loading and Defining all Parameters.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
%%%%%%%% Loads the folder data if LoadPara = Y on First Section  %%%%%%%%%% 
	if LoadPara == 'Y'
        
%%%%%%%% Converts unix file acqus into the acqusconv text file   %%%%%%%%%%
%%%%%%%% (Acqus contains parameters for the first dimension)     %%%%%%%%%%
	    res = unixtc('acqus','acqusconv');
        
%%%%%%%% List of parameters loaded:                              %%%%%%%%%%
%%%%%%%%                                                         %%%%%%%%%%
%%%%%%%%        NBL - ??not used here??                          %%%%%%%%%%
%%%%%%%%    PULPROG - name of the pulse program used (R-PDLF)    %%%%%%%%%%
%%%%%%%%         NS - number of scans used                       %%%%%%%%%%
%%%%%%%%         TD - number of data points of FID*              %%%%%%%%%%
%%%%%%%%         TE - temperature (K)                            %%%%%%%%%%
%%%%%%%%         DE - recicle delay                              %%%%%%%%%%
%%%%%%%%     DIGMOD - digital filtering on or off                %%%%%%%%%%
%%%%%%%%       SW_h - spectral width in hertz                    %%%%%%%%%%
%%%%%%%%     GRPDLY - group delay                                %%%%%%%%%%
%%%%%%%%        BF1 - Larmour frequency for the proton (Hz)      %%%%%%%%%%
%%%%%%%%      d(..) - delays                                     %%%%%%%%%%
%%%%%%%%      p(..) - pulses                                     %%%%%%%%%%
%%%%%%%%      l(..) - increments                                 %%%%%%%%%%
%%%%%%%% * FID - Free Induction Decay                            %%%%%%%%%%

	    nbl = jcampdx('acqusconv','NBL');
	    pulprog = jcampdx('acqusconv','PULPROG')
	    ns = jcampdx('acqusconv','NS');
	    td = jcampdx('acqusconv','TD');
        te = jcampdx('acqusconv','TE');
	    de = jcampdx('acqusconv','DE'); 
	    digmod = jcampdx('acqusconv','DIGMOD');
	    swh = jcampdx('acqusconv','SW_h');
        grpdly = jcampdx('acqusconv','GRPDLY');
        bf1 = jcampdx('acqusconv','BF1');
		d = jcampdx('acqusconv','D');
	    delays = [0:31];
	    for m = 1:length(delays)
	        eval(['d' num2str(delays(m)) ' = d(1+' num2str(delays(m)) ');']);
	    end
		p = jcampdx('acqusconv','P');
	    pulses = [0:32];
	    for m = 1:length(pulses)
	        eval(['p' num2str(pulses(m)) ' = p(1+' num2str(pulses(m)) ');']);
	    end
		l = jcampdx('acqusconv','L');
	    pulses = [0:31];
	    for m = 1:length(pulses)
	        eval(['l' num2str(pulses(m)) ' = l(1+' num2str(pulses(m)) ');']);
        end

%%%%%%%% Converts unix file acqu2s into the acqu2sconv text file %%%%%%%%%%
%%%%%%%% (Acqu2s contains parameters for the second dimension)   %%%%%%%%%%        
	    res = unixtc('acqu2s','acqu2sconv');
        
%%%%%%%% TD2 - number data points in the second dimension        %%%%%%%%%%
	    td2 = jcampdx('acqu2sconv','TD'); %set to odd number 


%%%%%%%% td is rounded to a multiple of 256%%%%%%%%%% 
        td1=td;
		td = 256*ceil(td/256);
        
%%%%%%%% in case si or si2 doesn't exist                         %%%%%%%%%%
%%%%%%%% si: td/2 real points + td/2 imaginary points            %%%%%%%%%%
        if exist('si') == 0, 
            powers=1:20;
            powers=2.^powers;
            sipos=find(powers>td/2);
            si=powers(sipos(3));end
            %si = td/2; end
		if exist('si2') == 0, 
        si2 = td2*8-1; %set to odd number 
        end

%%%%%%%%   dw - dwell time, time between consequent data points  %%%%%%%%%%
%%%%%%%%    t - time vector for the acquisition                  %%%%%%%%%%
%%%%%%%%   nu - frequency vector for the FT of the FID           %%%%%%%%%%
%%%%%%%%  ppm - nu converted into ppm values                     %%%%%%%%%%
	    dw = 1/swh/2;
        t = de*1e-6 + 2*dw*(0:(td/2-1))';
        nu = swh/2*linspace(-1, 1, si)';
        ppm = nu/bf1;

%%%%%%%%  dw2 - time between data points in the indirect dim.    %%%%%%%%%%
%%%%%%%%   t2 - time vector for the indirect dim.                %%%%%%%%%%
%%%%%%%% swh2 - spectral width for the indirect dim.             %%%%%%%%%%
%%%%%%%%  nu2 - frequency vector for FT of indirect dim.         %%%%%%%%%%        
        dw2 = 2*4*20*1e-6;
        t2 = dw2*(0:(td2-1));
        swh2 = 1/dw2;
	    nu2 = swh2/2*linspace(-1, 1, si2);

%%%%%%%% Saves all the parameters defined in file ParaDat        %%%%%%%%%%        
		save ParaDat
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Third Section.    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Signal Processing.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%% LoadRaw activates load of parameters                    %%%%%%%%%%
    if LoadRaw == 'Y'
		load ParaDat

%%%%%%%% Reads the FID's from ser file.                          %%%%%%%%%%
%%%%%%%% (ieee-le is the format used to store numbers in ser)    %%%%%%%%%%
        fid = fopen('ser','r','ieee-le');
        
        %return
%%%%%%%% fread splits the imaginary part from the real part of   %%%%%%%%%%
%%%%%%%% all the FID's                                           %%%%%%%%%%
%%%%%%%% S becomes a vector with real and imaginary components   %%%%%%%%%%          
        S = fread(fid,[2, td*td2/2],'long')';
        S = S(:,1) + i*S(:,2);
        %figure(1), clf, plot(real(S)), return
        %figure(1), clf, plot(imag(S)), return
        %figure(1), clf, plot(real(S(1:td/2))), return
        
%%%%%%%% reshape converts S into a matrix with columns as FID's  %%%%%%%%%% 
%%%%%%%% to each point in the indirect dimension                 %%%%%%%%%%           
        S = reshape(S,td/2,td2);
        S(round(td1/2)+1:td/2,:)=0;
        %figure(1), clf, plot(real(S(:,2))), return

%%%%%%%% First FID of the experiment                             %%%%%%%%%%        
        Stemp = S(:,1);
        %figure(1), clf, plot(real(Stemp)), return

%%%%%%%% Exponential multiplication for good aesthetics          %%%%%%%%%%
        lbfun = exp(-lb*pi*t);
        Stemp = lbfun.*Stemp;
        %figure(1), clf, plot(t,real(Stemp),t,lbfun*max(real(Stemp))), return

%%%%%%%% Gives Stemp the size (si) by adding zeros               %%%%%%%%%% 
        if si>(td/2), Stemp = [Stemp; zeros(si-td/2,1)]; end
        
%%%%%%%% Trick for removing artifact left by digital filtering   %%%%%%%%%%        
        Stemp = flipud(fftshift(fft(Stemp)));
        %count = (1:si)'/si - floor(si/2);
        %Stemp = Stemp .* exp(-i*(grpdly*2*pi*count));
%%%%%%%% Gives back the FID without digital filtering artifact   %%%%%%%%%%
        Stemp = ifft(fftshift(flipud(Stemp)));
        %figure(1), clf, plot(real(Stemp)), return
        
        Stemp(1) = 0.5*Stemp(1);%Why?????

%%%%%%%% Time domain to frequency domain, Stemp - Itemp          %%%%%%%%%%         
        Itemp = flipud(fftshift(fft(Stemp,si)));
        %figure(1), clf, plot(real(Itemp)), return

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


        %basl = [round(basl11*si):round(basl12*si) round(basl21*si):round(basl22*si)]';
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
    plot(I(:,1)/max(I(:,1))/4,'k'), hold on

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 2nd dimension Fourier transform                          %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  save Imatrix ppm I t2
  
  
  
  I2 = real(Itemp); I2 = I2/max(I2);
        
        Stemp2 = I;
        Stemp0=I;
        %FIDbasl = mean(I(:,round(.75*td2):td2),2);
        %[FIDbasl,t2array] = ndgrid(FIDbasl,t2);
        %lbfun2 = exp(-lb2*pi*t2array);
        %Stemp2 = Stemp2 - FIDbasl;
        %Stemp2 = lbfun2.*Stemp2;

        figure(2), clf, plot(Stemp2(:,1)), 
        Stemp2(:,1) = 0.5*Stemp2(:,1);
        I = fftshift(fft(Stemp2,si2,2),2);
        I = real(I); I = I./max(max(I));
        

        sliceposv=[3109 3384 3678 4882];
                figure(3), clf, 

        for nnnn=1:length(sliceposv)
        slicepos=sliceposv(nnnn);
        subplot(1,2,1),
        plot(t2,Stemp0(slicepos,:),'-o'), hold on
        %xlim([t2(1) t2(length(t2))]*1000)
        subplot(1,2,2),
        plot(nu2/1000,I(slicepos,:)),hold on
        %xlim([nu2(1) nu2(length(nu2))])
        
        end

        legend('C18','C17','C16','g2')

        save RPDLF_data I Stemp0 Stemp2 nu2 t2 ppm td2 si2 Itemp INEPT te
  
  
  
  return
  
        
end
end
cd ..
