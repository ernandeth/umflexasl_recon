% Script to simulate and optimize the timing of background suppression 
% pulses.  
% Assumes that the sequence is like this:
% presaturation - label_duration - interval1-  BGS1 - interval2 - BGS2 - interval3 - readout
% Post Label Delay = interval1 + interval2 + interval3

clear
close all

% script options:
doVSI = 1;  % this means that the label also inverts all the other species
do7t = 0;
do055t = 0;
showEvolution = 0;
doSPI = 0;

dt = 1e-3;    % time steps for the simulations (in seconds)

% Assumed Physical Constants - relaxation, efficiency , etc.
% Inversion efficiency of BGS pulses
alpha = 0.98;

%  R1 values for 3T
R1_a = 1/1.67;    % arterial blood from Lu paper

% For Brain
R1_1 = 1/1.3;
R1_2 = 1/0.9;
R1_3 = 1/2.4;

% For Placenta Imaging
% from https://ichgcp.net/clinical-trials-registry/publications/274711-quantitative-t1-and-t2-mapping-by-magnetic-resonance-fingerprinting-mrf-of-the-placenta-before-and
R1_1 = 1/1.8;    % placenta T1
R1_2 = 1/0.85;   % Liver
R1_3 = 1/2.5;    % amniotic fluid from https://pmc.ncbi.nlm.nih.gov/articles/PMC8511125/pdf/nihms-1739882.pdf 


% sequence timing parameters
% time between presat and VS label... OR ... PCASL duration
tag_length = 2.5;  
tag_time = tag_length;
% post label delay
PLD = 1.3;

% intervals before each BGS inversion pulse in seconds
interval1=[0.1:0.01:PLD*0.9];
interval2=[0.1:0.01:PLD*0.9];

% time for readout is the label plus the PLD (ms)
AQtime = ceil((tag_length + PLD) / dt);
AQdur = 500;

% arterial magnetization - control case
Mz_art_con= ones(AQtime + AQdur,1);
% arterial magnetization - tag case
Mz_art_tag = ones(AQtime+ AQdur,1);

% other magnetization species :  th stuff to be suppressed:
Mz_1 = ones(AQtime+ AQdur,1);
Mz_2 = ones(AQtime+ AQdur,1);
Mz_3 = ones(AQtime+ AQdur,1);

% intialise labeling:
Mz_art_con(1)= 1;
Mz_art_tag(1)= 1;

% Allocalte space  the blood signal observed after subtraction
asl = ones(length(interval1), length(interval1));
% Allocate space for the results:
bkgnd_1 = asl;
bkgnd_2 = asl;
bkgnd_3 = asl;

% Loop over combinations of intervals (BGS_time1 and BGS_time2)
for n1 = 1:length(interval1)    
    for n2=1:length(interval2)
        
        % Event Time markers for the different events.
        a = tag_length;             % end of labeling time
        b = a + interval1(n1);      % BGS1 time
        c = b + interval2(n2);      % BGS2 time

        % Reset magnetization saturate before labeling (Presat puls)
        Mz_1(1) = 0;
        Mz_2(1) = 0;
        Mz_3(1) = 0;

        % Make sure the BGS pulses fit inside the PLD
        if c < tag_length + PLD - 0.01
            
            % Bloch equation T1 decay ...  Euler approx.
            for n=2:AQtime + AQdur
                dMz_2 =  (1 - Mz_2(n-1))*R1_1;
                Mz_2(n) = Mz_2(n-1) + dMz_2*dt;
                
                dMz_3 =  (1 - Mz_3(n-1))*R1_2;
                Mz_3(n) = Mz_3(n-1) + dMz_3*dt;
                
                dMz_1 =  (1 - Mz_1(n-1))*R1_3;
                Mz_1(n) = Mz_1(n-1) + dMz_1*dt;
                
                dMz_art_con =  (1 - Mz_art_con(n-1))*R1_a;
                Mz_art_con(n) = Mz_art_con(n-1) + dMz_art_con*dt;
                
                dMz_art_tag =  (1 - Mz_art_tag(n-1))*R1_a;
                Mz_art_tag(n) = Mz_art_tag(n-1) + dMz_art_tag*dt;
                
                % (a) End of the tagging: Make sure the arterial spins are inverted
                if (abs((n*dt) - tag_time)) < (dt/2)
                    
                    % in the BIR8 case 
                    Mz_art_con(n) = 1;
                    Mz_art_tag(n) = 0;

                    % in VSI case, the stationary spins are inverted
                    if doVSI
                        Mz_1(n) = -alpha*Mz_1(n-1);
                        Mz_2(n) = -alpha*Mz_2(n-1);
                        Mz_3(n) = -alpha*Mz_3(n-1);

                        Mz_art_con(n) = 1;
                        Mz_art_tag(n) = -alpha;
                    end
                end

                % (b) BGS first inversion pulse: applied to the whole brain and the artery
                if abs(n*dt - b)< (dt/2)
                    Mz_1(n) = -alpha*Mz_1(n);
                    Mz_2(n) = -alpha*Mz_2(n);
                    Mz_3(n) = -alpha*Mz_3(n);

                    Mz_art_con(n) = -alpha*Mz_art_con(n);
                    Mz_art_tag(n) = -alpha*Mz_art_tag(n);
                end

                % (c) BGS second inversion pulse: applied to the whole brain and the artery
                if abs(n*dt - c)< (dt/2)
                    Mz_1(n) = -alpha*Mz_1(n);
                    Mz_2(n) = -alpha*Mz_2(n);
                    Mz_3(n) = -alpha*Mz_3(n);

                    Mz_art_con(n) = -alpha*Mz_art_con(n);
                    Mz_art_tag(n) = -alpha*Mz_art_tag(n);
                end

            end

            if showEvolution
                figure(1)
                plot(Mz_2,'k'); grid on
                hold on;
                plot(Mz_3,'y');
                plot(Mz_1,'g');
                plot(Mz_art_con,'b')
                plot(Mz_art_tag,'r')
                plot(Mz_art_tag-Mz_art_con,'c')
                legend('Gray','White','CSF','blood','tagged blood','difference','Location','NorthWest')
                rectangle('Position', [AQtime -1 AQdur  2] )
                xlabel('time'); ylabel('M_z');
                grid
                drawnow; %pause(0.2);
                hold off
            end
            
            % update the results for each combination of BGS times
            % depends on the readout type:
            if (doSPI)
                % for spiral projection: average over the readout
                asl(n1,n2) =  mean(Mz_art_tag(AQtime + AQdur)-Mz_art_con(AQtime+ AQdur));
                bkgnd_1(n1,n2) =  mean(Mz_1(AQtime: AQtime + AQdur));
                bkgnd_2(n1,n2) =  mean(Mz_2(AQtime: AQtime + AQdur));
                bkgnd_3(n1,n2) =  mean(Mz_3(AQtime: AQtime + AQdur));
            else
                % for Stack of spirals: center of k-space is at the beginning
                asl(n1,n2) =  Mz_art_tag(AQtime) - Mz_art_con(AQtime);
                bkgnd_1(n1,n2) =  Mz_1(AQtime);
                bkgnd_2(n1,n2) =  Mz_2(AQtime);
                bkgnd_3(n1,n2) =  Mz_3(AQtime);
            end
            
        end
    end

    if showEvolution
        figure (2)
        plot(interval2, abs(asl(n1,:))); hold on;
        plot(interval2, abs(bkgnd_1(n1,:)),'k');
        plot(interval2, abs(bkgnd_2(n1,:)),'y');
        plot(interval2, abs(bkgnd_3(n1,:)),'g');hold off
        
        xlabel('interval 2'); ylabel('M_z');
        legend('ASL signal','Grey', 'White','CSF','Location','NorthWest')
        title(sprintf('interval 1 = %0.3f',interval1(n1)));
    end
end
%%
% weights for placenta and brain
TotalBGsignal = 0.25*abs(bkgnd_1) + 0.15*abs(bkgnd_2) + 0.65*abs(bkgnd_3);

% Liver
%TotalBGsignal = 0.15*abs(bkgnd_1) + 0.8*abs(bkgnd_2) + 0.05*abs(bkgnd_3);

% identify combinations with small amounts of background (<10%)
smallSignal_inds = find(...
    abs(bkgnd_1(:)) < 0.1 & ...
    abs(bkgnd_2(:)) < 0.1 & ...
    abs(bkgnd_3(:)) < 0.1 ...
    ) ; 

% identify the combination with the least background signal
[best bestind] = min(TotalBGsignal(:));

[bn1, bn2] = ind2sub( size(bkgnd_1), bestind);

[sbn1, sbn2] = ind2sub( size(bkgnd_1), smallSignal_inds);


for n=1:length(sbn1)
    fprintf('\nLess than 15 percent signal from any of them  when: \n\tBS1 time:  %f sec., \n\tBS2 time: %f sec.', ...
        interval1(sbn1(n)), interval2(sbn2(n)) );
    fprintf('\n\tSignals for BG 1: %f \t BG 2: %f \t BG 3: %f', bkgnd_1(smallSignal_inds(n)), bkgnd_2(smallSignal_inds(n)) , bkgnd_3(smallSignal_inds(n)))
end

fprintf('\nLEAST total *WEIGHTED* signal (%f) when: \n\tBS1 time:  %f sec., \n\tBS2 time: %f sec.', ...
    TotalBGsignal(bestind), interval1(bn1), interval2(bn2) );
fprintf('\nT1 values :\t T1_1: %f,  \t T1_2: %f  \t T1_3: %f sec.', ...
    1/R1_1, 1/R1_2, 1/R1_3);
fprintf('\nSignals for BG 1: %f \t BG 2: %f \t BG 3: %f', bkgnd_1(bestind), bkgnd_2(bestind) , bkgnd_3(bestind));



% Show the combinations of BGS times for the different background species
figure(3)

subplot(221)
imagesc(TotalBGsignal); colormap gray
axis xy
title('Total BG signal')
xlabel('interval 2'); ylabel('interval 1');
hold on
plot(sbn2,sbn1, '*r')
plot(bn2,bn1, '*g')
colormap jet
colorbar
caxis([-1 1])

subplot(221)
imagesc(TotalBGsignal); colormap gray
axis xy
title('Total Weighted BG signal')
xlabel('interval 2'); ylabel('interval 1');
hold on
plot(sbn2,sbn1, '*r')
plot(bn2,bn1, '*g')
colormap jet
colorbar
caxis([-1 1])

subplot(222)
imagesc(bkgnd_1); colormap gray
axis xy
title('BG 1 signal')
xlabel('interval 2'); ylabel('interval 1');
hold on
plot(sbn2,sbn1, '*r')
plot(bn2,bn1, '*g')
colormap jet
colorbar
caxis([-1 1])

subplot(223)
imagesc(bkgnd_2); colormap gray
axis xy
title(' BG 2 signal')
xlabel('interval 2'); ylabel('interval 1');
hold on
plot(sbn2,sbn1, '*r')
plot(bn2,bn1, '*g')
colormap jet
colorbar
caxis([-1 1])


subplot(224)
imagesc(bkgnd_3); colormap gray
axis xy
title(' BG 3 signal')
xlabel('interval 2'); ylabel('interval 1');
hold on
plot(sbn2,sbn1, '*r')
plot(bn2,bn1, '*g')
colormap jet
colorbar
caxis([-1 1])

