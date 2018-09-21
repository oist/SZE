function[] = leCycle(z1,z2,z3,z4,ez1, ez2, ez3, ez4, beta)   

%% work
    %% expectation value
    expE1 = ez1;
    expE2 = ez2;
    expE3 = ez3; 
    expE4 = ez4;
    
    S1 = log(z1)+beta*expE1;
    S2 = log(z2)+beta*expE2;
    S3 = log(z3)+beta*expE3;
    S4 = log(z4)+beta*expE4;

    Dexp21 = expE2 - expE1; %insertion
    Dexp32 = expE3 - expE2; %measurement
    Dexp43 = expE4 - expE3; %expansion
    Dexp14 = expE1 - expE4; %removal

    %% partition function approach

    wz21 = 1/beta * (log(z2) - log(z1));
    wz32 = 1/beta * (log(z3) - log(z2));
    wz43 = 1/beta * (log(z4) - log(z3));
    wz14 = 1/beta * (log(z1) - log(z4));

    %% some plot

    figure(1);
    xv = [1,2,3,4,5];
    xw = [1.5, 2.5, 3.5, 4.5];
    zeroLine= zeros(1,length(xv));
    plot(xv,[expE1, expE2, expE3, expE4, expE1], '-.or'); hold on;
    plot(xv,[z1, z2, z3, z4, z1], '-.ob'); hold on;
    plot(xv,[S1, S2, S3, S4, S1], '-.og');
    legend('E', 'Z', 'S');
    xlabel('step in cycle');
    ylabel('E, Z, S');
    
    figure(2);
    plot(xw,[wz21, wz32, wz43, wz14],'-.ob'); hold on;
    plot(xw,[Dexp21, Dexp32, Dexp43, Dexp14],'-.og'); hold on;
    plot(xw,[(Dexp21+wz21), (Dexp32+wz32), (Dexp43+wz43), (Dexp14 +wz14)],'-.om'); hold on;
    plot(xw,[(S2-S1), (S3-S2), (S4-S3), (S1-S4)],'-.oc'); hold on;
    plot(xv, zeroLine, '--k'); 
    legend('W', 'Delta E', 'Q', 'DeltaS');
    xlabel('steps in cycle');
    ylabel('W, DeltaE, Q')

    disp('le Cycle')
    disp('total expectation value energy change')
    disp(string(Dexp21 + Dexp32 +Dexp43 + Dexp14))
    disp('total entropy change')
    disp(string((S2-S1) +(S3-S2) +  (S4-S3) +(S1-S4)))
    disp('total work')
    disp(string(wz21+ wz32 + wz43 + wz14))
    disp('total heat')
    disp(string(Dexp21+wz21 + Dexp32+wz32 + Dexp43+wz43 + Dexp14+wz14))
    disp('--------')
    disp('ignore measurement')
    disp('total expectation value energy change')
    disp(string(Dexp21 + Dexp43 +Dexp14))
    disp('total entropy change')
    disp(string((S2-S1) +  (S4-S3) +(S1-S4)))
    disp('total work')
    disp(string(wz21 + wz43 + wz14))
    disp('total heat')
    disp(string(Dexp21+wz21  + Dexp43+wz43 + Dexp14 + wz14))
    
    
