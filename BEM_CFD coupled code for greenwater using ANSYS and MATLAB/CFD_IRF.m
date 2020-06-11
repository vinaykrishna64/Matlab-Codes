function [] = CFD_IRF(fn,lamdabyL,CFD)
% fn = 0;
% lamdabyL = 1;
% CFD =1;
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%M
%Computational marine hydrodynamics #1

clc

[TP,nmax,del_t] = creat_inputs(fn,lamdabyL);


nirf= 100;
% del_t = 0.01;
% disp(strcat('delt == ',num2str(del_t)))%0.01; % time_step
n = nmax; % #time_steps
omegaF = 2*pi()/TP ;
M = 1000*163.9; %120*12*6*1025 ; 
Iyy = M*(0.25*31.39)^2;
c33 = 1025*9.81*31.39*4.343;%4.343*31.39;
c35=0;
c53=0;
c55 = M*(31.39*31.39/(12*1.737))*9.81;

% c55
% load('inputdata.mat');
fid1=fopen('excitingF.txt');  %'excitingF'
fid2=fopen('addedmass.txt'); %'addedmass.txt'
fid3=fopen('damp.txt');
num1=fscanf(fid1,'%4d',1);
num2=fscanf(fid2,'%4d',1);
num3=fscanf(fid3,'%4d',1);
F33 = zeros(num1,1);
F55 = zeros(num1,1);
time= zeros(num1,1);
for j = 1:num1
    A = fscanf(fid1,['%f' ''],3);  %3
    %B = fscanf(fid1,['%f' ''],3);
    F33(j,1) = A(2,1);
    time(j,1) = A(1,1);
    F55(j,1) = A(3,1) ; 
end
a33 = zeros(num2,1);
b33 = zeros(num3,1);
a35 = zeros(num2,1);
b35 = zeros(num3,1);
a53 = zeros(num2,1);
b53 = zeros(num3,1);
a55 = zeros(num2,1);
b55 = zeros(num3,1);
omega = zeros(num2,1);
v3 = zeros(nmax+1,1);
v5 = zeros(nmax+1,1);
z3 = zeros(nmax+1,1);
z5 = zeros(nmax+1,1);
for j = 1:num2
    C = fscanf(fid2,['%f' ''],5);  %5
    omega(j,1) = C(1,1); 
    a33(j,1) = C(2,1);
    a35(j,1) = C(3,1);
    a53(j,1) = C(4,1);
    a55(j,1) = C(5,1);
end

for j = 1:num3
B = fscanf(fid3,['%f' ''],5);  %5
b33(j,1) = B(2,1); % B(2,1)
b35(j,1) = B(3,1);
b53(j,1)= B(4,1);
b55(j,1) = B(5,1);
end

ndata = num3;
time2  = time(1:nirf);
B33 = IRF(ndata, nirf,time2, b33,del_t, omega);
B35 = IRF(ndata, nirf,time2, b35,del_t, omega);
B53 = IRF(ndata, nirf,time2, b53,del_t, omega);
B55 = IRF(ndata, nirf,time2, b55,del_t, omega);
A33 = GetAddedMass(omega, del_t,b33,a33,ndata);
A35 = GetAddedMass(omega, del_t,b35,a35,ndata);
A53 = GetAddedMass(omega, del_t,b53,a53,ndata);
A55 = GetAddedMass(omega, del_t,b55,a55,ndata);

z3(1) = 0;
z5(1)=0;
v3(1) = 0;
v5(1)=0;
t(1)=0;
temp_trapz=0;

% % figure
% plot(omega,a35)
% title('a35 vs w')
% % figure
% plot(omega,b35)
% title('b35 vs w')
% figure
% plot(time,F33)
% title('F3 vs t')
% figure
% plot(time,F55)
% title('F5 vs t')
% needed for radiation restoring forces

if fn~=0
    % needed for radiation restoring forces
    a33w = interp1(omega,a33,omegaF);  % calculating a33 freq domain for omegaF by interpolation
    a35w = interp1(omega,a35,omegaF);
    a55w = interp1(omega,a55,omegaF);
    a53w = interp1(omega,a53,omegaF);
    Crad33  = radiation_restoration(a33w,B33,omegaF,time2,A33);
    Crad35  = radiation_restoration(a35w,B35,omegaF,time2,A35); % no need to write twice as it doesnot take velocity as input
    Crad55  = radiation_restoration(a55w,B55,omegaF,time2,A55);
    Crad53 = radiation_restoration(a53w,B53,omegaF,time2,A53);
else
    Crad33  = 0;
    Crad35  = 0; % no need to write twice as it doesnot take velocity as input
    Crad55  = 0;
    Crad53 = 0;
end
%%
fidout = 'RESULTS_OUT.txt';
fopen(fidout ,'wt+');
dlmwrite(fidout,'t  z3  v3  z5  v5  P','newline', 'pc')


% xlswrite('results.xlsx',{'t','z3','v3','z5','v5','P'},strcat('fn=',num2str(fn),',lbyL=',num2str(lamdabyL)),'A1')

for i=1:n-1
    t(i) = (i-1)*del_t;
    Damp33 = DampingCoeff(i,nirf,B33,v3, del_t);
    Damp35 = DampingCoeff(i,nirf,B35,v5, del_t);
    Damp53 = DampingCoeff(i,nirf,B53,v3, del_t);
    Damp55 = DampingCoeff(i,nirf,B55,v5, del_t);
    
  

    

 %%
 if CFD ~= 1   
    cfdforce(i) = 0; % comment while running CFD coupling
 else
 !fluent 2ddp -gr -i 20.txt  
%      % fluent is called with f n closes
     filename = 'pp-1-0.001.dat';
      fresult = 'results1.dat';
     while(~exist(filename,'file')) 
  end
      delete('zerok.dat');
        delete('zerok.cas');
       status = system('cfdpost -session posters.cse &');
      while(~exist(fresult,'file'))
      end
      [status1]=system('"taskkill.exe" /F /im cmd.exe &');
      P(i) = csvread(fresult,5,1);
 % %  fclose(fresult);
      disp(P(i));
      cfdforce(i) = P(i);
 end   
 %%    
     %%%%%%%%%%
     
     
    b3 = F33(i)-c33*z3(i)+ A33*omegaF^2*z3(i)- Damp33+ A35*omegaF^2*z5(i) - del_t*Damp35-c35*z5(i)+(cfdforce(i)*4.343)    - Crad33*z3(i) - Crad35*z5(i);
    b5 = F55(i)-c55*z5(i)+ A55*omegaF^2*z5(i)- Damp55+ A53*omegaF^2*z3(i) - del_t*Damp53-c53*z3(i)+(cfdforce(i)*4.343*((31.39/2)-1))    - Crad55*z5(i) - Crad53*z3(i);
    
    v3(i+1) = v3(i) + del_t*b3/(M);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z3(i+1) = z3(i) + del_t*v3(i+1);
    
    v5(i+1) = v5(i) + del_t*b5/(Iyy);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z5(i+1) = z5(i) + del_t*v5(i+1);
    
    
    %%%%%%%%%%%%%%%%%%%
    
    Damp33 = DampingCoeff(i+1,nirf,B33,v3, del_t);
    Damp35 = DampingCoeff(i+1,nirf,B35,v5, del_t);
    Damp53 = DampingCoeff(i+1,nirf,B53,v3, del_t);
    Damp55 = DampingCoeff(i+1,nirf,B55,v5, del_t);
    
    b3 = F33(i+1)-c33*z3(i+1)+ A33*omegaF^2*z3(i+1) - Damp33+ A35*omegaF^2*z5(i) - del_t*Damp35-c35*z5(i)+(cfdforce(i)*4.343)  - Crad33*z3(i+1) - Crad35*z5(i+1);
    b5 = F55(i+1)-c55*z5(i+1)+ A55*omegaF^2*z5(i+1) - Damp55+ A53*omegaF^2*z3(i) - del_t*Damp53-c53*z3(i)+(cfdforce(i)*4.343*((31.39/2)-1))  - Crad55*z5(i+1) - Crad53*z3(i+1);
    v3(i+1) = v3(i) + del_t*b3/(M);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z3(i+1) = z3(i) + del_t*v3(i+1);
    
    v5(i+1) = v5(i) + del_t*b5/(Iyy);%  ((z(i)-z(i-1))/del_t)+del_t*a;
    z5(i+1) = z5(i) + del_t*v5(i+1);
    
   %if(i==99)
      % dummy = 1.0;
   %end
    
    a=z3(i+1);
    %disp(z3(i+1));
    
    b=z5(i+1);
    %disp(z5(i+1));  
% xlswrite('results.xlsx',[t(i),z3(i),v3(i),z5(i),v5(i),cfdforce(i)],strcat('fn=',num2str(fn),',lbyL=',num2str(lamdabyL)),strcat('A',num2str(i+1)))
dlmwrite(fidout,[t(i),z3(i),v3(i),z5(i),v5(i),cfdforce(i)],'delimiter','\t','newline', 'pc','precision',6,'-append')
disp(strcat('PERCENT COMPLETED == ',num2str((i/n)*100)))
%%
    if CFD == 1
    %     %%%%%
    %     
        fid = fopen('motion2.c','w+');
        fed = fopen('motion.c','r+');
        while ~feof(fed)
          tline = fgets(fed);
          fprintf(fid,tline);
        end
        fid = fopen('motion2.c','r+');
        fed = fopen('motion.c','w+');

        ta = ['a = ' num2str(a) ';\r\n'];
        tb = ['b = ' num2str(b) ';\r\n'];

        linenum = 0;
        frewind(fid);
         while ~feof(fid)

          tline = fgets(fid);
          linenum = linenum+1; 
          if(linenum==11)
                  fprintf(fed,ta);
           elseif(linenum==12)
               fprintf(fed,tb);
           else
                   fprintf(fed,tline);
           end
         end
         fclose(fid);
         fclose(fed);
    %     
         w = 'zerok.dat';
         e = 'zerok.cas';
    %    
         zz='pp-1-0.001.dat'; %num2str(i) '.dat'];
         y='pp-1-0.001.cas'; %num2str(i) '.cas'];
         eval(['!rename ' zz ' ' w]);
        eval(['!rename ' y ' ' e]);
           tic;pause(1);toc;
         delete('results1.dat');
         fclose all;
         
    end
end

figure(101)
plot(t,z3(1:length(t)))
title('z3')
figure(102)
plot(t,z5(1:length(t)))
title('z5')


end

