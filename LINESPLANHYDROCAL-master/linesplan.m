filename='MOULDED OFFSET.xlsx';
data = dataset('xlsfile',filename);
data=dataset2cell(data);
%stations
stn=[0 0.5 1 1.5 2 3 4 5 6 7 8 8.5 9 9.5 10 ]';
col='B':'S';
for i=1:length(col)
    Lloc(i)=xlsread(filename,'sheet1',strcat(col(i),'19'));
end
z=[0 0.5 1 2 3 4 5 6 7 8 9 10 10.45 12 14 16 18 20];
%forname
s1={'0','0_5','1','2','3','4','5','6','7','8','9','10','10_45','12','14','16','18','20'};

%HYDROSTATIC CALCULATIONS
%common columns
SM=[.25 1 .5 1 .75 2 1 2 1 2 .75 1 .5 1 .25]';
LLMS=[-5 -4.5 -4 -3.5 -3 -2 -1 0 1 2 3 3.5 4 4.5 5]';
h=17.4;
j=1;
%row identifier for station
id=[2 4 5 6 7 8 9 9 9 10 11 12 13 14 16];
while j<19
    sheetname=strcat('WaterLine',s1{j});
    for i=2:16
        WL(i-1)=data{id(i-1),j+1}(1);
    end
    if iscolumn(WL)==0
    WL=WL';
    end
    if j==17 || j==18
        for i=1:15
        if isnan(WL(i))==1
            WL(i)=0;
        end
        end
    end
    
   
    if j==1
        %WL0
Z=0;%change
secarea=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
fnvol=SM.*secarea;
fnlongmom=fnvol.*LLMS;
secareamom=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
fnofvermom=SM.*secareamom;
bdthwl=WL;%change
fnofwparea=bdthwl.*SM;
fnofwpmom=fnofwparea.*LLMS;
fnlongMI=fnofwpmom.*LLMS;
bdth3=bdthwl.*bdthwl.*bdthwl;
fntransMI=SM.*bdth3;
tabi=[stn SM secarea fnvol LLMS fnlongmom secareamom fnofvermom bdthwl fnofwparea fnofwpmom fnlongMI bdth3 fntransMI];
tabi=num2cell(tabi);
del=0;
for i=1:15
    if isnan(fnvol(i))==0
     del=del+4/3*h*fnvol(i);
    end
end
delsw=del*1.025;
delext=del*1.033;

cb=del/(Lloc(1)*bdthwl(9)*2*Z);%change
%amid
amid=0;
cm=amid/(bdthwl(9)*Z);
cp=cb/cm;
lmom=0;
for i=1:15
    if isnan(fnlongmom(i))==0
     lmom=lmom+4/3*h^2*fnlongmom(i);
    end
end
lcb=lmom/del;
vmom=0;
for i=1:15
    if isnan(fnofvermom(i))==0
     vmom=vmom+4/3*h^2*fnofvermom(i);
    end
end
vcb=vmom/del;
awp=0;
for i=1:15
    if isnan(fnofwparea(i))==0
     awp=awp+4/3*h*fnofwparea(i);
    end
end
tcp=1.025*awp/100;
cw=awp/(Lloc(1)*bdthwl(9)*2); %change
wpmom=0;
for i=1:15
    if isnan(fnofwpmom(i))==0
     wpmom=wpmom+4/3*h^2*fnofwpmom(i);
    end
end
lcf=wpmom/awp;
x=awp*lcf^2;
il=0;
for i=1:15
    if isnan(fnlongMI(i))==0
     il=il+4/3*h^3*fnlongMI(i);
    end
end
il=il-x;
bml=il/del;
it=0;
for i=1:15
    if isnan(fntransMI(i))==0
     it=it+4/9*h*fntransMI(i);
    end
end
bmt=it/del;
mct=(1.025*il)/(100*174);%check
tab2={del,delsw,delext,cb,cm,cp,lmom,lcb,vmom,vcb,awp,tcp,cw,wpmom,lcf,il,bml,it,bmt,mct};
xlswrite('hydrocal',tab2,sheetname,'A18')

header={'STATION','S.M','1/2SEC AREA','FN OF VOLUME','LONG LEVER FROM MIDSHIP','FN OF LONG. MOMENT','1/2 SEC AREA MOMENT','FN OF VERTICAL MOMENT','1/2 BREADTH OF W.L','FN OF W.P. AREA','FN OF W.P. MOMOMENT','FN OF LONG. MI','1/2 BREADTH^3','FN OF TRANS MI'};
xlswrite('hydrocal',header,sheetname,'A1');
xlswrite('hydrocal',tabi,sheetname,'A2');
header2={'vol','weightS.W','weightEXT','CB','CM','CP','L-MOM','LCB','V-MOM','VCB','AWP','TPCm','CW','W.P.MOM','LCF','IL','BML','IT','BMT','MCT1cm'};
xlswrite('hydrocal',header2,sheetname,'A17')
%%%       dummy         %%%
dum=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      

H=(z(j)-z(j-1));
secarea=(dum+WL)*(H/2)+secarea;
fnvol=SM.*secarea;
fnlongmom=fnvol.*LLMS;
%1/2secmoment
xx=((dum+WL)*(H/2))*((H/2)+z(j-1));
secareamom=secareamom+xx;
fnofvermom=SM.*secareamom;
bdthwl=WL;%change
if j==17 || j==18
        for i=1:15
        if bdthwl(i)==0
            bdthwl(i)=nan;
        end
        end
    end
fnofwparea=bdthwl.*SM;
fnofwpmom=fnofwparea.*LLMS;
fnlongMI=fnofwpmom.*LLMS;
bdth3=bdthwl.*bdthwl.*bdthwl;
fntransMI=SM.*bdth3;
tabi=[stn SM secarea fnvol LLMS fnlongmom secareamom fnofvermom bdthwl fnofwparea fnofwpmom fnlongMI bdth3 fntransMI];
tabi=num2cell(tabi);
Z=z(j);
del=0;
for i=1:15
    if isnan(fnvol(i))==0
     del=del+4/3*h*fnvol(i);
    end
end
delsw=del*1.025;
delext=del*1.033;

cb=del/(Lloc(j)*bdthwl(9)*2*z(j));%change
%amid

amid=secarea(9);
cm=amid/(bdthwl(9)*Z);
cp=cb/cm;
lmom=0;
for i=1:15
    if isnan(fnlongmom(i))==0
     lmom=lmom+4/3*h^2*fnlongmom(i);
    end
end
lcb=lmom/del;
cw=awp/(Lloc(j)*bdthwl(9)*2); %change
vmom=0;
for i=1:15
    if isnan(fnofvermom(i))==0
     vmom=vmom+4/3*h*fnofvermom(i);
    end
end
vcb=vmom/del;
awp=0;
for i=1:15
    if isnan(fnofwparea(i))==0
     awp=awp+4/3*h*fnofwparea(i);
    end
end
tcp=1.025*awp/100;
cw=awp/(184*22.9);
wpmom=0;
for i=1:15
    if isnan(fnofwpmom(i))==0
     wpmom=wpmom+4/3*h^2*fnofwpmom(i);
    end
end
lcf=wpmom/awp;
x=awp*lcf^2;
il=0;
for i=1:15
    if isnan(fnlongMI(i))==0
     il=il+4/3*h^3*fnlongMI(i);
    end
end
il=il-x;
bml=il/del;
it=0;
for i=1:15
    if isnan(fntransMI(i))==0
     it=it+4/9*h*fntransMI(i);
    end
end
bmt=it/del;
mct=(1.025*il)/(100*174);%check******************
tab2={del,delsw,delext,cb,cm,cp,lmom,lcb,vmom,vcb,awp,tcp,cw,wpmom,lcf,il,bml,it,bmt,mct};
xlswrite('hydrocal',tab2,sheetname,'A18')

header={'STATION','S.M','1/2SEC AREA','FN OF VOLUME','LONG LEVER FROM MIDSHIP','FN OF LONG. MOMENT','1/2 SEC AREA MOMENT','FN OF VERTICAL MOMENT','1/2 BREADTH OF W.L','FN OF W.P. AREA','FN OF W.P. MOMOMENT','FN OF LONG. MI','1/2 BREADTH^3','FN OF TRANS MI'};
xlswrite('hydrocal',header,sheetname,'A1');
xlswrite('hydrocal',tabi,sheetname,'A2');
header2={'vol','weightS.W','weightEXT','CB','CM','CP','L-MOM','LCB','V-MOM','VCB','AWP','TPCm','CW','W.P.MOM','LCF','IL','BML','IT','BMT','MCT1cm'};
xlswrite('hydrocal',header2,sheetname,'A17')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dum=WL;%dummy
 for i=1:15
     if isnan(dum(i))==1
         dum(i)=0;
     end
 end
    end
    for i=1:15
    if isnan(secarea(i))==1
        secarea(i)=0;
    end
    end
    for i=1:15
    if isnan(secareamom(i))==1
        secareamom(i)=0;
    end
    end
    
j=j+1;   
end

