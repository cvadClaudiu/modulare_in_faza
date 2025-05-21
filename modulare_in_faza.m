% MODULATIA IN FAZA
ap=input('ap=     ');
fp=input('fp=    ');
n=input('n=        ');
fn=input('fn=	');
tj=input('tj=     ');

VT = input('Alegeti viteza de transmisie, normal(1), rapid(2), incet(3)');
defj_a = input('Aplicati defazaj intre purtatori? Nu(0), Da(1): ');
if defj_a == 1
    alpha = input('valoarea alpha (in radiani): ');
    delta = pi/2;  
else
    defj_d = input('Aplicati defazaj diferit in blocul de defazaj: Nu(0), Da(1): ');
    if defj_d == 1
        delta = input('valoare delta (in radiani): ');
    else
        delta = pi/2;  
    end
end

dt=[1 1;1 1;0 0;0 1;1 0;1 0;0 1;0 1];


switch VT
    case 1
        %100 bps/20ms
        t1 = 0      :0.0001:0.02;
        t2 = 0.0201 :0.0001:0.03;
        t3 = 0.0301 :0.0001:0.04;
        t4 = 0.0401 :0.0001:0.06;
        t5 = 0.0601 :0.0001:0.08;
    case 2
        %200 bps/5ms
        t1 = 0      :0.0001:0.005;
        t2 = 0.0051 :0.0001:0.010;
        t3 = 0.0101 :0.0001:0.015;
        t4 = 0.0151 :0.0001:0.020;
        t5 = 0.0201 :0.0001:0.025;
    case 3
        %20 bps/50ms
        t1 = 0      :0.0001:0.05;
        t2 = 0.0501 :0.0001:0.10;
        t3 = 0.1001 :0.0001:0.15;
        t4 = 0.1501 :0.0001:0.20;
        t5 = 0.2001 :0.0001:0.25;
    otherwise
        error('Invalid option');
end

t=[t1 t2 t3 t4 t5];

smp1=ap*cos(2*pi*fp*t1+3*pi/4);
smp2=ap*cos(2*pi*fp*t2+7*pi/4);
smp3=ap*cos(2*pi*fp*t3+pi/4);
smp4=ap*cos(2*pi*fp*t4+5*pi/4);
smp5=ap*cos(2*pi*fp*t5+pi/4);
smp=[smp1 smp2 smp3 smp4 smp5];

sp=ap*cos(2*pi*fp*t);

figure(1)
subplot(211)
plot(t,sp)
xlabel('timp [s]')
title('SEMNAL PURTATOR')

subplot(212);
plot(t,smp);
title('SEMNAL MODULAT IN FAZA')	
xlabel('timp [s]')

figure(2)
p=fft(sp);
pp=p.*conj(p);
s=length(pp);
f=10000*(0:s/2-1)/s;
subplot(211)
plot(f(1:60),pp(1:60));
title('FCT. DENSIT. SPECTR.-SEMNAL PURTATOR')
xlabel('FRECVENTA-Hz')

mp=fft(smp);
mpmp=mp.*conj(mp);
s=length(mpmp);
subplot(212)
plot(f(1:60),mpmp(1:60));
title('FCT. DENSIT. SPECTR.-SEMNAL MODULAT')
xlabel('FRECVENTA-Hz')

eta=n*cos(2*pi*fn*t);
sr=smp+eta;
figure(3)
subplot(211)
plot(t,sr)
title('SEMNAL RECEPTIONAT')
xlabel('timp [s]')

r=fft(sr);
rr=r.*conj(r);
s=length(rr);
subplot(212)
plot(f(1:60),rr(1:60));
title('FCT. DENSIT. SPECTR.-SEMNAL RECEPTIONAT')
xlabel('FRECVENTA-Hz')

if defj_a == 1
    spi1 = sr .* cos(2*pi*fp*t + alpha);
    spi2 = sr .* sin(2*pi*fp*t + alpha);
else
    spi1 = sr .* cos(2*pi*fp*t);
    spi2 = sr .* sin(2*pi*fp*t);
end

figure(4)
subplot(211)
plot(t,spi1)
xlabel('timp [s]')
title('SEMNAL SPI1')
subplot(212)
plot(t,spi2)
xlabel('timp [s]')
title('SEMNAL SPI2')

figure(5)
pi1=fft(spi1);
pp1=pi1.*conj(pi1);
s=length(pp1);
subplot(211)
plot(f(1:60),pp1(1:60));
title('FCT. DENSIT. SPECTR.-SEMNAL SPI1')
xlabel('FRECVENTA-Hz')

pi2=fft(spi2);
pp2=pi2.*conj(pi2);
s=length(pp2);
subplot(212)
plot(f(1:60),pp2(1:60));
title('FCT. DENSIT. SPECTR.-SEMNAL SPI2')
xlabel('FRECVENTA-Hz')

num=[1];
den=[tj 1];
a=lsim(num,den,spi1,t);
b=lsim(num,den,spi2,t);
w=logspace(-1,4);
[mag,phase]=bode(num,den,w);
y=20*log10(mag);
figure(6)
subplot(211)
semilogx(w,y)
grid
title('CARACTERISTICA AMPLITUDINE-FRECVENTA FTJ')
ylabel('dB')
subplot(212)
semilogx(w,phase)
ylabel('rad')
grid
title('CARACTERISTICA FAZA-FRECVENTA FTJ')

figure(7)
subplot(211)
plot(t,a)
xlabel('timp [s]')
title('SEMNAL A')
subplot(212)
plot(t,b)
xlabel('timp [s]')
title('SEMNAL B')

sa=sign(a);
sb=sign(b);
figure(8)
subplot(211)
plot(t,sa)
xlabel('timp [s]')
title('SIGNUM(A)')
subplot(212)
plot(t,sb)
xlabel('timp [s]')
title('SIGNUM(B)')

l=length(t);
for i=0:1:7,
if sa(round(i*l/8+2/3*l/8))>0 & sb(round(i*l/8+2/3*l/8))>0
    dr(i+1,1)=0;
    dr(i+1,2)=0;
elseif sa(round(i*l/8+2/3*l/8))<0 & sb(round(i*l/8+2/3*l/8))>0
    dr(i+1,1)=1;
    dr(i+1,2)=0;
elseif sa(round(i*l/8+2/3*l/8))>0 & sb(round(i*l/8+2/3*l/8))<0
    dr(i+1,1)=0;
    dr(i+1,2)=1;
elseif sa(round(i*l/8+2/3*l/8))<0 & sb(round(i*l/8+2/3*l/8))<0
    dr(i+1,1)=1;
    dr(i+1,2)=1;
end
end

dt=dt
dr=dr

if dr==dt
   disp('TRANSMISIE CORECTA')
else
   disp('TRANSMISIE CU ERORI')
end