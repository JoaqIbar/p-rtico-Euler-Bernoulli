function [U,max]=porticoG4
%ELEMENTO DE PORTICO EULER BERNOULLI
%GEO Y PROP DE MATERIALES
E=30000000000;%N/m2
b=[0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2];%m 35 p2 y 28 v1
h=[0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];%m
%por como se escribe la matriz de conectividad de elementos, los primeros 35 elementos de b y h corresponden a P2, los restantes (últimos 28) a V1
I=[];%vector de inercia de cada elementos
A=[];%vector de áreas de cada elemento
for i=1:length(b)
I=[I,b(i)*h(i)^3/12];%m4
A=[A,b(i)*h(i)];%m^2
end
%MATRIZ DE COORD
%X(m) Y(m)
Mcoord= [0 0;
0 2.8;
0 5.4;
0 8;
0 10.6;
0 13.2;
0 15.8;
0 18.4;
5.6 0;
5.6 2.8;
5.6 5.4;
5.6 8;
5.6 10.6;
5.6 13.2;
5.6 15.8;
5.6 18.4;
12 0;
12 2.8;
12 5.4;
12 8;
12 10.6;
12 13.2;
12 15.8;
12 18.4;
19.2 0;
19.2 2.8;
19.2 5.4;
19.2 8;
19.2 10.6;
19.2 13.2;
19.2 15.8;
19.2 18.4;
25.6 0;
25.6 2.8;
25.6 5.4;
25.6 8;
25.6 10.6;
25.6 13.2;
25.6 15.8;
25.6 18.4];
%MATRIZ DE CONECTIVIDAD
Mconect=[1 2; % de aca para abajo pilares
         2 3;
         3 4;
         4 5;
         5 6;
         6 7;
         7 8;
         9 10;
         10 11;
         11 12;
         12 13;
         13 14;
         14 15;
         15 16;
         17 18;
         18 19;
         19 20;
         20 21;
         21 22;
         22 23;
         23 24;
         25 26;
         26 27;
         27 28;
         28 29;
         29 30;
         30 31;
         31 32;
         33 34;
         34 35;
         35 36;
         36 37;
         37 38;
         38 39;
         39 40;
         2 10; %de aca para abajo vigas
         10 18;
         18 26;
         26 34;
         3 11;
         11 19;
         19 27;
         27 35;
         4 12;
         12 20;
         20 28;
         28 36;
         5 13;
         13 21;
         21 29;
         29 37;
         6 14;
         14 22;
         22 30;
         30 38;
         7 15;
         15 23;
         23 31;
         31 39;
         8 16;
         16 24;
         24 32;
         32 40];

%MATRIZ DE RIGIDEZ
KG = zeros(3*40,3*40);%3 grados de libertad por elemento, 40 nodos
for elem = 1:size(Mconect)(1)
%Nodo inicial y final de cada elemento
Ni = Mconect(elem,1);
Nf = Mconect(elem,2);
%
%Coord del elemento y largo del mismo
X_elem = [Mcoord(Ni,1),Mcoord(Nf,1)];
Z_elem = [Mcoord(Ni,2),Mcoord(Nf,2)];
DeltaX=X_elem(2)-X_elem(1);
DeltaZ=Z_elem(2)-Z_elem(1);
l=sqrt((DeltaX^2) + (DeltaZ^2)); %largo del elem
%
%Grados de libertad del elemento
gdlelem = [3*Ni-2 3*Ni-1 3*Ni 3*Nf-2 3*Nf-1 3*Nf];
%
% se calcula la matriz cambio de base
Icb = [DeltaX/l -DeltaZ/l 0;
       DeltaZ/l  DeltaX/l 0;
         0          0      1];
%
%Matriz de giro
jj=zeros(3,3); %Matriz de ceros de 3x3
Q=[Icb jj;
    jj  Icb];
%
%Matriz de rigidez del elemento en coord locales:
KLelem = E*I(elem)/l^3*[A(elem)*l^2/I(elem)   0 0 -A(elem)*l^2/I(elem) 0 0;
0 12 6*l 0 -12 6*l;
0 6*l  4*l^2  0  -6*l  2*l^2;
-A(elem)*l^2/I(elem) 0 0 A(elem)*l^2/I(elem) 0 0;
0  -12  -6*l  0  12  -6*l;
0  6*l  2*l^2  0  -6*l  4*l^2];
%
%Matriz de rigidez del elemento en coord globales:
KGelem=Q*KLelem*Q.';
%Ubicamos la matriz KGelm en el lugar correspodiente en KG
KG(gdlelem, gdlelem) = KG(gdlelem,gdlelem) + KGelem;
%
end

%Condición de borde: se definen los grados de libertad fijos
%Los nodos 1 9 17 25 y 33 empotarados, se ubican estos 3 grados de libertad para cada nodo:
gdlfijos = [1 2 3 9*3-2 9*3-1 9*3 17*3-2 17*3-1 17*3 25*3-2 25*3-1 25*3 33*3-2 33*3-1 33*3 ];
%Los grados de libertad libres son los restantes
gdllibres=[1:3*size(Mcoord)(1)];
gdllibres(gdlfijos) = [];%Se eliminan los fijos
%Matriz reducida;
Kreducida= KG(gdllibres,gdllibres);
q=5500; %Se define una carga lineal cosntante en N/M2
%Se define el vector de fuerzas nodales externas sin reducir:
F=zeros(40*3,1);
for fs=2:7; %los nodos del 2 al 8 son donde se aplica q, Fx,Fy,Mz

F(fs*3-2)=F(fs*3-2)+q*2.6/2;
F(fs*3)=F(fs*3)+q*2.6*2.6/12;
F((fs+1)*3-2)=F((fs+1)*3-2)+q*2.6/2;
F((fs+1)*3)=F((fs+1)*3)-q*2.6*2.6/12;

endfor
F(gdlfijos)=[]; %Se reduce el vector de fuerzas

%Se resuelve el sistema reducido
U=Kreducida\F;
max=[0 0]; %defino vector auxiliar para ubicar desplazamiento máximo según x, y ubicación del nodo
for r=1:size(U)/3; %recorro solo las entradas de desplazamiento según x (esto se puede hacer gracias a que cada nodo de la matriz reducida no le falta ningún grado de libertad (x,y,thetaZ))
if max(1)<U(3*r-2);
    max=[U(3*r-2) 3*r-2];
endif
endfor
max % el primer valor indica desplazamiento máximo, el segundo el grado de libertad del nodo, este nodo se ubica en N= (ans + 2 )/3


