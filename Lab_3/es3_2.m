%Ciarlo Michele S5337477
%Giampietro Andrea S5208458

%%
%Esercizio 1
%Assegnazione valori
d0=8;
d1=5;
n=10*(d1+1)+d0;

A=diag(ones(1,n-1),1)+eye(n);
E=zeros(n,n); %matrice a elementi nulli

%1a) -> Calcolare autovalori di A e B
E(n,1)=pow2(-n); %assegnato come da richiesta
B=A+E; %matrice perturbata
%mostro le matrici di partenza
%disp(A);
%disp(B);

%calcolo degli autovalori
VA=eig(A);
VB=eig(B);
VBA=VA-VB;
%mostro per il confronto
st = "Autovalore di A:";
disp(st);
disp(VA);
st = "Autovalore di B";
disp(st);
disp(VB);
st = "VB-VA";
disp(st);
disp(VBA);

%confronto con norma
err_rel_in=norm(B-A)/norm(A); %errore relativo di input
err_rel_out=norm(VB-VA)/norm(A); %errore relativo di output
st="Errore di input";
disp(st);
disp(err_rel_in);
st="Errore di output";
disp(st);
disp(err_rel_out);

%1b) -> AtA e BtB 
At=transpose(A);
Bt=transpose(B);
Aa=A*At;
Bb=B*Bt;
VAa=eig(Aa);
VBb=eig(Bb);

%mostro i risultati
st = "Autovalore di Aa";
disp(st);
disp(VAa);
st = "Autovalore di B";
disp(st);
disp(VBb);


%norma
err_rel_in2=norm(Bb-Aa)/norm(Aa); %errore relativo di input
err_rel_out2=norm(VBb-VAa)/norm(Aa); %errore relativo di output
st="Errore di input";
disp(st);
disp(err_rel_in2);
st="Errore di output";
disp(st);
disp(err_rel_out2);

%%
%Esercizio 2
%2a)
%matrice di adiacenza ricavata dal grafo
A2= [1 1 1 1 1 1 1 0 0 0 0;
    1 1 0 0 0 0 0 0 0 0 0;
    1 0 1 0 0 0 0 0 0 1 0;
    1 0 0 1 1 0 0 0 0 1 0;
    1 0 0 1 1 1 0 1 0 0 0;
    1 0 0 0 1 1 0 1 0 0 0;
    1 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 1 1 0 1 1 0 0;
    0 0 0 0 0 0 0 1 1 0 0;
    0 0 1 1 0 0 0 0 0 1 1;
    0 0 0 0 0 0 0 0 0 1 1];
disp(A2);

%2b)
aux=sum(A2); %somma delle colonne della matrice
D=diag(aux); %matrice D
G=A2/D;
disp(G);
%autovettori e autovalori di G
[aVetG, aValG]=eig(G);
st="Autovettori e autovalori di G:";
disp(st);
disp(aValG);
disp(aVetG);

%%
%Esercizio 3
A3=[1 -1 2; -2 0 5; 6 -3 6];
aValA3=eig(A3);
disp(A3);
disp(aValA3);

%3a)
v1=[1; 1; 1];
v2=[3; 10; 4];
%applicazione del metodo delle potenza (tramite funzione ausiliaria)
%lv# -> autovalori dominanti
%k_v# -> numero di iterazioni richieste per raggiungere la convergenza
[l_v1,it_v1]=met_potenze(A3,v1);
[l_v2,it_v2]=met_potenze(A3,v2);

disp(l_v1);
disp(it_v1);
disp(l_v2);
disp(it_v2);

%3b)
pot=2; %potenza di due
B3=inv(A3-pot*eye(3));
disp(B3);
%applico il metodo delle potenze inverse
[l_B3, it_B3]=met_potenze(B3,v1);
l_B3=pot+1/l_B3;
disp(l_B3);

pot=9;
B4=inv(A3-pot*eye(3));
disp(B4);
%applico il metodo delle potenze inverse
[l_B4, it_B4]=met_potenze(B4,v1);
l_B4=pot+1/l_B4;
disp(l_B4);



%funzione ausiliaria per calcolare l'autovalore dominante 
% di una matrice A partendo da un vettore iniziale vett
%tramite il metodo delle potenze
function [l,it]=met_potenze(A,vett)
    it=1; %iterazioni necessarie per convergenza
    x_m=vett/norm(vett);
    y_m1=A*x_m;
    z_m1=y_m1/norm(y_m1);
    l_m1=(x_m'*y_m1)/(x_m'*x_m);
    l_m=0;
    while (abs(l_m-l_m1)>eps)  
        x_m=z_m1;
        y_m1=A*x_m;
        z_m1=y_m1/norm(y_m1);
        l_m=l_m1;
        l_m1=(x_m'*y_m1)/(x_m'*x_m);
        it=it+1;
    end
    l=l_m1; %autovalore dominante
end