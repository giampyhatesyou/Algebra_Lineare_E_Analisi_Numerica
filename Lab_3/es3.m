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
st = "Autovalore di A";
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
