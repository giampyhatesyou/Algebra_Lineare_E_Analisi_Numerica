//CIARLO MICHELE S5337477
//GIAMPIETRO ANDREA S5208458

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

//funzioni ausiliarie
void printM(const vector<vector<float>>& A){
    for(int i=0; i<A.size(); i++){
        for(int j=0; j<A[0].size(); j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

//funzione ausiliaria per il calcolo del fattoriale
float fact(int n){
    if(n==0 || n==1)
        return 1;
    else
        return n*fact(n-1);
}

//funzione ausiliaria per la creazione della matrice di Pascal
vector<vector<float> > generaMP(int size  ){
    vector<vector<float> > P(size , vector<float>(size ));
    for(int i=1; i<=size  ; i++)
        for(int j=1; j<=size  ; j++)
            P[i-1][j-1] = fact(i+j-2)/(fact(i-1)*fact(j-1));
    return P;
}

//funzione ausiliaria per la creazione della matrice tridiagonale
vector<vector<float> > generaMT(int size  ){
    vector<vector<float> > T(size , vector<float>(size));
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++){
            if(i==j)
                T[i][j] = 2;
            else if(abs(i-j)==1)
                T[i][j] = -1;
            else
                T[i][j] = 0;
        }
    return T;
}

//funzione ausiliaria per il calcolo della norma infinito
float n_infinito(const vector<vector<float> >& A){
    int norma_inf = 0;
    for(int i=0; i<A.size(); i++){
        int row_sum = 0;
        for(int j=0; j<A.size(); j++)
            row_sum += abs(A[i][j]);
        if(row_sum > norma_inf)
            norma_inf = row_sum;
    }
    return norma_inf;
}

//funzione ausliaria per il calcolo del vettore b
vector<float> vecB(const vector<float>& x, const vector<vector<float> >& A){
    float size = x.size(); //salvo per evitare di richiamare la funzione size() ogni volta
    vector<float> b(size  , 0);
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            b[i] += (A[i][j] * x[j]);
    return b;
}

void swap(vector<float>& a, vector<float>& b){
    vector<float> temp = a;
    a = b;
    b = temp;
}

//funzione ausiliaria per il calcolo di Gauss con pivoting parziale
vector<float> gausspiv(vector<vector<float>> A, vector<float> b){
    int n = A.size();
    vector<float> aux(n);
    for(int i=0; i<n-1; i++){
        //ricerca del massimo
        int max = i;
        for(int j=i+1; j<n; j++)
            if(abs(A[j][i]) > abs(A[max][i]))
                max = j;
        //scambio delle righe
        if(max != i){
            swap(A[i], A[max]);
            swap(b[i], b[max]);
        }
        //controllo che il pivot non sia troppo piccolo
        if(abs(A[i][i]) < 1.19209e-07){ //epsilon
            cout << "Pivot troppo piccolo" << endl;
            exit(1);
        }
        //normalizzazione della riga i
        float pivot = A[i][i];
        for(int j=i; j<n; j++)
            A[i][j] /= pivot;
        b[i] /= pivot;

        //sottrazione della riga i dalle righe successive
        for(int j=i+1; j<n; j++){
            float m = A[j][i];
            for(int k=i; k<n; k++)
                A[j][k] -= m*A[i][k];
            b[j] -= m*b[i];
        }
    }

    //risoluzione del sistema triangolare superiore
    for(int i=n-1; i>=0; i--){
        float sum = 0;
        for(int j=i+1; j<n; j++)
            sum += A[i][j]*aux[j];
        aux[i] = (b[i] - sum)/A[i][i];
    }
    return aux;
}

//funzione ausiliaria per il calcolo di Gauss senza pivoting
vector<float> gaussnopiv(vector<vector<float>> A, vector<float> b){
    int n = A.size();
    vector<float> aux(n);
    for(int i=0; i<n-1; i++){
        //controllo che il pivot non sia troppo piccolo
        if(abs(A[i][i]) < 1.19209e-07){ //epsilon
            cout << "Pivot troppo piccolo" << endl;
            exit(1);
        }
        //normalizzazione della riga i
        float pivot = A[i][i];
        for(int j=i; j<n; j++)
            A[i][j] /= pivot;
        b[i] /= pivot;

        //sottrazione della riga i dalle righe successive
        for(int j=i+1; j<n; j++){
            float m = A[j][i];
            for(int k=i; k<n; k++)
                A[j][k] -= m*A[i][k];
            b[j] -= m*b[i];
        }
    }
    //risoluzione del sistema
    for(int i=n-1; i>=0; i--){
        float sum = 0;
        for(int j=i+1; j<n; j++)
            sum += A[i][j]*aux[j];
        aux[i] = (b[i] - sum)/A[i][i];
    }
    return aux;
}

//funzione ausiliaria per stampare un vettore
void printV(const vector<float>& x){
    cout << "[";
    for(int i=0; i<x.size(); i++){
        cout << x[i];
        if(i != x.size()-1)
            cout << ", ";
    }
    cout << "]" << endl;
    cout << endl;
}

//funzione ausiliaria per il calcolo del vettore b perturbato
vector<float> Bpert(float n, vector<float> b){
    float size = b.size();
    vector<float> delta_b(size);
    for(int i=0; i<size; i++){
        if(i%2 == 0) 
            delta_b[i] = n * 0.01 * -1;
        else
            delta_b[i] = n * 0.01;
    }
    vector<float> b_pert(size);
    for(int i=0; i<size; i++)
        b_pert[i] = delta_b[i] + b[i];
    return b_pert;
}

int main(){
    //esercizio 1
    cout << "Esercizio 1:\n";
    //matrice A
    vector<vector<float>> A = { {3,1,-1,0}, {0,7,-3,0}, {0,-3,9,-2}, {0,0,4,-10} };
    //matrice B (nel testo ci sono due matrici A che potrebbero generare confusione nel codice)
    vector<vector<float>> B = { {2,4,-2,0}, {1,3,0,1}, {3,-1,1,2}, {0,-1,2,1} };
    
    cout << "Matrice A:\n";
    printM(A);
    cout << "Matrice B:\n";
    printM(B);
    //calcolo della norma infinito
    cout << "Norma infinito di A: " << n_infinito(A) << endl;
    cout << "Norma infinito di B: " << n_infinito(B) << endl << endl;
    
    //matrice di Pascal
    vector<vector<float> > P = generaMP(10);
    
    cout << "Matrice di Pascal:\n";
    printM(P);
    cout << "Norma infinito di P: " << n_infinito(P) << endl << endl;
    
    //matrice tridiagonale generata dalla matricola di Andrea Giampietro
    vector<vector<float> > T = generaMT(10*(5+1) + 8);
    
    cout << "Matrice tridiagonale:\n";
    printM(T);
    cout << "Norma infinito di T: " << n_infinito(T) << endl << endl;
    

    //esercizio 2
    cout << "Esercizio 2:\n";
    
    //costruzione vettore x per ogni matrice
    vector<float> x_A(A.size(), 1);
    vector<float> x_B(B.size(), 1);
    vector<float> x_P(P.size(), 1);
    vector<float> x_T(T.size(), 1);
    
    //vettore b per ogni matrice
    vector<float> b_A = vecB(x_A, A);
    vector<float> b_B = vecB(x_B, B);
    vector<float> b_P = vecB(x_P, P);
    vector<float> b_T = vecB(x_T, T);
    
    //applico gauss
    vector<float> x_A_gauss = gausspiv(A, b_A);
    cout << "Soluzione del sistema A:\n";
    printV(x_A_gauss);
    cout << "Soluzione del sistema A senza pivoting parziale:\n";
    vector<float> x_A_gaussnopiv = gaussnopiv(A, b_A);
    printV(x_A_gaussnopiv);

    vector<float> x_B_gauss = gausspiv(B, b_B); //con pivoting parziale
    cout << "Soluzione del sistema B:\n";
    printV(x_B_gauss);
    //applico gauss senza pivoting per confrontare i risultati
    vector<float> x_B_gaussnopiv = gaussnopiv(A, b_A);
    cout << "Soluzione del sistema B senza pivoting parziale:\n";
    printV(x_B_gaussnopiv);

    cout << "Soluzione del sistema P:\n";
    vector<float> x_P_gauss = gausspiv(P, b_P);
    printV(x_P_gauss);
    cout << "Soluzione del sistema P senza pivoting parziale:\n";
    vector<float> x_P_gaussnopiv = gaussnopiv(P, b_P);
    printV(x_P_gaussnopiv);

    cout << "Soluzione del sistema T:\n";
    vector<float> x_T_gauss = gausspiv(T, b_T);
    printV(x_T_gauss);
    cout << "Soluzione del sistema T senza pivoting parziale:\n";
    vector<float> x_T_gaussnopiv = gaussnopiv(T, b_T);
    printV(x_T_gaussnopiv);
    




    //esercizio 3
    cout << "Esercizio 3:\n";


    //creazione dei vettori con b perturbato 
    vector<float> b_A_pert = Bpert(n_infinito(A), b_A);
    vector<float> b_B_pert = Bpert(n_infinito(B), b_B);
    vector<float> b_P_pert = Bpert(n_infinito(P), b_P);
    vector<float> b_T_pert = Bpert(n_infinito(T), b_T);

    //stampa dei vettori b perturbati e non
    cout << "Vettore b perturbato (matrice A):\n";
    printV(b_A_pert);
    cout << "Vettore b (matrice A):\n";
    printV(b_A);

    cout << "Vettore b perturbato (matrice B):\n";
    printV(b_B_pert);
    cout << "Vettore b (matrice B):\n";
    printV(b_B);

    cout << "Vettore b perturbato (matrice P):\n";
    printV(b_P_pert);
    cout << "Vettore b (matrice P):\n";
    printV(b_P);

    cout << "Vettore b perturbato (matrice T):\n";
    printV(b_T_pert);
    cout << "Vettore b (matrice T):\n";
    printV(b_T);

    //applicazione di gauss con b perturbato e stampa dei risultati
    cout << "Soluzione del sistema A con b perturbato:\n";
    vector<float> x_A_pert = gausspiv(A, b_A_pert);
    printV(x_A_pert);

    cout << "Soluzione del sistema B con b perturbato:\n";
    vector<float> x_B_pert = gausspiv(B, b_B_pert);
    printV(x_B_pert);

    cout << "Soluzione del sistema P con b perturbato:\n";
    vector<float> x_P_pert = gausspiv(P, b_P_pert);
    printV(x_P_pert);

    cout << "Soluzione del sistema T con b perturbato:\n";
    vector<float> x_T_pert = gausspiv(T, b_T_pert);
    printV(x_T_pert);
    cout << "Soluzione del sistema T con b perturbato senza pivoting parziale:\n";
    vector<float> x_T_pertnopiv = gaussnopiv(T, b_T_pert);
    printV(x_T_pertnopiv);

    return 0;
}