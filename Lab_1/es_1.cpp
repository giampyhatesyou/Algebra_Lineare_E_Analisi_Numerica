//CIARLO MICHELE S5337477
//GIAMPIETRO ANDREA S5208458
#include <iostream>
#include <cmath>
#include <math.h>
#include <cfloat>

using namespace std;

//funzione ausiliaria ricorsiva per fattoriale
double fact(int n) {
  if (n==0 || n==1)
    return 1;
  else
    return n*fact(n-1);
}

//esercizio 1
void esercizio1(){
  //valori assegnati dalla matricola di Ciardo Michele
  double d0 = 7;
  double d1 = 7;
  double b = (d1+1)*pow(10,20);
  double c = -b;

  cout << "Esercizio 1:";
  cout << "d0 = " << d0 << "\t d1 = " << d1 << "\t b = " << b << "\t c = " << c << endl;
  //calcolo (a+b)+c e a+(b+c) al variare di a
  for (int i=0; i<=6; i++) {
    double a = (d0+1)*pow(10,i);
    double x = (a+b)+c;
    double y = a+(b+c);
    cout << "Con i = " << i << ":\t";
    cout << "a = " << a << "\t (a+b)+c = " << x << "\t a+(b+c) = " << y << endl;
  }

}

//funzioni ausiliarie per esercizio 2
//funzione ausiliaria ricvorsiva per polinomio di Taylor
double Taylor(double x, int n) {
  if (n<0)
    return 0;
  else
    return (pow(x,n)/fact(n)) + Taylor(x,n-1);
  }

//funzione ausiliaria per algoritmo 1
void algoritmo1(double x, int n) {
  double t = Taylor(x,n);
  cout << "\tfN in x=" << x << " di grado " << n << " vale circa " << t << endl;
  cout << "\tErroe assoluto: " << abs(t-exp(x)) << "\t Errore relativo: " << abs((t-exp(x))/exp(x)) << endl << endl;
}

//funzione ausiliaria per algoritmo 2
void algoritmo2(double x, int n) {
  double r = 1.0/Taylor(-x,n);
  cout << "\t1/fN in x=" << -x << " di grado " << n << " vale circa " << r << endl;
  cout << "\tErroe assoluto: " << abs(r-exp(x)) << "\t Errore relativo: " << abs((r-exp(x))/exp(x)) << endl << endl;
}

void esercizio2(){
  cout << "Esercizio 2:\n";
  //algoritmo 1
  cout << "Algoritmo 1:\n";
  cout << "--Con x=0.5, \n\tf vale: " << exp(0.5) << endl;
  algoritmo1(0.5, 3);
  algoritmo1(0.5, 10);
  algoritmo1(0.5, 50);
  algoritmo1(0.5, 100);
  algoritmo1(0.5, 150);

  cout << "--Con x=30, f vale: " << exp(30) << endl;
  algoritmo1(30, 3);
  algoritmo1(30, 10);
  algoritmo1(30, 50);
  algoritmo1(30, 100);
  algoritmo1(30, 150);

  cout << "--Con x=-0.5, f vale: " << exp(-0.5) << endl;
  algoritmo1(-0.5, 3);
  algoritmo1(-0.5, 10);
  algoritmo1(-0.5, 50);
  algoritmo1(-0.5, 100);
  algoritmo1(-0.5, 150);

  cout << "--Con x=-30, f vale: " << exp(-30) << endl;
  algoritmo1(-30, 3);
  algoritmo1(-30, 10);
  algoritmo1(-30, 50);
  algoritmo1(-30, 100);
  algoritmo1(-30, 150);

  //algoritmo 2
  cout << "Algoritmo 2:\n";
  cout << "--Con x=-0.5, f vale: " << exp(-0.5) << endl;
  algoritmo2(-0.5, 3);
  algoritmo2(-0.5, 10);
  algoritmo2(-0.5, 50);
  algoritmo2(-0.5, 100);
  algoritmo2(-0.5, 150);

  cout << "--Con x=-30, f vale: " << exp(-30) << endl;
  algoritmo2(-30, 3);
  algoritmo2(-30, 10);
  algoritmo2(-30, 50);
  algoritmo2(-30, 100);
  algoritmo2(-30, 150);
}

//funzioni ausiliere per esercizio 3
//funzione ausiliaria per calcolo precisione di macchina per float
float precisione_float(float x) {
  float eps = 1;
  while (1+eps != 1) //finchè eps viene percepito nella somma
    eps = eps/2;
  return eps*2;
}

//funzione ausiliaria per calcolo precisione di macchina per double
double precisione_double(double x) {
  double eps = 1;
  while (1+eps != 1) //finchè eps viene percepito nella somma
    eps = eps/2;
  return eps*2;
}

void esercizio3(){
  cout << "Esercizio 3:\n";
  //stampa della precisione di macchina ottenuta tramite funzione ausiliaria
  cout << "Precisione eps in float: " << precisione_float(0.5) << endl;
  cout << "Precisione eps in double: " <<  precisione_double(0.5) << endl << endl;
}


int main() {
  esercizio1();
  esercizio2();
  esercizio3();
  return 0;
}