# Risoluzione numerica di un sistema di equazioni differenziali ordinarie
Lo scopo di questo progetto è la risoluzione numerica di un sistema lineare di ODE (ordinary differential equations) caratterizzato dall'avere un alto coefficiente 
di stiffness ed essere pertanto soggetto a generare una forte instabilità numerica se approcciato con un metodo numerico non consono. 
Il problema è stato risolto con vari algoritmi numerici e con vari tipi di mesh. Come mostrano i risultati, la A-stabilità ed ancor più, la L-stabilità
sono due proprietà fortemente desiderabili in un algoritmo risolutore, infatti, in loro assenza si rende necessario utilizzare un numero esorbitante 
di nodi per approssimare numericamente la soluzione. I metodi numerici che posseggono solo la zero-stabilità si sono rivelati assolutamente non 
in grado di approssimare efficacemente la soluzione, coerentemente con le previsioni teoriche.  
I codici sono stati implementati in linguaggio MATLAB, tuttavia, essendo tutti gli algoritmi stati scritti from scratch (da zero) utilizzando 
unicamente la libreria standard di MATLAB, i codici risultano pienamente compatibili sia con le licenze base di MATLAB 
(senza la necessità di specifici toolbox a pagamento), sia con il software open-source GNU Octave.  
Inoltre, una relazione molto più dettagliata di questo file readme è stata inserita nel repository in formato PDF.  
Sia i codici che la relazione sono stati sviluppati come progetto individuale.

## Il problema

Si consideri il problema di Cauchy

$$
y \hspace{0.1cm}'(t)=M \cdot y(t), \qquad y(t_0)=y_0,
$$

dove M è definita come la matrice di dimensione $15 \times 15$, avente come j-esimo elemento della diagonale principale il valore $-j^2$ e come elementi 
della sovradiagonale i valori $150$:

$$
M :=
\begin{pmatrix}
-1     & 150    & 0      & \cdots    & 0 \\
0      & -4     & 150    & \cdots    & 0 \\
0      & 0      & -9     &           & \vdots \\
\vdots & \vdots &        & \ddots    & 150 \\
0      & 0      & \cdots & 0         &-225
\end{pmatrix}.
$$

e dove

$$
y_0 = \[ 1, 1, 1, ..., 1, 10 \]^T 
$$

da risolvere con $t \in \[ t_0 , T \]$, sia con $T=10$, sia con $T=10000$.  
Per calcolare numericamente la soluzione, l'intervallo temporale $[t_0,T]$ viene suddiviso in una successione di nodi:  

$$
t_0 < t_1 < \cdots < t_N=T ,
$$

chiamata *mesh* e su ogni nodo della mesh viene approssimato il valore della funzione (vettoriale).  
La distanza $h_n=t_{n+1}-t_n$ tra due nodi consecutivi viene detto *passo di discretizzazione*, se il passo di discretizzazione è lo stesso per 
tutti i nodi della mesh, la mesh viene detta *uniforme* oppure *omogenea*.


## Analisi teorica a priori

$M$ ha $15$ autovalori distinti ed ha spettro: $\sigma (M) = \\{ −1,−4,…,−225 \\} $, in particolare, tutti e $15$ i suoi autovalori sono reali negativi e 
pertanto il sistema di ODE in esame descrive una dinamica asintoticamente stabile per il teorema di Lyapunov.  
Tuttavia, la scrittura della soluzione come combinazione lineare di termini della forma $e^{\lambda_i \cdot (t- t_0) } \cdot v_i$ mostra l'esistenza di transitori veloci e transitori lenti, in altre parole, il problema forza l'utilizzo di tantissimi nodi negli istanti iniziali, mentre negli istanti finali sono necessari
molti meno nodi per approssimare efficacemente la soluzione. Conseguenza diretta di ciò è che utilizzare una mesh omogenea potrebbe non essere la scelta più 
vantaggiosa.  
Come noto dall'analisi, i sistemi lineari di equazioni differenziali ordinarie ammettono come soluzione esatta 
$y(t) = expm( M \cdot (t - t_0)) \cdot y_0$, dove:

$$
expm(A) := \sum_{k=0}^{+ \infty} \frac{A^k}{k!}
$$

La funzione *expm* è detta esponenziale di matrice, ed in MATLAB è calcolabile mediante la funzione di libreria standard *expm*.


## Il metodo di Radau IIA











