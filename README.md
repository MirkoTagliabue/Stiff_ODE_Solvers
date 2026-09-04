# Risoluzione numerica di un sistema di equazioni differenziali ordinarie
Lo scopo di questo progetto è la risoluzione numerica di un sistema lineare di ODE (ordinary differential equations) caratterizzato dall'avere un alto coefficiente 
di stiffness ed essere pertanto soggetto a generare una forte instabilità numerica se approcciato con un metodo numerico non consono. 
Il problema è stato risolto con vari algoritmi numerici e con vari tipi di mesh. Come mostrano i risultati, la A-stabilità ed ancor più, la L-stabilità
sono due proprietà fortemente desiderabili in un algoritmo risolutore, infatti, in loro assenza si rende necessario utilizzare un numero esorbitante 
di nodi per approssimare numericamente la soluzione. I metodi numerici che godono solo della zero-stabilità si sono rivelati assolutamente non 
in grado di approssimare efficacemente la soluzione, coerentemente con le previsioni teoriche.  
I codici sono stati implementati in linguaggio MATLAB, tuttavia, essendo tutti gli algoritmi stati scritti from scratch (da zero) utilizzando 
unicamente la libreria standard di MATLAB, i codici risultano pienamente compatibili sia con le licenze base di MATLAB 
(senza la necessità di specifici toolbox a pagamento), sia con il software open-source GNU Octave.  
Inoltre, una [relazione](./Analisi_di_stabilità.pdf) molto più dettagliata di questo file README è stata inserita nel repository in formato PDF.  
Sia i codici che la relazione sono stati sviluppati come progetto individuale.

## Il problema

Si consideri il problema di Cauchy

$$
y \hspace{0.1cm}'(t)=M \cdot y(t), \qquad y(t_0)=y_0,
$$

dove M è definita come la matrice di dimensione $15 \times 15$, avente come $j$-esimo elemento della diagonale principale il valore $-j^2$ e come elementi 
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
\underline{y_0} = \[ 1, 1, 1, ..., 1, 10 \]^T 
$$

da risolvere per $t \in \[ t_0 , T \]$,  con $T$ fissato e variabile nell'intervallo $\[ 10 , 10000 \]$.  
Per calcolare numericamente la soluzione, l'intervallo temporale $[t_0,T]$ viene suddiviso in una successione di nodi:  

$$
t_0 < t_1 < \cdots < t_N=T ,
$$

chiamata *mesh* e su ogni nodo della mesh viene approssimato il valore della funzione (vettoriale).  
La distanza $h_n=t_{n+1}-t_n$ tra due nodi consecutivi viene detta *passo di discretizzazione*, se il passo di discretizzazione è lo stesso per 
tutti i nodi della mesh, la mesh viene detta *uniforme* oppure *omogenea*.


## Analisi teorica a priori

$M$ ha $15$ autovalori distinti ed ha spettro: $\sigma (M) = \\{ −1,−4,…,−225 \\} $, in particolare, tutti e $15$ i suoi autovalori sono reali negativi e 
pertanto il sistema di ODE in esame descrive una dinamica asintoticamente stabile per il teorema di Lyapunov.  
Tuttavia, la scrittura della soluzione come combinazione lineare di termini della forma $e^{\lambda_i \cdot (t- t_0) } \cdot v_i$ mostra l'esistenza di transitori veloci e transitori lenti, in altre parole, il problema forza l'utilizzo di tantissimi nodi negli istanti iniziali, mentre negli istanti finali sono necessari
molti meno nodi per approssimare efficacemente la soluzione. Conseguenza diretta di ciò è che utilizzare una mesh omogenea potrebbe non essere la scelta più 
vantaggiosa.  
Come noto dall'analisi, i sistemi lineari di equazioni differenziali ordinarie ammettono come soluzione esatta 
$y(t) = exp( M \cdot (t - t_0)) \cdot y_0$, dove:

$$
exp(A) := \sum_{k=0}^{+ \infty} \frac{A^k}{k!}
$$

La funzione *exp(A)* è detta esponenziale di matrice della matrice $A$, ed in MATLAB è calcolabile mediante la funzione di libreria standard *expm*.  

Il coefficiente di stiffness del sistema vale:

$$
C_{stiff} = \frac{ \max\limits_{i=1,\ldots,15} \left| Re( \lambda_i ) \right| }{ \min\limits_{i=1,\ldots,15} \left| Re( \lambda_i ) \right| } = 225 \gg 1
$$

Dal momento che il sistema ha un coefficiente di stiffness così elevato, lo si può catalogare come problema stiff.

La soluzione esatta è la seguente:  
[INSERIRE FOTO SOLUZIONE ESATTA]  
possiamo osservare che essa presenta un picco iniziale nei primissimi valori di 𝑡, dove raggiunge anche valori di $10^7$ nella prima componente 
e $10^5$ nella seconda componente, tuttavia dopo il picco iniziale essa decade molto velocemente, anche solo per $t=30$ il vettore soluzione ha norma euclidea dell’ordine di $10^{-5}$, mentre per $t=10,000$ tutte le componenti della soluzione esatta hanno valore indistinguibile da zero in aritmetica finita.  
Dal momento che la funzione soluzione presenta valori così elevati, risulta significativo anche lo studio dell'errore relativo, oltre a quello dell'errore
assoluto.


## Il metodo di Radau IIA
Il metodo di [Radau_IIA](./Metodi_Numerici/Radau_IIA_5.m) è un metodo numerico di Runge-Kutta basato sul metodo di collocazione Radau IIA a tre stadi, 
si tratta di un metodo di ordine $5$ ed L-stabile. Il presente metodo è stato implementato su mesh omogenea.  
La L-stabilità permette di smorzare efficacemente le componenti associate agli autovalori grandi in modulo. Il prezzo da pagare è che il metodo è anche costoso da un punto di vista computazionale: per ogni nodo della mesh il metodo richiede di risolvere un sistema $45 \times 45$, tuttavia, in questo particolare 
caso dove il sistema di ODE è lineare e la mesh è omogenea, non solo i sistemi da risolvere numericamente sono lineari, ma la matrice $A$ che descrive 
il sistema $Ax_n = b_n$ è sempre la stessa per tutto il metodo, cambia solo il termine noto $b_n$, ciò consente di utilizzare la fattorizzazione $LU$
una volta per tutte all’inizio ed abbattere così il costo computazionale richiesto per risolvere i sistemi ad un costo quadratico di operazioni macchina (*flops*).

Seguono alcuni grafici per rappresentare la soluzione approssimata. I tre grafici seguenti mostrano la prima componente 
$y_1(t)$ del vettore soluzione $\underline{y}(t)$ e sono stati tutti e tre prodotti utilizzando come parametri $t_0 = 0$ e $T = 50$ , 
inoltre nella prima immagine $h = 1$, nella seconda $h = 0.1$, nella terza $h = 0.01$. È decisamente soddisfacente osservare come con anche solo $h = 0.1$
e dunque $N = 501$ nodi (un numero quasi sempre irrisorio quando si parla di cardinalità di una mesh) il metodo produca un errore relativo di $10^{-4}$, 
risultando quindi decisamente affidabile.

[INSERIRE GRAFICI]

Segue inoltre la tabella degli errori commessi:

[INSERIRE TABELLA]


## Il metodo di Gauss-Legendre su mesh non omogenea




## Il metodo di Gauss-Legendre adattivo





## Altri metodi implementati

Il repository comprende inoltre i seguenti metodi numerici:
- [BDF 2](./Metodi_Numerici/BDF_2.m) e [BDF 3](./Metodi_Numerici/BDF_3.m) a passo costante;  
- [Crank-Nicolson su mesh non omogenea](./Metodi_Numerici/Crank_Nicolson_mesh_non_omogenea.m);  
- [Eulero implicito](./Metodi_Numerici/Eulero_Implicito.m);  
- [Pareschi-Russo](./Metodi_Numerici/Pareschi_Russo.m);  
- [Runge-Kutta](./Metodi_Numerici/Runge_Kutta_4.m) classico di ordine 4;  
- [Runge-Kutta-Fehlberg 4(5)](./Metodi_Numerici/Runge_Kutta_Fehlberg_4_5.m) adattivo.  

La trattazione, l'analisi ed i risultati relativi a questi metodi vengono discussi nella [relazione completa](./Analisi_di_stabilità.pdf).


## Esecuzione

1. Aprire la cartella del repository in MATLAB o GNU Octave.
2. Aprire [main.m](./main.m) e impostare la variabile `metodo`, per esempio:

   ```matlab
   metodo = 'Radau_IIA_5';
   % metodo = 'Gauss_Legendre_6_mesh_non_omogenea';
   % metodo = 'Gauss_Legendre_2_6';
   ```

3. Modificare, se lo si desidera, i parametri `t_0`, `T`, il passo `h` uniforme, i tre passi della mesh non omogenea `passo_1`, `passo_2` ed `passo_3`,
    oppure la tolleranza `toll` del metodo adattivo.
4. Eseguire `main.m`.

Lo script genererà i dati iniziali, calcolerà la soluzione numerica e quella esatta, stamperà gli errori e produrrà i grafici delle componenti $y_1$, $y_2$ e $y_{15}$.


## Limiti dell'implementazione

I solver numerici sono stati implementati ed ottimizzati per questo specifico esperimento: assumono un sistema autonomo lineare $y  \hspace{0.1cm}'=My$ di
dimensione $15$ e non costituiscono una libreria ODE per un generale sistema di equazioni differenziali.  
Questa scelta permette di mettere in evidenza le proprietà numeriche dei metodi e di sfruttare direttamente la struttura lineare del problema per ottimizzare
il costo computazionale richiesto dai metodi numerici.





