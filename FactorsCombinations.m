RendFact = readmatrix("Dataset.xlsx", Sheet="Dataset");
%% 1. Definizione vettori
format long
nDays = input('Numero di giorni osservati: '); % y = 261; 6m = 129
ENPACL = RendFact(1:nDays,2);
R_tar = input('Settare R squared target: '); 
% Fattori di rischio (tutti della stessa lunghezza)
Nfact = input('Numero totale di fattori immessi: ');

TITOLI = zeros(nDays, Nfact);

col_par = input('Colonna di partenza: '); % colonna 4
for i = 1:Nfact
    TITOLI(:,i) = RendFact(1:nDays,col_par + i - 1);
end

%% 2. Verifica NaN
verifica = zeros(1, Nfact);
for i = 1:Nfact
    verifica(1,i) = sum(isnan(TITOLI(:,i)));
end

if sum(verifica) > 0
    disp("Ci sono "+sum(verifica)+" valori NaN");
else
    disp("Non ci sono valori NaN, procedere al punto 3.");
end
%% 3. Scelta del tipo di analisi
modo = input("1 per tutte le combinazioni per un dato numero di variabili (k), 2 per modello personalizzato: ");
if modo == 1
    k = input('Inserire il numero di variabili nel modello: ');
    AllCombinations(ENPACL, TITOLI, k, R_tar, Nfact, nDays);

elseif modo == 2
    k = input('Inserire il numero di variabili nel modello personalizzato: ');
    ModelloPersonalizzato(ENPACL, TITOLI, k, R_tar, Nfact, nDays);
else
    error('Scelta non valida');
end

%% 4. Funzioni locali
function AllCombinations(ENPACL, TITOLI, k, R_tar, Nfact, nDays)
    if Nfact ~= size(TITOLI,2)
        error("l'elenco dei Titoli è stato modificato, il numero di titoli e di fattori totali non corrsiponde");
    end
    
    if k < 1
        error("k deve essere almeno 1");

    elseif k > Nfact
        error("k non può superare il numero di fattori disponibili (%d)", Nfact);

    end

    indici = nchoosek(1:Nfact,k); % fa tutte le possibili combinazioni di k fattori
    count_out = 0;
    
    for c = 1:size(indici,1) % per ogni diversa combinazione
        idx = indici(c,:);
        fattori = TITOLI(:, idx); % serie storiche per quella data combinazione
        calcoli = Operazioni(ENPACL, fattori, nDays, k);

        if k == 1 || calcoli.R2_adj > R_tar
            printModello(idx, calcoli, k);
            count_out = 1;
        end
    end

    if count_out == 0
        fprintf("Non ci sono modelli che superino l'R squared target\n");
    end
end

function  ModelloPersonalizzato(ENPACL, TITOLI, k, R_tar, Nfact, nDays)
    if Nfact ~= size(TITOLI,2)
        error("l'elenco dei Titoli è stato modificato, il numero di titoli e di fattori totali non corrsiponde");
    end
    
    if k < 1
        error("k deve essere almeno 1");
    
    elseif k > Nfact
        error("k non può superare il numero di fattori disponibili (%d).", Nfact);
    
    elseif k > 7
        warning("il codice è pensato per un massimo di 7 variabili. k verrà impostato a 7");
        k = 7;
    end

    fattori = zeros(nDays, k);
    num_fattore = zeros(1,k);
    
    for i = 1:k
        switch mod(i,10)
            case 1
                suffix = "st";
            case 2
                suffix = "nd";
            case 3
                suffix = "rd";
            otherwise
                suffix = "th";
        end
        msg = i+suffix +" titolo (inserire la posizione del titolo all'interno della matrice): ";
        Fn = input(msg);
        fattori(:,i) = TITOLI(:,Fn);
        num_fattore(i) = Fn;
    end
    
    calcoli = Operazioni(ENPACL, fattori, nDays, k);
    printModello(num_fattore, calcoli, k);

    if calcoli.R2_adj < R_tar
        fprintf("Il modello richiesto non ha superato l'R squared target\n");
        decision = input("Reimpostare l'R squared target?\nDigitare 1 per Sì oppure 0 per No: ");
        
        if decision == 1
            R_tar = input('Inserire un nuovo valore per R squared target: ');
            fprintf("Nuovo R squared target impostato a %f, dovrai rilanciare la regressione\n", R_tar);
        else
            fprintf("Va bene! Ti conviene cambiare le componenti del modello allora\n");
        end
    end
end

function output = Operazioni(ENPACL, fattori, nDays, k)
    df = nDays - k - 1;

    medie = mean(fattori);
    A = [ENPACL fattori];
    V_tot = cov(A);

    % vettore delle covarianze tra ENPACL e ciascun fattore
    Cov = zeros(k, 1);
    for i = 1:k
        Cov(i) = V_tot(1, i+1);
    end

    % matrice covarianza dei fattori
    V = cov(fattori);
    inversa = inv(V);

    % coefficienti
    Beta = V\Cov;
    Alpha = mean(ENPACL) - medie*Beta;

    % fitted values e residui
    ENPACL_fit = Alpha + fattori*Beta;
    errori = ENPACL - ENPACL_fit;
    S_2 = (errori' * errori)/df;

    % std error Beta
    StdErr_Beta = zeros(k, 1);
    for j = 1:k
        StdErr_Beta(j) = sqrt(S_2/nDays * inversa(j,j));
    end

    % std error Alpha
    StdErr_Alpha = sqrt(S_2/nDays * (1 + medie / V * medie'));

    % t-stat e p-value
    t_Beta = Beta ./ StdErr_Beta;
    t_Alpha = Alpha / StdErr_Alpha;

    p_Beta = 2 * tcdf(-abs(t_Beta), df);
    p_Alpha = 2 * tcdf(-abs(t_Alpha), df);

    % intervalli di confidenza e significatività
    IC_Beta = zeros(k, 2);
    sB = strings(k, 1);
    for j = 1:k
        lower = Beta(j) - 1.96 * StdErr_Beta(j);
        upper = Beta(j) + 1.96 * StdErr_Beta(j);
        IC_Beta(j,1) = lower;
        IC_Beta(j,2) = upper;
        if lower * upper > 0
            sB(j) = "Sì";
        else
            sB(j) = "No";
        end
    end

    lowerA = Alpha - 1.96 * StdErr_Alpha;
    upperA = Alpha + 1.96 * StdErr_Alpha;
    IC_Alpha = [lowerA upperA];
    if lowerA * upperA > 0
        sA = "Sì";
    else
        sA = "No";
    end

    % bontà di adattamento
    R2 = var(ENPACL_fit) / var(ENPACL);
    R2_adj = 1 - (nDays - 1)/df * (1 - R2);
    F_stat = (R2/k)/((1-R2)/df);

    % struttura dell'output
    output.Alpha = Alpha;
    output.Beta = Beta;
    output.StdErr_Alpha = StdErr_Alpha;
    output.StdErr_Beta = StdErr_Beta;
    output.t_Alpha = t_Alpha;
    output.t_Beta = t_Beta;
    output.p_Alpha = p_Alpha;
    output.p_Beta = p_Beta;
    output.IC_Alpha = IC_Alpha;
    output.IC_Beta = IC_Beta;
    output.signif_Alpha = sA;
    output.signif_Beta = sB;
    output.R2 = R2;
    output.R2_adj = R2_adj;
    output.F = F_stat;
end

function printModello(idx, output, k)
   if k == 1
       fprintf("Fattore di rischio considerato: FATTORE %d\n", idx);
   else
        fprintf("Fattori di rischio considerati: ");
        for i = 1:k
            if i < k
                fprintf("FATTORE %d, ", idx(i));
            else
                fprintf("FATTORE %d\n", idx(i));
            end
        end
   end

   stringa1 = "Std. error";
   stringa2 = "t-stat";
   fprintf("\tcoefficienti%18s %15s\t\tp-value\t\tsignificatività\n", stringa1, stringa2);
   if output.p_Alpha < 0.01
       fprintf("(int)%13f%19f%18f%16f***%14s\n", ...
           output.Alpha, output.StdErr_Alpha, output.t_Alpha, output.p_Alpha, output.signif_Alpha);
   elseif output.p_Alpha < 0.05
       fprintf("(int)%13f%19f%18f%16f**%15s\n", ...
           output.Alpha, output.StdErr_Alpha, output.t_Alpha, output.p_Alpha, output.signif_Alpha);
   elseif output.p_Alpha < 0.1
       fprintf("(int)%13f%19f%18f%16f*%16s\n", ...
           output.Alpha, output.StdErr_Alpha, output.t_Alpha, output.p_Alpha, output.signif_Alpha);
   else
       fprintf("(int)%13f%19f%18f%16f%17s\n", ...
           output.Alpha, output.StdErr_Alpha, output.t_Alpha, output.p_Alpha, output.signif_Alpha);
   end

   for i = 1:k
       if output.p_Beta(i) < 0.01
           fprintf("%18f%19f%18f%16f***%14s\n", ...
               output.Beta(i), output.StdErr_Beta(i), output.t_Beta(i), output.p_Beta(i), output.signif_Beta(i));
       elseif output.p_Beta(i) < 0.05
            fprintf("%18f%19f%18f%16f**%15s\n", ...
               output.Beta(i), output.StdErr_Beta(i), output.t_Beta(i), output.p_Beta(i), output.signif_Beta(i));
       elseif output.p_Beta(i) < 0.1
           fprintf("%18f%19f%18f%16f*%16s\n", ...
               output.Beta(i), output.StdErr_Beta(i), output.t_Beta(i), output.p_Beta(i), output.signif_Beta(i));
       else
           fprintf("%18f%19f%18f%16f%17s\n", ...
               output.Beta(i), output.StdErr_Beta(i), output.t_Beta(i), output.p_Beta(i), output.signif_Beta(i));
       end
   end

   fprintf("\n");
   fprintf("R squared: %f, adjusted R squared: %f, F-stat: %f\n", output.R2, output.R2_adj, output.F);
   fprintf("\n=================================================================================================\n");
end
