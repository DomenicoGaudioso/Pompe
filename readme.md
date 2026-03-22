# Condotte – Grado di riempimento (Manning)

Web app professionale per la **determinazione del grado di riempimento** di una
**condotta circolare in moto uniforme** mediante la formula di Manning.

## Stato attuale dell'implementazione

**v2 – completata** (ultima iterazione 2026-03-22)

Funzionalità implementate:
- Risoluzione numerica y(Q) con metodo della bisezione nel ramo monotono (β ≤ 0.938).
- Calcolo di tutti i parametri idraulici: y, β, A, P, Rh, V.
- Portata massima a sezione libera (β ≈ 0.938) e portata a piena sezione (β = 1.0).
- Curva idraulica completa Q(β) e V(β) con 80 punti.
- Grafici Plotly interattivi con marcatura del punto di progetto.
- Validazione dati con messaggi esplicativi.
- Commenti progettuali automatici (β critico, V autocircolante, Manning elevato).
- Export CSV della curva idraulica.
- Layout `wide`, 5 metriche sintetiche, 3 tab (Risultati / Grafici / Note tecniche).

## Struttura del progetto

```text
Condotte/
├── app.py           # UI Streamlit (nessuna formula)
├── src.py           # Geometria, Manning, bisezione, curva, validazione, commenti
├── requirements.txt # streamlit, numpy, pandas, plotly
├── readme.md        # questo file
└── prompt.txt       # prompt per iterazioni future
```

## Input

| Parametro | Descrizione | Unità |
|-----------|-------------|-------|
| D         | Diametro interno | m |
| n         | Coefficiente di Manning | – |
| S         | Pendenza idraulica | – |
| Q         | Portata di progetto | m³/s |

## Output principali

- Tirante `y [m]`, grado di riempimento `β = y/D`, velocità `V [m/s]`
- Area bagnata `A [m²]`, perimetro bagnato `P [m]`, raggio idraulico `Rh [m]`
- Portata massima `Q_max [m³/s]` (β ≈ 0.938) e `Q_full [m³/s]` (β = 1.0)
- Curve Plotly `Q(β)` e `V(β)`
- CSV della curva idraulica completa

## Avvio rapido

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Nota metodologica

Per la sezione circolare in moto uniforme, la portata **non è monotona** in β:
raggiunge un massimo a β ≈ 0.938, poi diminuisce fino a β = 1.0 (piena sezione).
La risoluzione numerica avviene nel ramo monotono (β ∈ [0, 0.938]).
Se la portata di progetto supera Q_max, la condotta è sottodimensionata.

## Estensioni future consigliate

- Confronto interattivo fra più diametri o pendenze
- Sezioni non circolari (ovoidale, a uovo, rettangolare)
- Modalità inversa: dato β, calcolare Q o D
- Report PDF
- Preset materiali tipici (PVC, cls, grès)
