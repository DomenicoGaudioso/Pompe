# Pompe – Prevalenza (TDH) e potenza

Web app professionale per il calcolo della **prevalenza totale (TDH)** e della
**potenza richiesta** da una pompa idraulica, con analisi delle perdite distribuite
e concentrate e curva del sistema TDH(Q).

## Stato attuale dell'implementazione

**v2 – completata** (ultima iterazione 2026-03-22)

Funzionalità implementate:
- Fattore di attrito con Swamee–Jain (turbolento) e f = 64/Re (laminare, Re < 2300).
- Calcolo completo di TDH: ΔZ + hf_asp + hf_mand + hK_asp + hK_mand + ΔV²/2g + ΔP.
- Dettaglio numerico di tutti i contributi (Re, f, Darcy-Weisbach, K).
- Breakdown perdite in DataFrame per visualizzazione e export.
- Curva sistema TDH(Q) da 0.10·Q a 2.0·Q (50 punti).
- Validazione dati con messaggi esplicativi.
- Commenti progettuali automatici (regime laminare, cavitazione, perdite dominate).
- Layout `wide`, 5 metriche sintetiche, 3 tab (Risultati / Grafici / Note tecniche).
- Grafici Plotly: barchart perdite e curva TDH(Q) con marcatura del punto di lavoro.
- Export CSV (dettaglio perdite + curva TDH).

## Struttura del progetto

```text
Pompe/
├── app.py           # UI Streamlit (nessuna formula)
├── src.py           # Calcoli idraulici, TDH, curva sistema, validazione, commenti
├── requirements.txt # streamlit, numpy, pandas, plotly
├── readme.md        # questo file
└── prompt.txt       # prompt per iterazioni future
```

## Input

| Sezione | Parametro | Unità |
|---------|-----------|-------|
| Fluido  | ρ, ν | kg/m³, m²/s |
| Portata | Q | m³/s |
| Aspirazione | L, D, ε, K_tot, z | m, m, m, –, m |
| Mandata | L, D, ε, K_tot, z | m, m, m, –, m |
| Pompa | η | – |

## Output principali

- Prevalenza totale `TDH [m]`, potenza `P [kW]`
- Velocità `V_asp` e `V_mand [m/s]`, numeri di Reynolds `Re`
- Fattori di attrito `f` (Darcy-Weisbach)
- Breakdown di tutti i contributi alla prevalenza
- Curva sistema `TDH(Q)`
- CSV dettaglio perdite e curva TDH

## Avvio rapido

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Estensioni future consigliate

- Analisi NPSH: confronto NPSH_disponibile vs NPSH_richiesto
- Curva caratteristica pompa (import da CSV o parametrizzazione) → punto di lavoro
- Più pompe in serie o parallelo
- Transitori (colpo d'ariete semplificato)
- Report PDF
