# -*- coding: utf-8 -*-
import json
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from src import (DatiCondotta, valida_dati, stato_condotta,
                 curva_riempimento, tabella_passaggi, verifiche_idrauliche,
                 confronto_diametri, genera_pdf, commenti_progettuali,
                 MATERIALI_MANNING, V_MIN_AUTOCIRCOLANTE, classifica_moto)

st.set_page_config(page_title="Condotte - Grado di riempimento", layout="wide")
st.title("Condotta circolare - Moto uniforme (Manning)")
st.caption("Calcolo professionale: grado di riempimento, verifiche normative, confronto DN standard, report PDF.")

# ---------------------------------------------------------------------------
# Defaults e session state
# ---------------------------------------------------------------------------
_DEFAULTS = {
    "cond_D": 1.0, "cond_n": 0.013, "cond_S": 0.003, "cond_Q": 0.5,
    "cond_mat": "Calcestruzzo liscio (prefabbricato)",
    "cond_uso": "Fognatura nera (acque reflue)",
}
for k, v in _DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---------------------------------------------------------------------------
# Sidebar - Input
# ---------------------------------------------------------------------------
with st.sidebar:
    with st.expander("Salva / Carica parametri", expanded=False):
        uploaded = st.file_uploader("Carica parametri (JSON)", type=["json"], key="cond_upload")
        if uploaded is not None:
            try:
                loaded = json.loads(uploaded.read())
                if st.button("Applica parametri caricati", key="cond_apply"):
                    for k in _DEFAULTS:
                        if k in loaded:
                            st.session_state[k] = loaded[k]
                    st.rerun()
                st.caption(f"File: {uploaded.name}")
            except Exception:
                st.error("File JSON non valido.")
        params_json = json.dumps({k: st.session_state[k] for k in _DEFAULTS}, indent=2).encode()
        st.download_button("Scarica parametri JSON", params_json,
                           "condotta_parametri.json", "application/json")

    st.divider()
    st.header("Materiale e uso")
    lista_mat = list(MATERIALI_MANNING.keys())
    idx_mat = lista_mat.index(st.session_state["cond_mat"]) if st.session_state["cond_mat"] in lista_mat else 2
    materiale = st.selectbox("Materiale condotta", lista_mat, index=idx_mat, key="cond_mat")
    n_tipico = MATERIALI_MANNING[materiale]["n"]
    V_max_mat = MATERIALI_MANNING[materiale]["V_max"]
    st.info(f"n tipico: {n_tipico:.3f}  |  V_max: {V_max_mat:.1f} m/s\n\n_{MATERIALI_MANNING[materiale]['note']}_")

    lista_uso = list(V_MIN_AUTOCIRCOLANTE.keys())
    idx_uso = lista_uso.index(st.session_state["cond_uso"]) if st.session_state["cond_uso"] in lista_uso else 0
    uso = st.selectbox("Tipo di uso / servizio", lista_uso, index=idx_uso, key="cond_uso")

    st.divider()
    st.header("Parametri condotta")
    D = st.number_input("Diametro D [m]", min_value=0.05, step=0.05, key="cond_D")
    n = st.number_input("Manning n [-]", min_value=0.005, max_value=0.05,
                        step=0.001, format="%.3f", key="cond_n")
    S = st.number_input("Pendenza idraulica S [-]", min_value=1e-6,
                        step=0.0005, format="%.6f", key="cond_S")
    Q = st.number_input("Portata Q [m\u00b3/s]", min_value=0.001, step=0.01, key="cond_Q")

    if st.button("Usa n tipico del materiale"):
        st.session_state["cond_n"] = n_tipico
        st.rerun()

# ---------------------------------------------------------------------------
# Calcolo
# ---------------------------------------------------------------------------
dati = DatiCondotta(D=D, n=n, S=S, Q=Q)
errori = valida_dati(dati)
if errori:
    for e in errori:
        st.error(e)
    st.stop()

res = stato_condotta(dati)
df_curva = curva_riempimento(D, n, S)
df_pass = tabella_passaggi(dati, res)
df_ver = verifiche_idrauliche(dati, res, materiale, uso)
df_dn = confronto_diametri(n, S, Q, D, materiale, uso)
note = commenti_progettuali(dati, res, materiale, uso)

# ---------------------------------------------------------------------------
# Indicatori sintetici
# ---------------------------------------------------------------------------
n_ok = (df_ver["Esito"] == "OK").sum()
n_att = (df_ver["Esito"] == "ATTENZIONE").sum()
n_no = (df_ver["Esito"] == "NON OK").sum()

col1, col2, col3, col4, col5, col6 = st.columns(6)
col1.metric("Tirante y [m]", f"{res['y']:.4f}")
col2.metric("Grado riempimento \u03b2", f"{res['beta']:.3f}",
            delta="OK" if res['beta'] <= 0.80 else "ALTO",
            delta_color="normal" if res['beta'] <= 0.80 else "inverse")
col3.metric("Velocit\u00e0 V [m/s]", f"{res['V']:.3f}")
col4.metric("Froude Fr [-]", f"{res['Fr']:.3f}" if res['Fr'] == res['Fr'] else "n.d.")
col5.metric("Q / Q_max [-]", f"{Q / res['Q_max']:.3f}")
col6.metric("Verifiche OK / WARN / NO", f"{n_ok} / {n_att} / {n_no}",
            delta_color="off")

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab1, tab2, tab3, tab4 = st.tabs(["Risultati", "Grafici", "Verifiche avanzate", "Note tecniche"])

with tab1:
    st.subheader("Passaggi di calcolo (passo per passo)")
    st.dataframe(df_pass, use_container_width=True, hide_index=True)

    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("**Input**")
        st.markdown(f"- Diametro D: **{D:.3f} m** (DN {int(round(D*1000))} mm)")
        st.markdown(f"- Manning n: **{n:.4f}** — materiale: *{materiale}*")
        st.markdown(f"- Pendenza S: **{S:.6f}**")
        st.markdown(f"- Portata Q: **{Q:.5f} m\u00b3/s**")
    with col_b:
        st.markdown("**Risultati**")
        st.markdown(f"- Tirante y: **{res['y']:.5f} m**")
        st.markdown(f"- \u03b2 = y/D: **{res['beta']:.4f}**")
        st.markdown(f"- V media: **{res['V']:.4f} m/s**")
        st.markdown(f"- Froude Fr: **{res['Fr']:.4f}** — {classifica_moto(res['Fr'])}")
        st.markdown(f"- Energia specifica E: **{res['E']:.5f} m**")
        st.markdown(f"- Rh: **{res['Rh']:.5f} m**")
        st.markdown(f"- Q_full: **{res['Q_full']:.5f} m\u00b3/s**")
        st.markdown(f"- Q_max (\u03b2\u22480.94): **{res['Q_max']:.5f} m\u00b3/s**")

    st.divider()
    col_dl1, col_dl2, col_dl3 = st.columns(3)
    with col_dl1:
        st.download_button("Scarica passaggi CSV",
                           df_pass.to_csv(index=False).encode("utf-8"),
                           "condotta_passaggi.csv", "text/csv")
    with col_dl2:
        st.download_button("Scarica curva idraulica CSV",
                           df_curva.to_csv(index=False).encode("utf-8"),
                           "condotta_curva.csv", "text/csv")
    with col_dl3:
        try:
            pdf_bytes = genera_pdf(dati, res, note, materiale, uso)
            st.download_button("Scarica Report PDF", pdf_bytes,
                               "condotta_report.pdf", "application/pdf")
        except ImportError:
            st.warning("fpdf2 non installato: pip install fpdf2")

with tab2:
    st.subheader("Curva idraulica della sezione circolare")
    fig_q = px.line(df_curva, x="y/D [-]", y="Q [m3/s]",
                    title="Portata Q in funzione del grado di riempimento \u03b2 = y/D")
    fig_q.add_vline(x=res["beta"], line_dash="dash", line_color="red",
                    annotation_text=f"\u03b2={res['beta']:.3f}", annotation_position="top right")
    fig_q.add_hline(y=Q, line_dash="dot", line_color="orange",
                    annotation_text=f"Q={Q:.3f} m\u00b3/s", annotation_position="bottom right")
    fig_q.add_vline(x=0.80, line_dash="dot", line_color="green",
                    annotation_text="\u03b2=0.80 (limite)", annotation_position="top left")
    fig_q.update_layout(xaxis_title="y/D [-]", yaxis_title="Q [m\u00b3/s]")
    st.plotly_chart(fig_q, use_container_width=True)

    fig_v = px.line(df_curva, x="y/D [-]", y="V [m/s]",
                    title="Velocit\u00e0 media V in funzione del grado di riempimento")
    fig_v.add_vline(x=res["beta"], line_dash="dash", line_color="red",
                    annotation_text=f"\u03b2={res['beta']:.3f}", annotation_position="top right")
    fig_v.add_hline(y=V_MIN_AUTOCIRCOLANTE.get(uso, 0.6), line_dash="dot", line_color="green",
                    annotation_text=f"V_min={V_MIN_AUTOCIRCOLANTE.get(uso,0.6):.2f} m/s",
                    annotation_position="bottom right")
    fig_v.update_layout(xaxis_title="y/D [-]", yaxis_title="V [m/s]")
    st.plotly_chart(fig_v, use_container_width=True)

    st.subheader("Numero di Froude in funzione del grado di riempimento")
    df_fr = df_curva.dropna(subset=["Fr [-]"])
    fig_fr = px.line(df_fr, x="y/D [-]", y="Fr [-]",
                     title="Numero di Froude lungo la curva di riempimento")
    fig_fr.add_hline(y=1.0, line_dash="dash", line_color="red",
                     annotation_text="Fr = 1 (critico)", annotation_position="bottom right")
    fig_fr.update_layout(xaxis_title="y/D [-]", yaxis_title="Fr [-]")
    st.plotly_chart(fig_fr, use_container_width=True)

with tab3:
    st.subheader("Verifiche normative")
    _colori_esito = {"OK": "background-color: #d4edda",
                     "ATTENZIONE": "background-color: #fff3cd",
                     "NON OK": "background-color: #f8d7da",
                     "INFO": "background-color: #d1ecf1",
                     "SUPERCRITICO": "background-color: #f8d7da",
                     "CRITICO": "background-color: #fff3cd"}

    def _colora_righe(row):
        colore = _colori_esito.get(row["Esito"], "")
        return [colore] * len(row)

    styled = df_ver.style.apply(_colora_righe, axis=1)
    st.dataframe(styled, use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("Confronto diametri DN standard")
    st.caption(f"Diametri nell'intervallo 0.5\u00d7D \u00f7 2.5\u00d7D con stesso n={n:.3f}, S={S:.5f}, Q={Q:.4f} m\u00b3/s")

    def _colora_dn(row):
        if row.get("Nota", "") and "<<<" in str(row.get("Nota", "")):
            return ["background-color: #cce5ff"] * len(row)
        g = row.get("Stato globale", "")
        c = _colori_esito.get(g, "")
        return [c] * len(row)

    styled_dn = df_dn.style.apply(_colora_dn, axis=1)
    st.dataframe(styled_dn, use_container_width=True, hide_index=True)

    st.divider()
    st.download_button("Scarica verifiche CSV",
                       df_ver.to_csv(index=False).encode("utf-8"),
                       "condotta_verifiche.csv", "text/csv")
    st.download_button("Scarica confronto DN CSV",
                       df_dn.to_csv(index=False).encode("utf-8"),
                       "condotta_confronto_dn.csv", "text/csv")

with tab4:
    st.subheader("Note tecniche e commenti di progetto")
    for item in note:
        st.markdown(f"- {item}")
    with st.expander("Legenda verifiche e riferimenti normativi"):
        st.markdown("""
**Riferimenti normativi utilizzati:**
- **DIN EN 1671:2011** — Fognature sotto pressione: limite \u03b2 \u2264 0.80, V_min autocircolante.
- **EN 476:2011** — Requisiti generali componenti di sistemi fognari.
- **ASCE MOP 36** — Progetto delle reti di fognatura: velocita minima 0.60 m/s.
- **UNI EN 752** — Sistemi di drenaggio delle acque reflue.

**Interpretazione del Numero di Froude:**
- Fr < 1.0 \u2192 Moto subcritico (tranquillo) — condizione preferita per fognature.
- Fr \u2248 1.0 \u2192 Condizione critica — instabilita e transizioni di regime.
- Fr > 1.0 \u2192 Moto supercritico (veloce) — rischio risalti idraulici.

**Grado di riempimento \u03b2:**
- \u03b2 \u2264 0.80: condizione normale (DIN EN 1671).
- 0.80 < \u03b2 \u2264 0.94: alta occupazione — verificare ventilazione.
- \u03b2 > 0.94: oltre il massimo teorico di portata — sezione insufficiente.
        """)
