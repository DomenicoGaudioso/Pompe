# -*- coding: utf-8 -*-
import json
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from src import (Fluido, Linea, valida_dati, tdh_pump, potenza_pompa,
                 breakdown_perdite, curva_tdh_vs_Q, tabella_passaggi,
                 genera_pdf, commenti_progettuali,
                 verifiche_pompa, npsh_disponibile, velocita_specifica_ns,
                 tipo_pompa_da_ns, sovrapressione_joukowsky, costo_energetico_annuo,
                 FLUIDI_PRESET, RUGOSITA_TUBI, VELOCITA_ONDA_MATERIALE)

st.set_page_config(page_title="Pompe - Prevalenza TDH", layout="wide")
st.title("Calcolo della prevalenza (TDH) e potenza della pompa")
st.caption("Software professionale: Darcy-Weisbach, NPSH, velocit\u00e0 specifica, colpo d'ariete, verifiche normative, report PDF.")

# ---------------------------------------------------------------------------
# Defaults e session state
# ---------------------------------------------------------------------------
_DEFAULTS = {
    "pmp_fluido_preset": "Acqua 20 deg C (std)",
    "pmp_rho": 998.2, "pmp_nu": 1.004e-6, "pmp_p_vap": 2338.0,
    "pmp_Q": 0.05,
    "pmp_mat_asp": "Acciaio senza saldatura (lappato)",
    "pmp_Ls": 8.0, "pmp_Ds": 0.15, "pmp_epss": 0.000046,
    "pmp_Ks": 2.0, "pmp_zs": 0.0,
    "pmp_mat_mand": "Acciaio senza saldatura (lappato)",
    "pmp_Ld": 40.0, "pmp_Dd": 0.10, "pmp_epsd": 0.000046,
    "pmp_Kd": 6.0, "pmp_zd": 12.0,
    "pmp_eta": 0.75,
    "pmp_n_rpm": 1450.0,
    "pmp_z_serb": 0.0,
    "pmp_npsh_r": 3.0,
    "pmp_mat_onda": "Acciaio",
    "pmp_ore_annue": 2000.0,
    "pmp_costo_kwh": 0.15,
    "pmp_eta_motore": 0.92,
}
for k, v in _DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---------------------------------------------------------------------------
# Sidebar - Input
# ---------------------------------------------------------------------------
with st.sidebar:
    with st.expander("Salva / Carica parametri", expanded=False):
        uploaded = st.file_uploader("Carica parametri (JSON)", type=["json"], key="pmp_upload")
        if uploaded is not None:
            try:
                loaded = json.loads(uploaded.read())
                if st.button("Applica parametri caricati", key="pmp_apply"):
                    for k in _DEFAULTS:
                        if k in loaded:
                            st.session_state[k] = loaded[k]
                    st.rerun()
                st.caption(f"File: {uploaded.name}")
            except Exception:
                st.error("File JSON non valido.")
        params_json = json.dumps({k: st.session_state[k] for k in _DEFAULTS}, indent=2).encode()
        st.download_button("Scarica parametri JSON", params_json,
                           "pompa_parametri.json", "application/json")

    st.divider()
    st.header("Fluido")
    lista_preset = list(FLUIDI_PRESET.keys())
    idx_preset = lista_preset.index(st.session_state["pmp_fluido_preset"]) \
        if st.session_state["pmp_fluido_preset"] in lista_preset else 3
    fluido_preset = st.selectbox("Fluido / temperatura", lista_preset,
                                 index=idx_preset, key="pmp_fluido_preset")
    preset_data = FLUIDI_PRESET[fluido_preset]
    st.info(f"\u03c1={preset_data['rho']:.1f} kg/m\u00b3  |  "
            f"\u03bd={preset_data['nu']*1e6:.3f}\u00d710\u207b\u2076 m\u00b2/s  |  "
            f"p_vap={preset_data['p_vap']:.0f} Pa\n\n_{preset_data['note']}_")
    if st.button("Usa propriet\u00e0 preset"):
        st.session_state["pmp_rho"] = preset_data["rho"]
        st.session_state["pmp_nu"] = preset_data["nu"]
        st.session_state["pmp_p_vap"] = preset_data["p_vap"]
        st.rerun()

    rho = st.number_input("Densit\u00e0 \u03c1 [kg/m\u00b3]", min_value=600.0, max_value=1200.0,
                          step=10.0, key="pmp_rho")
    nu = st.number_input("Viscosit\u00e0 cinematica \u03bd [m\u00b2/s]", min_value=0.2e-6,
                         max_value=50.0e-6, step=0.1e-6, format="%.7f", key="pmp_nu")
    p_vap = st.number_input("Pressione di vapore p_vap [Pa]", min_value=0.0,
                             step=100.0, key="pmp_p_vap")
    Q = st.number_input("Portata Q [m\u00b3/s]", min_value=0.001, step=0.001, key="pmp_Q")

    st.divider()
    st.header("Linea di aspirazione")
    lista_mat = list(RUGOSITA_TUBI.keys())
    idx_asp = lista_mat.index(st.session_state["pmp_mat_asp"]) \
        if st.session_state["pmp_mat_asp"] in lista_mat else 1
    mat_asp = st.selectbox("Materiale tubo asp.", lista_mat, index=idx_asp, key="pmp_mat_asp")
    eps_asp_tipico = RUGOSITA_TUBI[mat_asp]
    st.caption(f"\u03b5 tipica: {eps_asp_tipico*1000:.4f} mm")
    if st.button("Usa \u03b5 tipica asp."):
        st.session_state["pmp_epss"] = eps_asp_tipico
        st.rerun()
    Ls = st.number_input("Lunghezza Ls [m]", min_value=0.1, step=0.5, key="pmp_Ls")
    Ds = st.number_input("Diametro Ds [m]", min_value=0.02, step=0.01, key="pmp_Ds")
    epss = st.number_input("Scabrezza \u03b5s [m]", min_value=0.0,
                            step=0.000001, format="%.6f", key="pmp_epss")
    Ks = st.number_input("Somma K concentrati (asp.)", min_value=0.0, step=0.1, key="pmp_Ks")
    zs = st.number_input("Quota pompa z_s [m]", step=0.1, key="pmp_zs")

    st.divider()
    st.header("Linea di mandata")
    idx_mand = lista_mat.index(st.session_state["pmp_mat_mand"]) \
        if st.session_state["pmp_mat_mand"] in lista_mat else 1
    mat_mand = st.selectbox("Materiale tubo mand.", lista_mat, index=idx_mand, key="pmp_mat_mand")
    eps_mand_tipico = RUGOSITA_TUBI[mat_mand]
    st.caption(f"\u03b5 tipica: {eps_mand_tipico*1000:.4f} mm")
    if st.button("Usa \u03b5 tipica mand."):
        st.session_state["pmp_epsd"] = eps_mand_tipico
        st.rerun()
    Ld = st.number_input("Lunghezza Ld [m]", min_value=0.1, step=1.0, key="pmp_Ld")
    Dd = st.number_input("Diametro Dd [m]", min_value=0.02, step=0.01, key="pmp_Dd")
    epsd = st.number_input("Scabrezza \u03b5d [m]", min_value=0.0,
                            step=0.000001, format="%.6f", key="pmp_epsd")
    Kd = st.number_input("Somma K concentrati (mand.)", min_value=0.0, step=0.1, key="pmp_Kd")
    zd = st.number_input("Quota scarico z_d [m]", step=0.1, key="pmp_zd")

    st.divider()
    st.header("Pompa e motore")
    eta = st.number_input("Rendimento pompa \u03b7 [-]", min_value=0.2, max_value=0.98,
                          step=0.01, key="pmp_eta")
    n_rpm = st.number_input("Velocit\u00e0 di rotazione n [giri/min]",
                             min_value=100.0, max_value=6000.0, step=50.0, key="pmp_n_rpm")

    st.divider()
    st.header("NPSH")
    z_serb = st.number_input("Quota pelo libero serbatoio z_serb [m]",
                              step=0.1, key="pmp_z_serb")
    npsh_r = st.number_input("NPSH richiesto (dalla pompa) [m]",
                              min_value=0.1, step=0.1, key="pmp_npsh_r")

    st.divider()
    st.header("Colpo d'ariete")
    lista_onda = list(VELOCITA_ONDA_MATERIALE.keys())
    idx_onda = lista_onda.index(st.session_state["pmp_mat_onda"]) \
        if st.session_state["pmp_mat_onda"] in lista_onda else 0
    mat_onda = st.selectbox("Materiale condotta (velocit\u00e0 d'onda)",
                             lista_onda, index=idx_onda, key="pmp_mat_onda")
    a_onda = VELOCITA_ONDA_MATERIALE[mat_onda]
    st.caption(f"a_onda = {a_onda:.0f} m/s")

    st.divider()
    st.header("Costo energetico")
    ore_annue = st.number_input("Ore di funzionamento annuo [h/anno]",
                                 min_value=100.0, max_value=8760.0, step=100.0, key="pmp_ore_annue")
    costo_kwh = st.number_input("Costo energia [EUR/kWh]", min_value=0.01,
                                 max_value=1.0, step=0.01, format="%.3f", key="pmp_costo_kwh")
    eta_motore = st.number_input("Rendimento motore elettrico [-]", min_value=0.5,
                                  max_value=0.99, step=0.01, key="pmp_eta_motore")

# ---------------------------------------------------------------------------
# Calcolo
# ---------------------------------------------------------------------------
fluido = Fluido(rho=rho, nu=nu, p_vap=p_vap)
suction = Linea(L=Ls, D=Ds, eps=epss, K_tot=Ks, z=zs)
discharge = Linea(L=Ld, D=Dd, eps=epsd, K_tot=Kd, z=zd)

errori = valida_dati(Q, suction, discharge, fluido, eta)
if errori:
    for e in errori:
        st.error(e)
    st.stop()

res = tdh_pump(Q=Q, suction=suction, discharge=discharge, fluido=fluido)
P_kW = potenza_pompa(Q=Q, H=res["H"], rho=fluido.rho, eta=eta) / 1000.0
df_perdite = breakdown_perdite(res)
df_curva = curva_tdh_vs_Q(suction, discharge, fluido, Q)
df_pass = tabella_passaggi(Q, suction, discharge, fluido, eta, res)
note = commenti_progettuali(Q, suction, discharge, fluido, res, eta)

hf_asp_tot = res["hf_s"] + res["hK_s"]
npsh_d = npsh_disponibile(z_serb, zs, hf_asp_tot, fluido.rho, p_vap=fluido.p_vap)
Ns = velocita_specifica_ns(n_rpm, Q, res["H"])
tipo_pompa = tipo_pompa_da_ns(Ns)
ariete = sovrapressione_joukowsky(res["V_d"], fluido.rho, a_onda)
costi = costo_energetico_annuo(P_kW, ore_annue, costo_kwh, eta_motore)
df_ver = verifiche_pompa(Q, suction, discharge, fluido, eta, res,
                          n_rpm=n_rpm, z_serbatoio=z_serb,
                          npsh_richiesto=npsh_r)

# ---------------------------------------------------------------------------
# Indicatori sintetici
# ---------------------------------------------------------------------------
n_ok  = (df_ver["Esito"] == "OK").sum()
n_att = (df_ver["Esito"] == "ATTENZIONE").sum()
n_no  = (df_ver["Esito"] == "NON OK").sum()

col1, col2, col3, col4, col5, col6, col7 = st.columns(7)
col1.metric("TDH [m]", f"{res['H']:.2f}")
col2.metric("Potenza asse [kW]", f"{P_kW:.2f}")
col3.metric("NPSH disp. [m]", f"{npsh_d:.2f}",
            delta="OK" if npsh_d >= npsh_r else "CRIT.",
            delta_color="normal" if npsh_d >= npsh_r else "inverse")
col4.metric("Ns [-]", f"{Ns:.0f}" if Ns == Ns else "n.d.")
col5.metric("V asp [m/s]", f"{res['V_s']:.2f}")
col6.metric("V mand [m/s]", f"{res['V_d']:.2f}")
col7.metric("Verif. OK / WARN / NO", f"{n_ok} / {n_att} / {n_no}", delta_color="off")

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
        st.markdown("**Input fluido e linee**")
        st.markdown(f"- Fluido: **{fluido_preset}**  (\u03c1={rho:.1f} kg/m\u00b3, \u03bd={nu*1e6:.3f}\u00d710\u207b\u2076 m\u00b2/s)")
        st.markdown(f"- Portata Q: **{Q:.5f} m\u00b3/s**")
        st.markdown(f"- Asp.: L={Ls:.1f} m, D={Ds:.3f} m, \u03b5={epss*1e3:.4f} mm, K={Ks:.2f}")
        st.markdown(f"- Mand.: L={Ld:.1f} m, D={Dd:.3f} m, \u03b5={epsd*1e3:.4f} mm, K={Kd:.2f}")
        st.markdown(f"- Rendimento pompa \u03b7: **{eta:.2f}**")
    with col_b:
        st.markdown("**Risultati principali**")
        st.markdown(f"- Re asp.: **{res['Re_s']:.0f}** | f_s: **{res['f_s']:.5f}**")
        st.markdown(f"- Re mand.: **{res['Re_d']:.0f}** | f_d: **{res['f_d']:.5f}**")
        st.markdown(f"- Perdite distribuite tot.: **{res['hf_s']+res['hf_d']:.3f} m**")
        st.markdown(f"- Perdite concentrate tot.: **{res['hK_s']+res['hK_d']:.3f} m**")
        st.markdown(f"- Dislivello geodesico: **{res['dZ']:.3f} m**")
        st.markdown(f"- TDH: **{res['H']:.4f} m**")
        st.markdown(f"- Potenza asse: **{P_kW:.3f} kW**")
        st.markdown(f"- Costo annuo: **{costi['Costo annuo [EUR]']:,.0f} EUR/anno** "
                    f"({ore_annue:.0f} h/anno, {costo_kwh:.3f} EUR/kWh)")

    st.divider()
    st.subheader("Dettaglio componenti della prevalenza")
    st.dataframe(df_perdite, use_container_width=True, hide_index=True)

    st.divider()
    col_dl1, col_dl2, col_dl3, col_dl4 = st.columns(4)
    with col_dl1:
        st.download_button("Scarica passaggi CSV",
                           df_pass.to_csv(index=False).encode("utf-8"),
                           "pompa_passaggi.csv", "text/csv")
    with col_dl2:
        st.download_button("Scarica dettaglio perdite CSV",
                           df_perdite.to_csv(index=False).encode("utf-8"),
                           "pompa_perdite.csv", "text/csv")
    with col_dl3:
        st.download_button("Scarica curva TDH(Q) CSV",
                           df_curva.to_csv(index=False).encode("utf-8"),
                           "pompa_tdh_Q.csv", "text/csv")
    with col_dl4:
        try:
            pdf_bytes = bytes(genera_pdf(Q, suction, discharge, fluido, eta, res, note,
                                        n_rpm=n_rpm, z_serbatoio=z_serb,
                                        npsh_richiesto=npsh_r,
                                        ore_annue=ore_annue, costo_kwh=costo_kwh,
                                        eta_motore=eta_motore))
            st.download_button("Scarica Report PDF", pdf_bytes,
                               "pompa_report.pdf", "application/pdf")
        except ImportError:
            st.warning("fpdf2 non installato. Eseguire: pip install fpdf2")

with tab2:
    st.subheader("Distribuzione delle perdite")
    df_bar = df_perdite[df_perdite["Componente"] != "TDH totale"].copy()
    fig_bar = px.bar(df_bar, x="Componente", y="Valore [m]",
                     title="Contributi alla prevalenza totale",
                     color="Componente", text_auto=".2f")
    fig_bar.update_layout(showlegend=False, xaxis_title="", yaxis_title="m c.a.")
    st.plotly_chart(fig_bar, use_container_width=True)

    st.subheader("Curva sistema TDH(Q)")
    fig_tdh = px.line(df_curva, x="Q [m3/s]", y="TDH [m]",
                      title="Curva del sistema TDH in funzione della portata")
    fig_tdh.add_vline(x=Q, line_dash="dash", line_color="red",
                      annotation_text=f"Q={Q:.3f} m\u00b3/s", annotation_position="top right")
    fig_tdh.add_hline(y=res["H"], line_dash="dot", line_color="orange",
                      annotation_text=f"TDH={res['H']:.2f} m", annotation_position="bottom right")
    fig_tdh.update_layout(xaxis_title="Q [m\u00b3/s]", yaxis_title="TDH [m]")
    st.plotly_chart(fig_tdh, use_container_width=True)

    st.subheader("Scomposizione perdite lungo la curva sistema")
    fig_stack = go.Figure()
    fig_stack.add_trace(go.Scatter(x=df_curva["Q [m3/s]"], y=df_curva["dZ [m]"],
                                   name="dZ (geodesico)", stackgroup="one", mode="none",
                                   fillcolor="rgba(31,119,180,0.5)"))
    fig_stack.add_trace(go.Scatter(x=df_curva["Q [m3/s]"], y=df_curva["Perdite tot [m]"],
                                   name="Perdite totali", stackgroup="one", mode="none",
                                   fillcolor="rgba(255,127,14,0.5)"))
    fig_stack.update_layout(title="Scomposizione TDH: geodesico + perdite",
                             xaxis_title="Q [m\u00b3/s]", yaxis_title="m c.a.")
    st.plotly_chart(fig_stack, use_container_width=True)

with tab3:
    st.subheader("Verifiche normative pompa")
    _colori = {"OK": "background-color: #d4edda",
               "ATTENZIONE": "background-color: #fff3cd",
               "NON OK": "background-color: #f8d7da",
               "INFO": "background-color: #d1ecf1"}

    def _colora_righe(row):
        c = _colori.get(row["Esito"], "")
        return [c] * len(row)

    styled = df_ver.style.apply(_colora_righe, axis=1)
    st.dataframe(styled, use_container_width=True, hide_index=True)

    st.divider()
    col_npsh, col_ns = st.columns(2)

    with col_npsh:
        st.subheader("NPSH - Analisi cavitazione")
        st.markdown(f"- **NPSH disponibile**: {npsh_d:.3f} m")
        st.markdown(f"- **NPSH richiesto**: {npsh_r:.2f} m")
        margine = npsh_d - npsh_r
        stato_npsh = "OK" if margine >= 0.5 else ("ATTENZIONE" if margine >= 0 else "NON OK")
        st.markdown(f"- **Margine NPSH**: {margine:.3f} m  \u2192  **{stato_npsh}**")
        st.markdown(f"- Perdite aspirazione (hf+hK): {hf_asp_tot:.4f} m")
        st.markdown(f"- Quota pelo libero z_serb: {z_serb:.2f} m  |  z_pompa: {zs:.2f} m")
        st.caption("Rif.: ISO 9906:2012, EN 12845. Margine minimo consigliato: 0.5 m.")

    with col_ns:
        st.subheader("Velocit\u00e0 specifica e tipo di pompa")
        st.markdown(f"- **Ns** = n \u00b7 Q\u00b0\u00b7\u2075 / H\u2070\u00b7\u2077\u2075 = **{Ns:.1f}** (rpm, m\u00b3/s, m)")
        st.markdown(f"- Tipo pompa: **{tipo_pompa}**")
        st.caption("Rif.: Stepanoff, KSB Pump Handbook.")

    st.divider()
    col_ariete, col_costo = st.columns(2)

    with col_ariete:
        st.subheader("Colpo d'ariete (Joukowsky)")
        st.markdown(f"- Materiale condotta: **{mat_onda}**")
        st.markdown(f"- Velocit\u00e0 d'onda a = **{a_onda:.0f} m/s**")
        st.markdown(f"- Velocit\u00e0 flusso in mandata V_d = **{res['V_d']:.3f} m/s**")
        st.markdown(f"- Sovrapressione \u0394P = **{ariete['dP [bar]']:.2f} bar** "
                    f"= **{ariete['dP [m c.a.]']:.1f} m c.a.**")
        st.caption("Rif.: EN 805:2000, Wylie & Streeter (1993). "
                   "Formula: \u0394P = \u03c1 \u00b7 a \u00b7 V (chiusura istantanea).")

    with col_costo:
        st.subheader("Costo energetico annuo")
        st.markdown(f"- Potenza idraulica: **{costi['P_idraulica [kW]']:.3f} kW**")
        st.markdown(f"- Potenza asse pompa: **{costi['P_asse [kW]']:.3f} kW**")
        st.markdown(f"- Potenza elettrica assorbita: **{costi['P_elettrica [kW]']:.3f} kW**")
        st.markdown(f"  (\u03b7 motore = {eta_motore:.2f})")
        st.markdown(f"- Ore/anno: **{ore_annue:.0f} h**")
        st.markdown(f"- Energia annua: **{costi['Energia annua [kWh]']:,.0f} kWh/anno**")
        st.markdown(f"- **Costo annuo: {costi['Costo annuo [EUR]']:,.0f} EUR/anno**")
        st.caption(f"Tariffa: {costo_kwh:.3f} EUR/kWh")

    st.divider()
    st.download_button("Scarica verifiche CSV",
                       df_ver.to_csv(index=False).encode("utf-8"),
                       "pompa_verifiche.csv", "text/csv")

with tab4:
    st.subheader("Note tecniche e commenti di progetto")
    for item in note:
        st.markdown(f"- {item}")
    with st.expander("Legenda e riferimenti normativi"):
        st.markdown("""
**TDH (Total Dynamic Head):**
Prevalenza totale = dislivello geodesico + perdite distribuite + perdite concentrate + variazione cinetica.

**Fattore di attrito Darcy-Weisbach:**
- Re < 2300: regime laminare, f = 64/Re
- Re \u2265 2300: regime turbolento, formula di Swamee-Jain (approssimazione Colebrook-White)

**NPSH disponibile (ISO 9906):**
NPSH_d = (p_atm - p_vap) / (\u03c1g) + (z_serb - z_pompa) - hf_asp

**Velocit\u00e0 specifica Ns:**
Ns = n \u00b7 Q\u00b0\u00b7\u2075 / H\u2070\u00b7\u2077\u2075 [rpm, m\u00b3/s, m]
- Ns < 400: centrifuga radiale (alta prevalenza)
- 400 \u2264 Ns < 1200: centrifuga mista (uso comune)
- 1200 \u2264 Ns < 3000: flusso misto (alta portata)
- Ns \u2265 3000: elicoidale / assiale

**Colpo d'ariete (Joukowsky):**
\u0394P = \u03c1 \u00b7 a \u00b7 V  (chiusura istantanea valvola)
La velocit\u00e0 d'onda a dipende dal materiale e dallo spessore della condotta.

**Riferimenti normativi:**
- ISO 9906:2012 — Prove di accettazione pompe rotodinamiche
- EN 12845 — Sistemi fissi antincendio (portata, NPSH)
- EN 805:2000 — Sistemi di adduzione acqua (colpo d'ariete)
- UNI EN ISO 17769 — Pompe a liquido e unita di pompaggio
        """)
