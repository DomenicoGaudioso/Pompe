# -*- coding: utf-8 -*-
"""
POMPE - Calcolo della prevalenza (TDH) e potenza richiesta.
Versione professionale:
- Perdite distribuite (Darcy-Weisbach, Swamee-Jain; laminare: f = 64/Re)
- Perdite concentrate (somma K)
- NPSH disponibile e verifica cavitazione
- Velocita specifica Ns e classificazione tipo pompa
- Colpo d'ariete: sovrapressione Joukowsky
- Costo energetico annuo
- Curva sistema TDH(Q)
- Tabella passaggi di calcolo
- Verifiche normative complete
- Generazione report PDF
"""
from __future__ import annotations

import datetime
import math
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

G = 9.81
P_ATM_STD = 101_325.0  # Pa (pressione atmosferica standard)


# ---------------------------------------------------------------------------
# Librerie preset fluidi (rho, nu, p_vap a T di riferimento)
# ---------------------------------------------------------------------------

FLUIDI_PRESET: Dict[str, dict] = {
    "Acqua 5 deg C": {
        "rho": 999.9, "nu": 1.518e-6, "p_vap": 872,
        "note": "Acqua di falda fredda",
    },
    "Acqua 10 deg C": {
        "rho": 999.7, "nu": 1.307e-6, "p_vap": 1228,
        "note": "Acqua tipica invernale",
    },
    "Acqua 15 deg C": {
        "rho": 999.1, "nu": 1.140e-6, "p_vap": 1705,
        "note": "Temperatura ambiente fresca",
    },
    "Acqua 20 deg C (std)": {
        "rho": 998.2, "nu": 1.004e-6, "p_vap": 2338,
        "note": "Condizione standard di riferimento",
    },
    "Acqua 40 deg C": {
        "rho": 992.2, "nu": 0.658e-6, "p_vap": 7384,
        "note": "Acqua calda (impianti termici)",
    },
    "Acqua 60 deg C": {
        "rho": 983.2, "nu": 0.474e-6, "p_vap": 19_940,
        "note": "Acqua calda sanitaria / riscaldamento",
    },
    "Acqua 80 deg C": {
        "rho": 971.8, "nu": 0.364e-6, "p_vap": 47_360,
        "note": "Impianti riscaldamento alta temperatura",
    },
    "Gasolio (diesel)": {
        "rho": 840.0, "nu": 3.5e-6, "p_vap": 100,
        "note": "Carburante industriale, stazioni di pompaggio",
    },
    "Olio minerale leggero": {
        "rho": 870.0, "nu": 20.0e-6, "p_vap": 10,
        "note": "Lubrificanti industriali, idraulica",
    },
}

# Rugosita assoluta per materiale tubo [m]
RUGOSITA_TUBI: Dict[str, float] = {
    "PVC / PE liscio": 0.000_002,
    "Acciaio senza saldatura (lappato)": 0.000_046,
    "Acciaio saldato": 0.000_046,
    "Ghisa sferoidale rivestita": 0.000_100,
    "Ghisa grigia vecchia": 0.000_500,
    "Calcestruzzo liscio": 0.000_300,
    "Calcestruzzo rugoso": 0.001_000,
    "Ferro zincato": 0.000_150,
}


# ---------------------------------------------------------------------------
# Dataclass input
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Fluido:
    rho: float = 998.2     # densita [kg/m³]
    nu: float = 1.004e-6   # viscosita cinematica [m²/s]
    p_vap: float = 2338.0  # pressione di vapore [Pa]


@dataclass(frozen=True)
class Linea:
    L: float       # lunghezza [m]
    D: float       # diametro interno [m]
    eps: float     # rugosita assoluta [m]
    K_tot: float   # somma coefficienti perdite concentrate [-]
    z: float       # quota di riferimento [m]


# ---------------------------------------------------------------------------
# Validazione
# ---------------------------------------------------------------------------

def valida_dati(Q: float, suction: Linea, discharge: Linea,
                fluido: Fluido, eta: float) -> List[str]:
    errori: List[str] = []
    if Q <= 0:
        errori.append("La portata Q deve essere positiva.")
    if fluido.rho <= 0:
        errori.append("La densita' del fluido deve essere positiva.")
    if fluido.nu <= 0:
        errori.append("La viscosita' cinematica deve essere positiva.")
    for nome, linea in [("aspirazione", suction), ("mandata", discharge)]:
        if linea.L <= 0:
            errori.append(f"La lunghezza della linea di {nome} deve essere positiva.")
        if linea.D <= 0:
            errori.append(f"Il diametro della linea di {nome} deve essere positivo.")
        if linea.eps < 0:
            errori.append(f"La rugosita' della linea di {nome} non puo' essere negativa.")
        if linea.K_tot < 0:
            errori.append(f"La somma K concentrati ({nome}) non puo' essere negativa.")
    if not (0 < eta <= 1.0):
        errori.append("Il rendimento eta deve essere compreso tra 0 e 1.")
    return errori


# ---------------------------------------------------------------------------
# Calcoli idraulici di base
# ---------------------------------------------------------------------------

def velocity(Q: float, D: float) -> float:
    return Q / (math.pi * D ** 2 / 4.0)


def reynolds(V: float, D: float, nu: float) -> float:
    if D <= 0 or nu <= 0:
        return float("nan")
    return V * D / nu


def friction_factor(eps_rel: float, Re: float) -> float:
    """f = 64/Re (laminare, Re<2300) oppure Swamee-Jain (turbolento)."""
    if Re <= 0:
        return float("nan")
    if Re < 2300:
        return 64.0 / Re
    return 0.25 / (math.log10(eps_rel / 3.7 + 5.74 / (Re ** 0.9)) ** 2)


def head_loss_friction(L: float, D: float, V: float, f: float) -> float:
    return f * (L / D) * (V * V) / (2.0 * G)


def head_loss_minor(K: float, V: float) -> float:
    return K * (V * V) / (2.0 * G)


# ---------------------------------------------------------------------------
# Calcolo TDH
# ---------------------------------------------------------------------------

def tdh_pump(Q: float, suction: Linea, discharge: Linea,
             fluido: Fluido, p_s: float = 0.0, p_d: float = 0.0) -> dict:
    V_s = velocity(Q, suction.D)
    V_d = velocity(Q, discharge.D)
    Re_s = reynolds(V_s, suction.D, fluido.nu)
    Re_d = reynolds(V_d, discharge.D, fluido.nu)
    f_s = friction_factor(suction.eps / suction.D, Re_s)
    f_d = friction_factor(discharge.eps / discharge.D, Re_d)
    hf_s = head_loss_friction(suction.L, suction.D, V_s, f_s)
    hf_d = head_loss_friction(discharge.L, discharge.D, V_d, f_d)
    hK_s = head_loss_minor(suction.K_tot, V_s)
    hK_d = head_loss_minor(discharge.K_tot, V_d)
    dZ = discharge.z - suction.z
    dVel = (V_d ** 2 - V_s ** 2) / (2.0 * G)
    dP = (p_d - p_s) / (fluido.rho * G)
    H = dZ + hf_s + hf_d + hK_s + hK_d + dVel + dP
    return dict(H=H, V_s=V_s, V_d=V_d, Re_s=Re_s, Re_d=Re_d,
                f_s=f_s, f_d=f_d, hf_s=hf_s, hf_d=hf_d,
                hK_s=hK_s, hK_d=hK_d, dZ=dZ, dVel=dVel, dP=dP)


def potenza_pompa(Q: float, H: float, rho: float = 998.2, eta: float = 0.75) -> float:
    if eta <= 0:
        return float("nan")
    return rho * G * Q * H / eta


# ---------------------------------------------------------------------------
# NPSH disponibile
# ---------------------------------------------------------------------------

def npsh_disponibile(z_serbatoio: float, z_pompa: float,
                     hf_asp_tot: float,
                     rho: float = 998.2,
                     p_atm: float = P_ATM_STD,
                     p_vap: float = 2338.0) -> float:
    """
    NPSH disponibile [m c.a.]:
    NPSH_d = (p_atm - p_vap) / (rho * g) + (z_serbatoio - z_pompa) - hf_asp
    dove hf_asp = perdite totali nel tratto di aspirazione.
    Riferimento: ISO 9906:2012, EN 12845.
    """
    p_term = (p_atm - p_vap) / (rho * G)
    z_term = z_serbatoio - z_pompa
    return p_term + z_term - hf_asp_tot


# ---------------------------------------------------------------------------
# Velocita specifica e tipo di pompa
# ---------------------------------------------------------------------------

def velocita_specifica_ns(n_rpm: float, Q: float, H: float) -> float:
    """
    Velocita specifica Ns = n * Q^0.5 / H^0.75
    [giri/min, m³/s, m] -> valore adimensionale o in rpm-m3/4/s1/2.
    Riferimento: Stepanoff, Kaplan, KSB Pump Handbook.
    """
    if H <= 0 or Q <= 0 or n_rpm <= 0:
        return float("nan")
    return n_rpm * (Q ** 0.5) / (H ** 0.75)


def tipo_pompa_da_ns(Ns: float) -> str:
    """Classificazione del tipo di pompa dalla velocita specifica Ns."""
    if math.isnan(Ns):
        return "n.d."
    if Ns < 400:
        return "Centrifuga radiale (alta prevalenza, bassa portata)"
    if Ns < 1200:
        return "Centrifuga mista (normale)"
    if Ns < 3000:
        return "Centrifuga a flusso misto (alta portata)"
    if Ns < 8000:
        return "Elicoidale / assiale (bassa prevalenza, alta portata)"
    return "Assiale pura (pompe fluviali / di drenaggio)"


# ---------------------------------------------------------------------------
# Colpo d'ariete (Joukowsky)
# ---------------------------------------------------------------------------

def sovrapressione_joukowsky(V: float, rho: float = 998.2,
                              a_onda: float = 1200.0) -> dict:
    """
    Sovrapressione da colpo d'ariete (chiusura istantanea valvola):
    dP = rho * a_onda * V  [Pa]
    Velocita d'onda a_onda dipende dal materiale e dalla condotta.
    Riferimento: EN 805:2000, Wylie & Streeter (1993).
    """
    dP_Pa = rho * a_onda * V
    dP_bar = dP_Pa / 1e5
    dP_mcol = dP_Pa / (rho * G)
    return {
        "dP [Pa]": dP_Pa,
        "dP [bar]": dP_bar,
        "dP [m c.a.]": dP_mcol,
        "a_onda [m/s]": a_onda,
    }


# Velocita d'onda tipiche per materiale [m/s]
VELOCITA_ONDA_MATERIALE: Dict[str, float] = {
    "Acciaio": 1300.0,
    "Ghisa sferoidale": 1200.0,
    "PVC rigido": 350.0,
    "PE/HDPE": 300.0,
    "Calcestruzzo armato": 1000.0,
    "Vetroresina (GRP)": 700.0,
}


# ---------------------------------------------------------------------------
# Costo energetico annuo
# ---------------------------------------------------------------------------

def costo_energetico_annuo(P_kW: float, ore_annue: float = 2000.0,
                            costo_kwh: float = 0.15,
                            rendimento_motore: float = 0.92) -> dict:
    """
    Costo energetico annuo considerando anche il rendimento del motore elettrico.
    P_assorbita_motore = P_idraulica / (eta_pompa gia' inclusa in P_kW) / eta_motore
    -> In realta P_kW e' gia' la potenza all'asse; divido per eta_motore per la potenza elettrica.
    """
    P_elettrica_kW = P_kW / rendimento_motore
    energia_annua_kWh = P_elettrica_kW * ore_annue
    costo_annuo = energia_annua_kWh * costo_kwh
    return {
        "P_idraulica [kW]": P_kW * rendimento_motore,  # potenza utile all'acqua
        "P_asse [kW]": P_kW,
        "P_elettrica [kW]": P_elettrica_kW,
        "Ore annue": ore_annue,
        "Energia annua [kWh]": energia_annua_kWh,
        "Costo unitario [EUR/kWh]": costo_kwh,
        "Costo annuo [EUR]": costo_annuo,
        "Rendimento motore": rendimento_motore,
    }


# ---------------------------------------------------------------------------
# Breakdown perdite per componente
# ---------------------------------------------------------------------------

def breakdown_perdite(res: dict) -> pd.DataFrame:
    """Tabella riassuntiva dei contributi alla prevalenza totale."""
    rows = [
        {"Componente": "Dislivello geodesico dZ",        "Valore [m]": round(res["dZ"], 4)},
        {"Componente": "Perdite distribuite asp. hf_s",  "Valore [m]": round(res["hf_s"], 4)},
        {"Componente": "Perdite concentrate asp. hK_s",  "Valore [m]": round(res["hK_s"], 4)},
        {"Componente": "Perdite distribuite mand. hf_d", "Valore [m]": round(res["hf_d"], 4)},
        {"Componente": "Perdite concentrate mand. hK_d", "Valore [m]": round(res["hK_d"], 4)},
        {"Componente": "Variazione cinetica dV^2/(2g)",  "Valore [m]": round(res["dVel"], 4)},
        {"Componente": "TDH totale",                      "Valore [m]": round(res["H"], 4)},
    ]
    df = pd.DataFrame(rows)
    
    # RISOLUZIONE ERRORE: 
    # Usiamo un costrutto if-else logico invece di un .where() di pandas su uno scalare
    if res["H"] != 0:
        df["% su TDH"] = (df["Valore [m]"] / res["H"] * 100).round(1)
    else:
        df["% su TDH"] = None
        
    return df


# ---------------------------------------------------------------------------
# Curva sistema TDH(Q)
# ---------------------------------------------------------------------------

def curva_tdh_vs_Q(suction: Linea, discharge: Linea, fluido: Fluido,
                   Q_ref: float, n_punti: int = 50) -> pd.DataFrame:
    """Curva del sistema TDH(Q) da 0.10*Q_ref a 2.0*Q_ref."""
    records = []
    for Q in np.linspace(0.10 * Q_ref, 2.0 * Q_ref, n_punti):
        r = tdh_pump(float(Q), suction, discharge, fluido)
        perdite_tot = r["hf_s"] + r["hf_d"] + r["hK_s"] + r["hK_d"]
        P = potenza_pompa(float(Q), r["H"], fluido.rho, 0.75)  # eta=0.75 per curva
        records.append({
            "Q [m3/s]":           round(float(Q), 5),
            "TDH [m]":            round(r["H"], 3),
            "V asp [m/s]":        round(r["V_s"], 3),
            "V mand [m/s]":       round(r["V_d"], 3),
            "Perdite tot [m]":    round(perdite_tot, 3),
            "dZ [m]":             round(r["dZ"], 3),
        })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Tabella passaggi di calcolo
# ---------------------------------------------------------------------------

def tabella_passaggi(Q: float, suction: Linea, discharge: Linea,
                     fluido: Fluido, eta: float, res: dict) -> pd.DataFrame:
    """Tabella con tutti i passaggi intermedi del calcolo pompe."""
    A_s = math.pi * suction.D ** 2 / 4.0
    A_d = math.pi * discharge.D ** 2 / 4.0
    eps_rel_s = suction.eps / suction.D
    eps_rel_d = discharge.eps / discharge.D
    reg_s = "laminare (f=64/Re)" if res["Re_s"] < 2300 else "turbolento (Swamee-Jain)"
    reg_d = "laminare (f=64/Re)" if res["Re_d"] < 2300 else "turbolento (Swamee-Jain)"
    P_kW = potenza_pompa(Q, res["H"], fluido.rho, eta) / 1000.0

    rows = [
        (1,  "Sezione asp. (area)",        "A_s",    "pi*Ds^2/4",
             f"{A_s:.6f}",       "m^2",  "Area della sezione circolare della condotta di aspirazione"),
        (2,  "Sezione mand. (area)",       "A_d",    "pi*Dd^2/4",
             f"{A_d:.6f}",       "m^2",  "Area della sezione circolare della condotta di mandata"),
        (3,  "Velocita' aspirazione",      "V_s",    "Q / A_s",
             f"{res['V_s']:.4f}","m/s",  "Velocita' media del fluido nella linea di aspirazione"),
        (4,  "Velocita' mandata",          "V_d",    "Q / A_d",
             f"{res['V_d']:.4f}","m/s",  "Velocita' media del fluido nella linea di mandata"),
        (5,  "Reynolds asp.",              "Re_s",   "V_s*Ds/nu",
             f"{res['Re_s']:.0f}","-",   "Numero di Reynolds: discrimina regime laminare (Re<2300) da turbolento"),
        (6,  "Reynolds mand.",             "Re_d",   "V_d*Dd/nu",
             f"{res['Re_d']:.0f}","-",   "Numero di Reynolds: discrimina regime laminare (Re<2300) da turbolento"),
        (7,  "Rugosita' relativa asp.",    "eps/Ds", "eps_s/Ds",
             f"{eps_rel_s:.7f}", "-",    "Rugosita' relativa: ingresso nell'equazione di Swamee-Jain"),
        (8,  "Rugosita' relativa mand.",   "eps/Dd", "eps_d/Dd",
             f"{eps_rel_d:.7f}", "-",    "Rugosita' relativa: ingresso nell'equazione di Swamee-Jain"),
        (9,  "Regime asp.",                "-",      "-",
             reg_s,              "-",    "Classificazione del regime di moto nella linea di aspirazione"),
        (10, "Regime mand.",               "-",      "-",
             reg_d,              "-",    "Classificazione del regime di moto nella linea di mandata"),
        (11, "Fattore attrito asp.",       "f_s",    "64/Re o Swamee-Jain",
             f"{res['f_s']:.6f}","-",   "Fattore di attrito di Darcy-Weisbach (approx. Colebrook-White)"),
        (12, "Fattore attrito mand.",      "f_d",    "64/Re o Swamee-Jain",
             f"{res['f_d']:.6f}","-",   "Fattore di attrito di Darcy-Weisbach (approx. Colebrook-White)"),
        (13, "Perdite distr. asp.",        "hf_s",   "f_s*(Ls/Ds)*V_s^2/(2g)",
             f"{res['hf_s']:.5f}","m",  "Perdite per attrito distribuite nell'intera lunghezza di aspirazione"),
        (14, "Perdite distr. mand.",       "hf_d",   "f_d*(Ld/Dd)*V_d^2/(2g)",
             f"{res['hf_d']:.5f}","m",  "Perdite per attrito distribuite nell'intera lunghezza di mandata"),
        (15, "Perdite conc. asp.",         "hK_s",   "Ks*V_s^2/(2g)",
             f"{res['hK_s']:.5f}","m",  "Perdite concentrate per valvole, curve, raccordi in aspirazione"),
        (16, "Perdite conc. mand.",        "hK_d",   "Kd*V_d^2/(2g)",
             f"{res['hK_d']:.5f}","m",  "Perdite concentrate per valvole, curve, raccordi in mandata"),
        (17, "Dislivello geodesico",       "dZ",     "z_d - z_s",
             f"{res['dZ']:.4f}", "m",   "Componente statica del TDH: dislivello tra punto di consegna e aspirazione"),
        (18, "Variazione cinetica",        "dVel",   "(V_d^2-V_s^2)/(2g)",
             f"{res['dVel']:.5f}","m",  "Differenza di energia cinetica tra mandata e aspirazione"),
        (19, "TDH totale",                 "H",      "dZ+hf_s+hf_d+hK_s+hK_d+dVel",
             f"{res['H']:.4f}",  "m",   "Prevalenza totale dinamica: parametro di selezione della pompa"),
        (20, "Potenza all'asse pompa",     "P",      "rho*g*Q*H/eta",
             f"{P_kW:.4f}",      "kW",  "Potenza meccanica richiesta alla pompa (divisa per il rendimento idraulico)"),
    ]
    return pd.DataFrame(rows, columns=["Passo", "Grandezza", "Simbolo",
                                       "Formula", "Valore", "Unita", "Descrizione"])


# ---------------------------------------------------------------------------
# Verifiche normative complete
# ---------------------------------------------------------------------------

def verifiche_pompa(Q: float, suction: Linea, discharge: Linea,
                    fluido: Fluido, eta: float, res: dict,
                    n_rpm: float = 1450.0,
                    z_serbatoio: Optional[float] = None,
                    npsh_richiesto: float = 3.0,
                    p_atm: float = P_ATM_STD) -> pd.DataFrame:
    """
    Tabella di verifiche ingegneristiche complete per il sistema pompa.
    Colonne: N., Verifica, Valore calcolato, Limite / soglia, Esito, Riferimento normativo.
    """
    P_kW = potenza_pompa(Q, res["H"], fluido.rho, eta) / 1000.0
    Ns = velocita_specifica_ns(n_rpm, Q, res["H"])

    z_suct = z_serbatoio if z_serbatoio is not None else suction.z
    hf_asp_tot = res["hf_s"] + res["hK_s"]
    npsh_disp = npsh_disponibile(z_suct, suction.z, hf_asp_tot, fluido.rho,
                                  p_atm, fluido.p_vap)
    margine_npsh = npsh_disp - npsh_richiesto

    ar = sovrapressione_joukowsky(res["V_d"], fluido.rho)
    dP_bar = ar["dP [bar]"]

    hK_tot = res["hK_s"] + res["hK_d"]
    hf_tot = res["hf_s"] + res["hf_d"]
    H = res["H"]

    def esito(cond: bool, fallback: str = "ATTENZIONE") -> str:
        return "OK" if cond else fallback

    checks = [
        (1, "Velocita asp. V_s <= 2.0 m/s",
         f"{res['V_s']:.3f} m/s", "<= 2.0 m/s",
         esito(res["V_s"] <= 2.0, "NON OK"), "Buona pratica progettuale / ISO 9906"),
        (2, "Velocita mand. V_d <= 3.5 m/s",
         f"{res['V_d']:.3f} m/s", "<= 3.5 m/s",
         esito(res["V_d"] <= 3.5, "NON OK"), "Buona pratica progettuale"),
        (3, "Regime asp. turbolento (Re_s > 4000)",
         f"{res['Re_s']:.0f}", "> 4000",
         esito(res["Re_s"] > 4000, "ATTENZIONE"), "Idraulica condotte"),
        (4, "Regime mand. turbolento (Re_d > 4000)",
         f"{res['Re_d']:.0f}", "> 4000",
         esito(res["Re_d"] > 4000, "ATTENZIONE"), "Idraulica condotte"),
        (5, "NPSH disponibile",
         f"{npsh_disp:.2f} m",
         f">= NPSH_r = {npsh_richiesto:.1f} m",
         "OK" if margine_npsh >= 0.5 else ("ATTENZIONE" if margine_npsh >= 0 else "NON OK"),
         "ISO 9906:2012, EN 12845:2015"),
        (6, "Margine NPSH_d - NPSH_r >= 0.5 m",
         f"{margine_npsh:.2f} m",
         ">= 0.5 m",
         esito(margine_npsh >= 0.5, "ATTENZIONE"), "ISO 9906:2012"),
        (7, "Rendimento pompa eta >= 0.65",
         f"{eta:.3f}", ">= 0.65",
         esito(eta >= 0.65, "ATTENZIONE"), "ErP 2019/1781 UE"),
        (8, "Perdite concentrate < 50% TDH",
         f"{hK_tot:.3f} m ({hK_tot/H*100:.1f}% se H>0)" if H > 0 else f"{hK_tot:.3f} m",
         "< 50% TDH",
         esito(H <= 0 or hK_tot < 0.5 * H, "ATTENZIONE"), "Buona pratica - ridurre accessori"),
        (9, "Perdite distribuite < 70% TDH",
         f"{hf_tot:.3f} m ({hf_tot/H*100:.1f}%)" if H > 0 else f"{hf_tot:.3f} m",
         "< 70% TDH",
         esito(H <= 0 or hf_tot < 0.7 * H, "ATTENZIONE"), "Buona pratica - aumentare DN"),
        (10, "TDH > 0 (configurazione valida)",
         f"{H:.3f} m", "> 0",
         esito(H > 0, "NON OK"), "Verifica configurazione impianto"),
        (11, "Velocita specifica Ns",
         f"{Ns:.0f}" if not math.isnan(Ns) else "n.d.",
         "300 -- 5000 (centrifuga)",
         "OK" if (not math.isnan(Ns) and 200 < Ns < 6000) else "ATTENZIONE",
         "KSB Pump Handbook / Stepanoff"),
        (12, "Sovrapressione Joukowsky (colpo ariete)",
         f"{dP_bar:.2f} bar",
         "Verif. vs PN condotta",
         "INFO", "EN 805:2000, Wylie & Streeter"),
    ]
    df = pd.DataFrame(checks, columns=["N.", "Verifica", "Valore calcolato",
                                        "Limite / soglia", "Esito",
                                        "Riferimento normativo"])
    return df


# ---------------------------------------------------------------------------
# Generazione report PDF
# ---------------------------------------------------------------------------

def _pdf_sezione(pdf, titolo: str) -> None:
    pdf.set_font("Helvetica", "B", 11)
    pdf.set_fill_color(30, 70, 140)
    pdf.set_text_color(255, 255, 255)
    pdf.cell(0, 7, titolo, ln=True, fill=True)
    pdf.set_text_color(0, 0, 0)
    pdf.ln(1)


def _pdf_riga_kv(pdf, chiave: str, valore: str) -> None:
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(72, 5, chiave + ":", border="B")
    pdf.set_font("Helvetica", "", 8)
    pdf.cell(0, 5, valore, border="B", ln=True)


def _pdf_tabella_gen(pdf, df: pd.DataFrame,
                     larghezze: Optional[dict] = None) -> None:
    cols = list(df.columns)
    lw = larghezze or {}
    default_w = max(10, int(190 / max(len(cols), 1)))
    row_h = 5
    pdf.set_font("Helvetica", "B", 7)
    pdf.set_fill_color(200, 220, 245)
    for col in cols:
        w = lw.get(col, default_w)
        pdf.cell(w, row_h + 1, str(col), border=1, align="C", fill=True)
    pdf.ln()
    pdf.set_font("Helvetica", "", 7)
    for _, row in df.iterrows():
        for col in cols:
            w = lw.get(col, default_w)
            val = row[col]
            txt = str(val) if val is not None else ""
            align = "C" if col in ("N.", "Passo", "Simbolo", "Valore", "Unita", "Esito",
                                   "Valore [m]", "% su TDH") else "L"
            max_c = max(4, int(w / 1.85))
            if len(txt) > max_c:
                txt = txt[: max_c - 2] + ".."
            pdf.cell(w, row_h, txt, border=1, align=align)
        pdf.ln()


def genera_pdf(Q: float, suction: Linea, discharge: Linea,
               fluido: Fluido, eta: float, res: dict, note: List[str],
               n_rpm: float = 1450.0, z_serbatoio: Optional[float] = None,
               npsh_richiesto: float = 3.0, ore_annue: float = 2000.0,
               costo_kwh: float = 0.15, eta_motore: float = 0.92) -> bytes:
    """Genera un report PDF completo e restituisce i bytes."""
    from fpdf import FPDF

    df_pass = tabella_passaggi(Q, suction, discharge, fluido, eta, res)
    df_bkd = breakdown_perdite(res)
    df_ver = verifiche_pompa(Q, suction, discharge, fluido, eta, res,
                              n_rpm, z_serbatoio, npsh_richiesto)
    P_kW = potenza_pompa(Q, res["H"], fluido.rho, eta) / 1000.0
    Ns = velocita_specifica_ns(n_rpm, Q, res["H"])
    en = costo_energetico_annuo(P_kW, ore_annue, costo_kwh, eta_motore)
    ar = sovrapressione_joukowsky(res["V_d"], fluido.rho)

    lw_pass = {"Passo": 9, "Grandezza": 35, "Simbolo": 18,
               "Formula": 45, "Valore": 24, "Unita": 14, "Descrizione": 45}
    lw_bkd = {"Componente": 100, "Valore [m]": 30, "% su TDH": 25}
    lw_ver = {"N.": 8, "Verifica": 65, "Valore calcolato": 28,
              "Limite / soglia": 28, "Esito": 18, "Riferimento normativo": 43}

    pdf = FPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    pdf.set_font("Helvetica", "B", 15)
    pdf.set_fill_color(20, 50, 110)
    pdf.set_text_color(255, 255, 255)
    pdf.cell(0, 12, "Report - Prevalenza Pompa (TDH) e Verifiche Sistema",
             ln=True, align="C", fill=True)
    pdf.set_text_color(0, 0, 0)
    pdf.set_font("Helvetica", "", 8)
    pdf.cell(0, 6, f"Generato il {datetime.date.today().strftime('%d/%m/%Y')}", ln=True, align="C")
    pdf.ln(4)

    _pdf_sezione(pdf, "1. Parametri di input")
    _pdf_riga_kv(pdf, "Portata Q", f"{Q:.5f} m3/s")
    _pdf_riga_kv(pdf, "Densita fluido rho", f"{fluido.rho:.2f} kg/m3")
    _pdf_riga_kv(pdf, "Viscosita nu", f"{fluido.nu:.4e} m2/s")
    _pdf_riga_kv(pdf, "Pressione vapore p_vap", f"{fluido.p_vap:.0f} Pa")
    _pdf_riga_kv(pdf, "Rendimento pompa eta", f"{eta:.3f}")
    _pdf_riga_kv(pdf, "Velocita di rotazione n", f"{n_rpm:.0f} giri/min")
    _pdf_riga_kv(pdf, "Asp.: L/D/eps/K/z",
                 f"{suction.L:.1f}m / {suction.D:.3f}m / {suction.eps:.5f}m / {suction.K_tot:.1f} / {suction.z:.2f}m")
    _pdf_riga_kv(pdf, "Mand.: L/D/eps/K/z",
                 f"{discharge.L:.1f}m / {discharge.D:.3f}m / {discharge.eps:.5f}m / {discharge.K_tot:.1f} / {discharge.z:.2f}m")
    pdf.ln(4)

    _pdf_sezione(pdf, "2. Risultati principali")
    _pdf_riga_kv(pdf, "TDH totale", f"{res['H']:.4f} m")
    _pdf_riga_kv(pdf, "Potenza all'asse pompa P", f"{P_kW:.3f} kW")
    _pdf_riga_kv(pdf, "V_s (aspirazione)", f"{res['V_s']:.4f} m/s  |  Re_s = {res['Re_s']:.0f}")
    _pdf_riga_kv(pdf, "V_d (mandata)", f"{res['V_d']:.4f} m/s  |  Re_d = {res['Re_d']:.0f}")
    _pdf_riga_kv(pdf, "f_s (Darcy-Weisbach asp.)", f"{res['f_s']:.6f}")
    _pdf_riga_kv(pdf, "f_d (Darcy-Weisbach mand.)", f"{res['f_d']:.6f}")
    _pdf_riga_kv(pdf, "Velocita specifica Ns", f"{Ns:.0f}  ->  {tipo_pompa_da_ns(Ns)}" if not math.isnan(Ns) else "n.d.")
    _pdf_riga_kv(pdf, "NPSH disponibile", f"{npsh_disponibile(z_serbatoio or suction.z, suction.z, res['hf_s']+res['hK_s'], fluido.rho, P_ATM_STD, fluido.p_vap):.3f} m")
    _pdf_riga_kv(pdf, "Sovrapressione Joukowsky", f"{ar['dP [bar]']:.2f} bar  ({ar['dP [m c.a.]']:.1f} m c.a.)")
    _pdf_riga_kv(pdf, "Costo energetico annuo", f"{en['Costo annuo [EUR]']:.0f} EUR/anno  ({en['Energia annua [kWh]']:.0f} kWh/anno)")
    pdf.ln(4)

    _pdf_sezione(pdf, "3. Dettaglio componenti TDH")
    _pdf_tabella_gen(pdf, df_bkd, lw_bkd)
    pdf.ln(4)

    _pdf_sezione(pdf, "4. Passaggi di calcolo")
    _pdf_tabella_gen(pdf, df_pass, lw_pass)

    pdf.add_page()
    _pdf_sezione(pdf, "5. Verifiche normative")
    _pdf_tabella_gen(pdf, df_ver, lw_ver)
    pdf.ln(4)

    _pdf_sezione(pdf, "6. Note tecniche")
    pdf.set_font("Helvetica", "", 8)
    for item in note:
        pdf.multi_cell(0, 5, "- " + item)
        pdf.ln(1)

    return pdf.output()


# ---------------------------------------------------------------------------
# Commenti progettuali automatici
# ---------------------------------------------------------------------------

def commenti_progettuali(Q: float, suction: Linea, discharge: Linea,
                          fluido: Fluido, res: dict, eta: float) -> List[str]:
    note: List[str] = []
    Re_s, Re_d = res["Re_s"], res["Re_d"]
    V_s, V_d = res["V_s"], res["V_d"]
    H = res["H"]
    hK_tot = res["hK_s"] + res["hK_d"]
    hf_tot = res["hf_s"] + res["hf_d"]

    if Re_s < 4000:
        note.append(
            f"Linea asp. in regime laminare/transitorio (Re = {Re_s:.0f}): f = 64/Re. "
            "Aumentare il diametro per passare a regime turbolento."
        )
    if Re_d < 4000:
        note.append(
            f"Linea mand. in regime laminare/transitorio (Re = {Re_d:.0f}): f = 64/Re."
        )
    if V_s > 2.0:
        note.append(
            f"Velocita asp. elevata ({V_s:.2f} m/s > 2.0 m/s): rischio cavitazione. "
            "Verificare NPSH disponibile e aumentare il diametro di aspirazione."
        )
    if V_d > 3.5:
        note.append(
            f"Velocita mand. elevata ({V_d:.2f} m/s > 3.5 m/s): verificare abrasione "
            "e colpo d'ariete alla chiusura valvola."
        )
    if H > 0 and hK_tot > 0.5 * H:
        note.append(
            f"Perdite concentrate ({hK_tot:.2f} m) > 50% del TDH: ridurre accessori, "
            "aumentare i diametri, ridistribuire i carichi locali."
        )
    if H > 0 and hf_tot > 0.7 * H:
        note.append(
            f"Perdite distribuite ({hf_tot:.2f} m) > 70% del TDH: valutare l'aumento "
            "del diametro delle tubazioni per ridurre le perdite."
        )
    if eta < 0.65:
        note.append(
            f"Rendimento basso (eta = {eta:.2f} < 0.65): verificare il punto di lavoro "
            "sulla curva caratteristica e l'invecchiamento della macchina."
        )
    if H < 0:
        note.append(
            "TDH negativo: verificare le quote z_s, z_d, i segni delle perdite e la "
            "configurazione dell'impianto."
        )
    if not note:
        note.append(
            "I parametri rientrano in un intervallo operativo tipico. "
            "Completare la verifica con la curva caratteristica della pompa, l'analisi NPSH "
            "e la stima del colpo d'ariete."
        )
    return note
