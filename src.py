# -*- coding: utf-8 -*-
"""
CONDOTTE - Calcolo idraulico sezione circolare (moto uniforme, Manning)
Versione professionale:
- Geometria e idraulica completa della sezione circolare
- Libreria materiali con n di Manning e velocita massima
- Diametri DN standard (EN 476 / ISO)
- Confronto automatico diametri adiacenti
- Numero di Froude, classificazione del moto, energia specifica
- Tabella verifiche normative (beta, V_min, V_max, Fr, margine)
- Curva idraulica completa Q(beta), V(beta), A(beta)
- Tabella passaggi di calcolo
- Generazione report PDF
"""
from __future__ import annotations

import datetime
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

G = 9.81  # accelerazione di gravita [m/s²]


# ---------------------------------------------------------------------------
# Libreria materiali (n Manning e V_max [m/s])
# ---------------------------------------------------------------------------

MATERIALI_MANNING: Dict[str, dict] = {
    "PVC / PE liscio (HDPE)": {
        "n": 0.010, "V_max": 5.0,
        "note": "Condotte plastiche nuove. Massima efficienza idraulica.",
    },
    "Gres ceramico smaltato": {
        "n": 0.011, "V_max": 4.5,
        "note": "Fognatura mista o nera; eccellente resistenza all'abrasione.",
    },
    "Calcestruzzo liscio (prefabbricato)": {
        "n": 0.012, "V_max": 6.0,
        "note": "Tubi in calcestruzzo prefabbricato, buone condizioni.",
    },
    "Ghisa sferoidale (duttile)": {
        "n": 0.013, "V_max": 4.0,
        "note": "Condotte in pressione e gravita. Lunga vita utile.",
    },
    "Acciaio rivestito (epossidico)": {
        "n": 0.012, "V_max": 5.0,
        "note": "Acciaio con rivestimento interno in resina epossidica.",
    },
    "Calcestruzzo armato (scatolare)": {
        "n": 0.015, "V_max": 4.0,
        "note": "Culvert, scatolari, canali in calcestruzzo armato.",
    },
    "HDPE corrugato (doppia parete)": {
        "n": 0.018, "V_max": 4.0,
        "note": "HDPE a parete corrugata esterna / liscia interna.",
    },
    "Acciaio non rivestito / vecchio": {
        "n": 0.020, "V_max": 3.0,
        "note": "Condotte metalliche degradate, rugosita elevata.",
    },
    "Muratura in mattoni": {
        "n": 0.020, "V_max": 2.5,
        "note": "Fognature storiche. Rugosita variabile. Da ispezionare.",
    },
    "Fibrocemento (storico)": {
        "n": 0.011, "V_max": 4.0,
        "note": "ATTENZIONE: materiale contenente amianto. Non usare in nuove opere.",
    },
}

# ---------------------------------------------------------------------------
# Diametri nominali DN standard (EN 476 / ISO) [m]
# ---------------------------------------------------------------------------

DN_STANDARD_M: List[float] = [
    0.100, 0.125, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500,
    0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.400, 1.600, 1.800,
    2.000, 2.200, 2.400, 2.600, 2.800, 3.000,
]

# Velocita minima autocircolante per tipo di uso [m/s]
V_MIN_AUTOCIRCOLANTE: Dict[str, float] = {
    "Fognatura nera (acque reflue)": 0.60,
    "Fognatura mista": 0.70,
    "Acque meteoriche / pluviali": 0.50,
    "Irrigazione / acquedotto": 0.30,
    "Generico": 0.40,
}


# ---------------------------------------------------------------------------
# Dataclass input
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DatiCondotta:
    D: float    # diametro [m]
    n: float    # coeff. di Manning [-]
    S: float    # pendenza idraulica [-]
    Q: float    # portata di progetto [m³/s]


# ---------------------------------------------------------------------------
# Validazione
# ---------------------------------------------------------------------------

def valida_dati(dati: DatiCondotta) -> List[str]:
    errori: List[str] = []
    if dati.D <= 0:
        errori.append("Il diametro D deve essere maggiore di 0 m.")
    if dati.n <= 0 or dati.n > 0.05:
        errori.append("Il coefficiente di Manning n deve essere compreso tra 0.005 e 0.05.")
    if dati.S <= 0:
        errori.append("La pendenza idraulica S deve essere positiva.")
    if dati.Q <= 0:
        errori.append("La portata Q deve essere positiva.")
    return errori


# ---------------------------------------------------------------------------
# Geometria sezione circolare
# ---------------------------------------------------------------------------

def theta_from_y(y: float, D: float) -> float:
    """Angolo al centro theta [rad] corrispondente al tirante y."""
    if y <= 0:
        return 0.0
    if y >= D:
        return math.pi
    cos_theta = max(-1.0, min(1.0, 1.0 - 2.0 * y / D))
    return math.acos(cos_theta)


def area_perimetro(y: float, D: float) -> Tuple[float, float]:
    """Area bagnata A [m²] e perimetro bagnato P [m] per sezione circolare."""
    R = D / 2.0
    theta = theta_from_y(y, D)
    A = 0.5 * R * R * (2.0 * theta - math.sin(2.0 * theta))
    P = 2.0 * R * theta
    return A, P


def larghezza_superficie_libera(y: float, D: float) -> float:
    """Larghezza della superficie libera T [m] per sezione circolare."""
    if y <= 0 or y >= D:
        return 0.0
    R = D / 2.0
    return 2.0 * math.sqrt(max(0.0, R ** 2 - (R - y) ** 2))


# ---------------------------------------------------------------------------
# Formula di Manning e Froude
# ---------------------------------------------------------------------------

def manning_Q(A: float, P: float, n: float, S: float) -> float:
    if P <= 0 or n <= 0 or S <= 0:
        return 0.0
    Rh = A / P
    return (1.0 / n) * A * (Rh ** (2.0 / 3.0)) * (S ** 0.5)


def numero_froude(V: float, A: float, T: float) -> float:
    """Numero di Froude per sezione parzialmente piena: Fr = V / sqrt(g * A/T)."""
    if T <= 0 or A <= 0:
        return float("nan")
    D_idraulico = A / T
    return V / math.sqrt(G * D_idraulico)


def classifica_moto(Fr: float) -> str:
    """Classificazione del regime del moto."""
    if math.isnan(Fr):
        return "n.d."
    if abs(Fr - 1.0) <= 0.05:
        return "Critico (Fr ~ 1)"
    return "Subcritico (Fr < 1)" if Fr < 1.0 else "Supercritico (Fr > 1)"


def energia_specifica(V: float, y: float) -> float:
    """Energia specifica E = y + V²/(2g) [m]."""
    return y + V ** 2 / (2.0 * G)


# ---------------------------------------------------------------------------
# Risoluzione numerica y(Q) — bisezione
# ---------------------------------------------------------------------------

def risolvi_y_per_Q(D: float, n: float, S: float, Q_target: float,
                    tol: float = 1e-8, max_iter: int = 300) -> float:
    """Risolve y dato Q_target per una sezione circolare con Manning (bisezione)."""
    y_max_search = 0.938 * D
    a, b = 1e-9, y_max_search
    Q_a = manning_Q(*area_perimetro(a, D), n, S)
    Q_b = manning_Q(*area_perimetro(b, D), n, S)
    if Q_target <= Q_a:
        return a
    if Q_target >= Q_b:
        return b
    ya, yb = a, b
    for _ in range(max_iter):
        ym = 0.5 * (ya + yb)
        Qm = manning_Q(*area_perimetro(ym, D), n, S)
        if abs(Qm - Q_target) < tol:
            return ym
        if Qm < Q_target:
            ya = ym
        else:
            yb = ym
    return 0.5 * (ya + yb)


def portata_massima(D: float, n: float, S: float, n_punti: int = 500) -> float:
    """Portata massima nella sezione circolare (beta ~ 0.938)."""
    Q_max = 0.0
    for i in range(1, n_punti):
        y = i / n_punti * D
        A, P = area_perimetro(y, D)
        Q = manning_Q(A, P, n, S)
        if Q > Q_max:
            Q_max = Q
    return Q_max


def portata_piena_sezione(D: float, n: float, S: float) -> float:
    """Portata a sezione piena (beta = 1.0)."""
    A, P = area_perimetro(D, D)
    return manning_Q(A, P, n, S)


# ---------------------------------------------------------------------------
# Stato della condotta al punto di progetto
# ---------------------------------------------------------------------------

def stato_condotta(dati: DatiCondotta) -> dict:
    y = risolvi_y_per_Q(dati.D, dati.n, dati.S, dati.Q)
    A, P = area_perimetro(y, dati.D)
    Rh = A / P if P > 0 else float("nan")
    V = dati.Q / A if A > 0 else float("nan")
    beta = y / dati.D
    T = larghezza_superficie_libera(y, dati.D)
    Fr = numero_froude(V, A, T)
    E = energia_specifica(V, y)
    Q_max = portata_massima(dati.D, dati.n, dati.S)
    Q_full = portata_piena_sezione(dati.D, dati.n, dati.S)
    return dict(y=y, beta=beta, A=A, P=P, Rh=Rh, V=V, T=T,
                Fr=Fr, E=E, Q_max=Q_max, Q_full=Q_full)


# ---------------------------------------------------------------------------
# Curva idraulica completa
# ---------------------------------------------------------------------------

def curva_riempimento(D: float, n: float, S: float, n_punti: int = 80) -> pd.DataFrame:
    """Curva idraulica completa per la sezione circolare."""
    records = []
    for i in range(1, n_punti + 1):
        beta = i / n_punti
        y = beta * D
        A, P = area_perimetro(y, D)
        Q = manning_Q(A, P, n, S)
        V = Q / A if A > 0 else 0.0
        Rh = A / P if P > 0 else 0.0
        T = larghezza_superficie_libera(y, D)
        Fr = numero_froude(V, A, T)
        E = energia_specifica(V, y)
        records.append({
            "y/D [-]": round(beta, 4),
            "y [m]": round(y, 5),
            "Q [m3/s]": round(Q, 6),
            "V [m/s]": round(V, 4),
            "A [m2]": round(A, 6),
            "P [m]": round(P, 5),
            "Rh [m]": round(Rh, 5),
            "T [m]": round(T, 5),
            "Fr [-]": round(Fr, 4) if not math.isnan(Fr) else None,
            "E [m]": round(E, 5),
        })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Confronto diametri DN standard
# ---------------------------------------------------------------------------

def confronto_diametri(n: float, S: float, Q: float,
                       D_ref: float,
                       materiale: str = "Calcestruzzo liscio (prefabbricato)",
                       uso: str = "Fognatura nera (acque reflue)") -> pd.DataFrame:
    """
    Confronta i diametri DN standard nell'intervallo 0.5*D_ref -- 2.5*D_ref.
    Per ciascun diametro calcola beta, V, Q_max, Fr e verifica le condizioni.
    """
    D_list = [d for d in DN_STANDARD_M if 0.5 * D_ref <= d <= 2.5 * D_ref]
    V_max = MATERIALI_MANNING.get(materiale, {}).get("V_max", 5.0)
    V_min = V_MIN_AUTOCIRCOLANTE.get(uso, 0.60)
    records = []
    for D in D_list:
        try:
            dati_d = DatiCondotta(D=D, n=n, S=S, Q=Q)
            if valida_dati(dati_d):
                continue
            res_d = stato_condotta(dati_d)
            beta = res_d["beta"]
            V = res_d["V"]
            Q_max = res_d["Q_max"]
            Fr = res_d["Fr"]

            ok_beta = beta <= 0.80
            ok_V = V_min <= V <= V_max
            ok_Q = Q <= Q_max
            stato = "OK" if (ok_beta and ok_V and ok_Q) else (
                "ATTENZIONE" if (beta <= 0.94 and ok_Q) else "NON OK")
            sel = "<<< DN progetto" if abs(D - D_ref) < 0.001 else ""

            records.append({
                "DN [mm]": int(round(D * 1000)),
                "beta [-]": round(beta, 3),
                "V [m/s]": round(V, 3),
                "Fr [-]": round(Fr, 3) if not math.isnan(Fr) else None,
                "Q_max [m3/s]": round(Q_max, 4),
                "Q/Q_max [-]": round(Q / Q_max, 3),
                "OK beta": "SI" if ok_beta else "NO",
                "OK V": "SI" if ok_V else "NO",
                "Stato globale": stato,
                "Nota": sel,
            })
        except Exception:
            continue
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Verifiche normative (tabella semaforo)
# ---------------------------------------------------------------------------

def verifiche_idrauliche(dati: DatiCondotta, res: dict,
                         materiale: str = "Calcestruzzo liscio (prefabbricato)",
                         uso: str = "Fognatura nera (acque reflue)") -> pd.DataFrame:
    """
    Tabella di verifica normativa completa della condotta.
    Restituisce colonne: Verifica, Valore, Limite, Esito, Riferimento.
    Esito puo' essere: OK | ATTENZIONE | NON OK | INFO.
    """
    mat = MATERIALI_MANNING.get(materiale, {})
    V_max_mat = mat.get("V_max", 5.0)
    V_min = V_MIN_AUTOCIRCOLANTE.get(uso, 0.60)

    beta = res["beta"]
    V = res["V"]
    Q = dati.Q
    Q_max = res["Q_max"]
    Q_full = res["Q_full"]
    Fr = res["Fr"]
    E = res["E"]
    Rh = res["Rh"]
    y = res["y"]
    margine = Q_full / Q if Q > 0 else float("nan")
    regime = classifica_moto(Fr)

    Fr_str = f"{Fr:.4f}" if not math.isnan(Fr) else "n.d."
    esito_Fr = (("OK" if Fr < 0.95 else "ATTENZIONE") if Fr < 1.0 else
                ("ATTENZIONE" if Fr < 1.05 else "NON OK")) if not math.isnan(Fr) else "INFO"

    checks = [
        {
            "N.": 1,
            "Verifica": "Grado di riempimento beta <= 0.80",
            "Valore calcolato": f"{beta:.4f}",
            "Limite / soglia": "<= 0.80",
            "Esito": "OK" if beta <= 0.80 else ("ATTENZIONE" if beta <= 0.94 else "NON OK"),
            "Riferimento normativo": "DIN EN 1671:2011, EN 476:2011",
        },
        {
            "N.": 2,
            "Verifica": f"Velocita min. autocircolante V >= {V_min:.2f} m/s",
            "Valore calcolato": f"{V:.4f} m/s",
            "Limite / soglia": f">= {V_min:.2f} m/s",
            "Esito": "OK" if V >= V_min else "NON OK",
            "Riferimento normativo": "DIN EN 1671:2011 - ASCE MOP 36",
        },
        {
            "N.": 3,
            "Verifica": f"Velocita max. materiale V <= {V_max_mat:.1f} m/s",
            "Valore calcolato": f"{V:.4f} m/s",
            "Limite / soglia": f"<= {V_max_mat:.1f} m/s",
            "Esito": "OK" if V <= V_max_mat else "NON OK",
            "Riferimento normativo": f"Scheda tecnica materiale ({materiale})",
        },
        {
            "N.": 4,
            "Verifica": "Q progetto <= Q_max (beta = 0.94)",
            "Valore calcolato": f"{Q:.5f} m3/s",
            "Limite / soglia": f"<= {Q_max:.5f} m3/s",
            "Esito": "OK" if Q <= Q_max else "NON OK",
            "Riferimento normativo": "Idraulica sezione circolare",
        },
        {
            "N.": 5,
            "Verifica": "Regime di moto subcritico (Fr < 1.0)",
            "Valore calcolato": Fr_str,
            "Limite / soglia": "Fr < 1.0",
            "Esito": esito_Fr,
            "Riferimento normativo": "Idraulica canali in pressione",
        },
        {
            "N.": 6,
            "Verifica": "Margine capacita' Q_full / Q >= 1.25",
            "Valore calcolato": f"{margine:.3f}",
            "Limite / soglia": ">= 1.25",
            "Esito": "OK" if margine >= 1.25 else "ATTENZIONE",
            "Riferimento normativo": "Linee guida ISPRA / buona pratica",
        },
        {
            "N.": 7,
            "Verifica": "Rapporto Q/Q_max [-]",
            "Valore calcolato": f"{Q / Q_max:.4f}",
            "Limite / soglia": "<= 1.0",
            "Esito": "OK" if Q <= Q_max else "NON OK",
            "Riferimento normativo": "Verifica di portata",
        },
        {
            "N.": 8,
            "Verifica": "Energia specifica E = y + V2/(2g)",
            "Valore calcolato": f"{E:.5f} m",
            "Limite / soglia": "INFO",
            "Esito": "INFO",
            "Riferimento normativo": "Idraulica di base",
        },
        {
            "N.": 9,
            "Verifica": "Regime di moto classificato",
            "Valore calcolato": regime,
            "Limite / soglia": "Subcritico preferito",
            "Esito": "INFO",
            "Riferimento normativo": "Classificazione idraulica",
        },
        {
            "N.": 10,
            "Verifica": "Pendenza minima autocircolante (V = V_min)",
            "Valore calcolato": f"S progetto = {dati.S:.6f}",
            "Limite / soglia": f"S >= S_min per V={V_min:.2f} m/s",
            "Esito": "OK" if V >= V_min else "ATTENZIONE",
            "Riferimento normativo": "DIN EN 1671:2011",
        },
    ]
    return pd.DataFrame(checks)


# ---------------------------------------------------------------------------
# Tabella passaggi di calcolo
# ---------------------------------------------------------------------------

def tabella_passaggi(dati: DatiCondotta, res: dict) -> pd.DataFrame:
    """Tabella con tutti i passaggi intermedi del calcolo idraulico."""
    R = dati.D / 2.0
    y = res["y"]
    beta = res["beta"]
    theta = theta_from_y(y, dati.D)
    A = res["A"]
    P = res["P"]
    Rh = res["Rh"]
    V = res["V"]
    T = res["T"]
    Fr = res["Fr"]
    E = res["E"]
    Q_calc = manning_Q(A, P, dati.n, dati.S)
    Q_max = res["Q_max"]
    Q_full = res["Q_full"]

    D_idr_val = A / T if T > 0 else float("nan")
    Q_Qmax  = dati.Q / Q_max  if Q_max  > 0 else float("nan")
    Q_Qfull = dati.Q / Q_full if Q_full > 0 else float("nan")

    rows = [
        (1,  "Raggio interno",             "R",       "D / 2",
             f"{R:.5f}",               "m",
             "Raggio geometrico della sezione circolare"),
        (2,  "Tirante idraulico",          "y",       "Manning^-1(Q)",
             f"{y:.5f}",               "m",
             "Altezza del pelo libero, trovata per bisezione sull'equazione di Manning"),
        (3,  "Grado di riempimento",       "beta",    "y / D",
             f"{beta:.5f}",            "-",
             "Frazione del diametro occupata dall'acqua (0=vuota, 1=piena)"),
        (4,  "Angolo al centro",           "theta",   "arccos(1 - 2*y/D)",
             f"{theta:.5f}",           "rad",
             "Semi-angolo al centro del settore circolare bagnato"),
        (5,  "Area bagnata",               "A",       "R^2*(theta-sin(2*theta)/2)",
             f"{A:.6f}",               "m^2",
             "Sezione trasversale occupata dall'acqua: ingresso in Manning"),
        (6,  "Perimetro bagnato",          "P",       "2 * R * theta",
             f"{P:.6f}",               "m",
             "Lunghezza del contorno a contatto con l'acqua: ingresso in Manning"),
        (7,  "Raggio idraulico",           "Rh",      "A / P",
             f"{Rh:.6f}",              "m",
             "Parametro chiave di Manning: esprime l'efficienza della sezione"),
        (8,  "Larghezza superficie libera","T",       "2*sqrt(R^2-(R-y)^2)",
             f"{T:.6f}",               "m",
             "Larghezza del pelo libero: necessaria per il calcolo del numero di Froude"),
        (9,  "Portata Manning (verifica)", "Q_calc",  "(1/n)*A*Rh^(2/3)*S^(1/2)",
             f"{Q_calc:.6f}",          "m^3/s",
             "Portata ricalcolata dai parametri di sezione: deve coincidere con Q input"),
        (10, "Velocita' media",            "V",       "Q / A",
             f"{V:.5f}",               "m/s",
             "Velocita' media nella sezione bagnata: verificare vs V_min e V_max"),
        (11, "Diametro idraulico",         "D_idr",   "A / T",
             f"{D_idr_val:.5f}" if not math.isnan(D_idr_val) else "n.d.", "m",
             "Profondita' idraulica media della sezione parzialmente piena"),
        (12, "Numero di Froude",           "Fr",      "V / sqrt(g * A/T)",
             f"{Fr:.5f}" if not math.isnan(Fr) else "n.d.", "-",
             "Fr<1: moto subcritico (corretto); Fr=1: critico; Fr>1: supercritico"),
        (13, "Energia specifica",          "E",       "y + V^2 / (2g)",
             f"{E:.5f}",               "m",
             "Energia totale per unita' di peso (utile per analisi del moto)"),
        (14, "Portata piena sezione",      "Q_full",  "Manning(A_piena, P_piena)",
             f"{Q_full:.6f}",          "m^3/s",
             "Portata a sezione completamente piena (beta=1): limite teorico"),
        (15, "Portata massima (beta~0.94)","Q_max",   "max Q(beta) in [0,1]",
             f"{Q_max:.6f}",           "m^3/s",
             "Portata massima effettiva: il picco della curva Q-beta si ha a beta~0.94"),
        (16, "Rapporto Q / Q_max",         "Q/Q_max", "Q / Q_max",
             f"{Q_Qmax:.5f}",          "-",
             "Indice di utilizzo della condotta: deve essere <= 1.0"),
        (17, "Rapporto Q / Q_full",        "Q/Q_full","Q / Q_full",
             f"{Q_Qfull:.5f}",         "-",
             "Riempimento normalizzato sulla portata a sezione piena"),
    ]
    return pd.DataFrame(rows, columns=["Passo", "Grandezza", "Simbolo",
                                       "Formula", "Valore", "Unita", "Descrizione"])


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
                     larghezze: Optional[Dict[str, int]] = None) -> None:
    """Tabella generica da DataFrame nel PDF."""
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
            align = "C" if col in ("N.", "Passo", "Simbolo", "Valore", "Unita",
                                   "Esito", "OK beta", "OK V", "Stato globale",
                                   "DN [mm]") else "L"
            max_c = max(4, int(w / 1.85))
            if len(txt) > max_c:
                txt = txt[: max_c - 2] + ".."
            pdf.cell(w, row_h, txt, border=1, align=align)
        pdf.ln()


def genera_pdf(dati: DatiCondotta, res: dict, note: List[str],
               materiale: str = "", uso: str = "") -> bytes:
    """Genera un report PDF completo e restituisce i bytes."""
    from fpdf import FPDF

    df_pass = tabella_passaggi(dati, res)
    df_ver = verifiche_idrauliche(dati, res, materiale or list(MATERIALI_MANNING.keys())[2], uso)
    df_dn = confronto_diametri(dati.n, dati.S, dati.Q, dati.D, materiale or list(MATERIALI_MANNING.keys())[2], uso)

    lw_pass = {"Passo": 8, "Grandezza": 36, "Simbolo": 16,
               "Formula": 48, "Valore": 24, "Unita": 14, "Descrizione": 44}
    lw_ver = {"N.": 8, "Verifica": 68, "Valore calcolato": 30,
              "Limite / soglia": 28, "Esito": 18, "Riferimento normativo": 38}
    lw_dn = {"DN [mm]": 15, "beta [-]": 15, "V [m/s]": 14, "Fr [-]": 12,
             "Q_max [m3/s]": 22, "Q/Q_max [-]": 20, "OK beta": 15,
             "OK V": 12, "Stato globale": 25, "Nota": 30}

    pdf = FPDF(orientation="P", unit="mm", format="A4")
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    # Titolo
    pdf.set_font("Helvetica", "B", 16)
    pdf.set_fill_color(20, 50, 110)
    pdf.set_text_color(255, 255, 255)
    pdf.cell(0, 12, "Report - Condotta Circolare (Moto Uniforme - Manning)",
             ln=True, align="C", fill=True)
    pdf.set_text_color(0, 0, 0)
    pdf.set_font("Helvetica", "", 8)
    pdf.cell(0, 6, f"Generato il {datetime.date.today().strftime('%d/%m/%Y')}  |  "
             f"Materiale: {materiale}  |  Uso: {uso}", ln=True, align="C")
    pdf.ln(4)

    # Input
    _pdf_sezione(pdf, "1. Parametri di input")
    _pdf_riga_kv(pdf, "Diametro D", f"{dati.D:.4f} m  (DN {int(round(dati.D*1000))} mm)")
    _pdf_riga_kv(pdf, "Coefficiente Manning n", f"{dati.n:.4f}  (materiale: {materiale})")
    _pdf_riga_kv(pdf, "Pendenza idraulica S", f"{dati.S:.6f}")
    _pdf_riga_kv(pdf, "Portata di progetto Q", f"{dati.Q:.5f} m3/s")
    _pdf_riga_kv(pdf, "Uso previsto", uso)
    pdf.ln(4)

    # Risultati
    _pdf_sezione(pdf, "2. Risultati principali al punto di progetto")
    _pdf_riga_kv(pdf, "Tirante y", f"{res['y']:.5f} m")
    _pdf_riga_kv(pdf, "Grado riempimento beta = y/D", f"{res['beta']:.5f}")
    _pdf_riga_kv(pdf, "Area bagnata A", f"{res['A']:.6f} m2")
    _pdf_riga_kv(pdf, "Perimetro bagnato P", f"{res['P']:.6f} m")
    _pdf_riga_kv(pdf, "Raggio idraulico Rh", f"{res['Rh']:.6f} m")
    _pdf_riga_kv(pdf, "Velocita media V", f"{res['V']:.4f} m/s")
    _pdf_riga_kv(pdf, "Numero di Froude Fr", f"{res['Fr']:.4f}  ({classifica_moto(res['Fr'])})"
                 if not math.isnan(res["Fr"]) else "n.d.")
    _pdf_riga_kv(pdf, "Energia specifica E", f"{res['E']:.5f} m")
    _pdf_riga_kv(pdf, "Q piena sezione Q_full", f"{res['Q_full']:.5f} m3/s")
    _pdf_riga_kv(pdf, "Q massima Q_max (beta~0.94)", f"{res['Q_max']:.5f} m3/s")
    _pdf_riga_kv(pdf, "Q / Q_max", f"{dati.Q / res['Q_max']:.4f}")
    pdf.ln(4)

    # Passaggi
    _pdf_sezione(pdf, "3. Passaggi di calcolo (passo per passo)")
    _pdf_tabella_gen(pdf, df_pass, lw_pass)
    pdf.ln(4)

    # Verifiche
    pdf.add_page()
    _pdf_sezione(pdf, "4. Verifiche normative (tabella semaforo)")
    _pdf_tabella_gen(pdf, df_ver, lw_ver)
    pdf.ln(4)

    # Confronto DN
    _pdf_sezione(pdf, "5. Confronto diametri DN standard adiacenti")
    _pdf_tabella_gen(pdf, df_dn, lw_dn)
    pdf.ln(4)

    # Note
    _pdf_sezione(pdf, "6. Note tecniche e commenti progettuali")
    pdf.set_font("Helvetica", "", 8)
    for item in note:
        pdf.multi_cell(0, 5, "- " + item)
        pdf.ln(1)

    return bytes(pdf.output())


# ---------------------------------------------------------------------------
# Commenti progettuali automatici
# ---------------------------------------------------------------------------

def commenti_progettuali(dati: DatiCondotta, res: dict,
                         materiale: str = "", uso: str = "") -> List[str]:
    note: List[str] = []
    beta = res["beta"]
    V = res["V"]
    Q_max = res["Q_max"]
    Q_full = res["Q_full"]
    Fr = res["Fr"]

    mat = MATERIALI_MANNING.get(materiale, {})
    V_max_mat = mat.get("V_max", 5.0)
    V_min = V_MIN_AUTOCIRCOLANTE.get(uso, 0.60)
    n_materiale = mat.get("n", None)

    if beta > 0.94:
        note.append(
            f"Grado di riempimento critico (beta = {beta:.3f} > 0.94): la portata supera il massimo "
            "teorico della sezione. Incrementare il diametro o la pendenza."
        )
    elif beta > 0.80:
        note.append(
            f"Grado di riempimento elevato (beta = {beta:.3f} > 0.80): verificare margine di sicurezza, "
            "ventilazione e condizioni di piena eccezionale."
        )
    elif beta < 0.10:
        note.append(
            f"Grado di riempimento molto basso (beta = {beta:.3f} < 0.10): la sezione e' "
            "sovradimensionata. Verificare le condizioni di moto a bassa portata."
        )

    if not math.isnan(V):
        if V < V_min:
            note.append(
                f"Velocita media V = {V:.3f} m/s inferiore al minimo autocircolante {V_min:.2f} m/s "
                f"per uso '{uso}'. Rischio deposito sedimenti; aumentare la pendenza o ridurre il diametro."
            )
        elif V > V_max_mat:
            note.append(
                f"Velocita media V = {V:.3f} m/s supera il limite del materiale ({V_max_mat:.1f} m/s per "
                f"{materiale}). Rischio di abrasione/erosione interna."
            )

    if n_materiale and abs(dati.n - n_materiale) > 0.003:
        note.append(
            f"Il valore di Manning n = {dati.n:.4f} differisce dal valore tipico del materiale "
            f"'{materiale}' (n_tipico = {n_materiale:.3f}). Verificare lo stato di conservazione."
        )

    if not math.isnan(Fr) and Fr > 1.0:
        note.append(
            f"Moto supercritico (Fr = {Fr:.3f} > 1.0): attenzione ai risalti idraulici e alle "
            "transizioni di sezione. Progettare con cautela gli elementi dissipatori."
        )

    ratio_Q = dati.Q / Q_max if Q_max > 0 else float("nan")
    note.append(
        f"Portata di progetto Q = {dati.Q:.5f} m3/s  |  "
        f"Q_max(beta=0.94) = {Q_max:.5f} m3/s  |  "
        f"Q_full(beta=1.0) = {Q_full:.5f} m3/s  |  Q/Q_max = {ratio_Q:.3f}."
    )

    if not any(k in " ".join(note[:-1]) for k in ("critico", "elevato", "basso", "inferiore",
                                                    "supera", "differisce", "supercritico")):
        note.insert(0,
                    "Tutti i parametri idraulici rientrano nei limiti normativi tipici. "
                    "Completare la verifica con l'analisi di piena e le condizioni al contorno.")
    return note
