from __future__ import annotations

import importlib.util
import json
import math
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]


def load_src():
    spec = importlib.util.spec_from_file_location("pompe_src", ROOT / "src.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def assert_close(name: str, actual: float, expected: float, tol: float) -> None:
    if math.isnan(actual) or abs(actual - expected) > tol:
        raise AssertionError(f"{name}: actual={actual!r}, expected={expected!r}, tol={tol}")


def main() -> None:
    src = load_src()
    bench = json.loads((ROOT / "test" / "benchmark" / "base.json").read_text(encoding="utf-8"))
    data = bench["input"]
    fluido = src.Fluido(**data["fluido"])
    suction = src.Linea(**data["suction"])
    discharge = src.Linea(**data["discharge"])
    res = src.tdh_pump(data["Q"], suction, discharge, fluido)
    actual = dict(res)
    actual["P_kW"] = src.potenza_pompa(data["Q"], res["H"], fluido.rho, data["eta"]) / 1000.0
    actual["NPSH_d"] = src.npsh_disponibile(
        data["z_serbatoio"], suction.z, res["hf_s"] + res["hK_s"],
        fluido.rho, p_vap=fluido.p_vap
    )
    actual["Ns"] = src.velocita_specifica_ns(data["n_rpm"], data["Q"], res["H"])
    tol = float(bench["abs_tolerance"])
    for key, expected in bench["expected"].items():
        assert_close(key, float(actual[key]), float(expected), tol)
    print("OK Pompe benchmark: base")


if __name__ == "__main__":
    main()
