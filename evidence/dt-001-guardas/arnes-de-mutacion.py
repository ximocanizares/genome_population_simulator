#!/usr/bin/env python3
"""Arnés de mutación: una mutación cada vez, con huella antes y después.

Uso: arnes-de-mutacion.py MUTACIONES.json ARBOL COPIA_DE_SEGURIDAD [SOLO_ESTE_ID]

Cada mutación:
  1. huella SHA-256 del fichero antes;
  2. sustitución única (falla si no aparece exactamente una vez);
  3. huella después, que tiene que diferir;
  4. ejecución de la prueba, que tiene que FALLAR;
  5. restauración desde la copia de seguridad fuera de git y comprobación de
     que la huella vuelve a ser la de antes.
"""

import hashlib
import os
import json
import pathlib
import shutil
import subprocess
import sys

ROOT = None
BACKUP = None

# Las dos rutas son **argumentos**, no constantes: este arnés muta el árbol que
# se le da y guarda la copia de seguridad donde se le diga. Con la ruta del
# repositorio escrita a fuego, quien lo ejecutara mañana mutaría el árbol real
# sin saberlo.


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if len(sys.argv) < 4:
        print(
            "uso: arnes-de-mutacion.py MUTACIONES.json ARBOL COPIA_DE_SEGURIDAD [SOLO_ESTE_ID]\n"
            "  ARBOL                 raíz del repositorio que se va a mutar\n"
            "  COPIA_DE_SEGURIDAD    directorio **fuera de git** donde guardar el original",
            file=sys.stderr,
        )
        return 2
    spec = json.loads(pathlib.Path(sys.argv[1]).read_text())
    global ROOT, BACKUP
    ROOT = pathlib.Path(sys.argv[2]).resolve()
    BACKUP = pathlib.Path(sys.argv[3]).resolve()
    if BACKUP.is_relative_to(ROOT):
        print("la copia de seguridad no puede vivir dentro del árbol", file=sys.stderr)
        return 2
    only = sys.argv[4] if len(sys.argv) > 4 else None
    BACKUP.mkdir(parents=True, exist_ok=True)
    touched = {m["file"] for m in spec["mutations"]}
    for name in touched:
        shutil.copy2(ROOT / name, BACKUP / pathlib.Path(name).name)

    results = []
    for mutation in spec["mutations"]:
        if only and mutation["id"] != only:
            continue
        path = ROOT / mutation["file"]
        before = sha(path)
        text = path.read_text()
        count = text.count(mutation["old"])
        if count != 1:
            results.append((mutation["id"], "ANCLA-NO-UNICA", f"apariciones={count}", "", ""))
            continue
        path.write_text(text.replace(mutation["old"], mutation["new"], 1))
        after = sha(path)
        assert after != before, f"{mutation['id']}: la mutación no mutó"
        command = [
            "cargo", "test", "--locked",
            "--test", mutation["target"], mutation["test"], "--", "--exact",
        ] if mutation.get("exact", True) else [
            "cargo", "test", "--locked", "--test", mutation["target"], mutation["test"],
        ]
        if mutation["target"] == "lib":
            command = ["cargo", "test", "--locked", "--lib", mutation["test"]]
        run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        passed = run.returncode == 0
        compiled = "error: could not compile" not in run.stderr and "error[E" not in run.stderr
        # Restaurar SIEMPRE, y desde fuera de git.
        # `copy2` conserva la fecha del original, y entonces `cargo` cree que
        # su binario sigue al día y **no recompila**: el árbol queda restaurado
        # pero la biblioteca compilada sigue mutada. Ya pasó una vez. Por eso
        # se toca la fecha después de restaurar.
        shutil.copy2(BACKUP / pathlib.Path(mutation["file"]).name, path)
        os.utime(path, None)
        restored = sha(path)
        assert restored == before, f"{mutation['id']}: la restauración no volvió al original"
        if not compiled:
            # Una mutación que no compila no demuestra nada: el fallo no es de
            # la prueba, es del compilador.
            verdict = "NO COMPILA (nula)"
        else:
            verdict = "SOBREVIVE (MALA)" if passed else "MUERE (buena)"
        evidence = ""
        if not passed and compiled:
            for line in run.stdout.splitlines():
                if "panicked at" in line or line.strip().startswith("test ") and "FAILED" in line:
                    evidence = line.strip()[:160]
                    break
            if not evidence:
                for line in run.stdout.splitlines():
                    if "assertion" in line or "left:" in line:
                        evidence = line.strip()[:160]
                        break
        results.append((mutation["id"], verdict, before[:16], after[:16], evidence))

    print(f"{'id':<52} {'veredicto':<18} {'sha-antes':<18} {'sha-después':<18}")
    bad = 0
    for row in results:
        print(f"{row[0]:<52} {row[1]:<18} {row[2]:<18} {row[3]:<18} {row[4]}")
        if "SOBREVIVE" in row[1] or "ANCLA" in row[1] or "NO COMPILA" in row[1]:
            bad += 1
    print(f"\ntotal={len(results)} fallidas={bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
