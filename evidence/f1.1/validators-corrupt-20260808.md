# Evidencia A7: rechazo de artefactos corruptos

- `evidence_version`: `1`
- Acta iniciada: 8 de agosto de 2026
- Repetición correctora: 9 de agosto de 2026
- Plataforma: Ubuntu 24.04 sobre WSL2, Linux
  `6.18.33.2-microsoft-standard-WSL2`, x86_64
- Sistema de ficheros: `ext4`, fuente `/dev/sdc`, ruta bajo
  `/home/ximo/devel/genome_population_simulator` y fuera de `/mnt/c`
- Comando: `./scripts/check-independent-validators.sh evidence/f1.1/validators-corrupt-20260809`
- Resultado: código de salida 0; las cinco corrupciones fueron rechazadas por
  la herramienta independiente prevista.
- Salidas y metadatos conservados:
  `evidence/f1.1/validators-corrupt-20260809/`. Cada clase contiene sus ficheros
  `.stdout`, `.stderr` y `.metadata`; las entradas ejecutadas también se
  conservan en `inputs/`.

| Clase | Herramienta | Código | SHA-256 de stdout | SHA-256 de stderr | SHA-256 del artefacto | Resultado observado |
| --- | --- | ---: | --- | --- | --- | --- |
| `fasta` | samtools 1.21 / htslib 1.21 | 1 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `41f1d20f55375b412fc2d245eecc9b9255ab9dbf6b7fd6f4c091d834955c7606` | `45a22e7909c678743900bb02cb1e3f45924e46bc685bbcc8f68d6501b78318f5` | Rechazo por estructura FASTA inválida |
| `sha256` | GNU sha256sum 9.4 | 1 | `516484dcc7c8874d5d653e31e07cd71ab7101b2d1086b19ea170b7acc5472d1b` | `c21162752e0e9447063c54afadaabcd769aa943d1dc1e9641847e990d83020b9` | `57eeabf40951f15e90223334c1cc380ccec58a791823f6b777ab6e312cae9f22` | Rechazo por discordancia de digest |
| `tar_gz` | GNU tar 1.35 | 2 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `04010bb3b0c6463e6b05a20987c8a12623f120f8318a2b6c60e3b6816b888d4a` | `9bc6b13f4bf21eeda475a6091c419a63c5b58efd4fe438a9d3921640b7fa219d` | Rechazo por no ser gzip |
| `json_schema` | CPython 3.12.3 / jsonschema 4.10.3 | 1 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `e777e75e5082976321d2a3e8c4e5976ecb9975981fb04f3f7204f3a674ac97d8` | `751dd6a80e5fba2a18df6ff84db84c09381a8348948ceed516385770beba31b2` | Rechazo del campo superior desconocido |
| `oci` | Docker CLI 29.6.1, build 8900f1d | 1 | `37517e5f3dc66819f61f5a7bb8ace1921282415f10551d2defa5c3eb0985b570` | `584919538cff73f6f0fe28764c29b7564702a15858c94cb6a456e4852f5090c6` | `0000000000000000000000000000000000000000000000000000000000000000` | Rechazo como imagen inexistente con daemon disponible |

Salida resumida conservada:

```text
PASS class=fasta status=1 stdout_sha256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 stderr_sha256=41f1d20f55375b412fc2d245eecc9b9255ab9dbf6b7fd6f4c091d834955c7606 artifact_sha256=45a22e7909c678743900bb02cb1e3f45924e46bc685bbcc8f68d6501b78318f5
PASS class=sha256 status=1 stdout_sha256=516484dcc7c8874d5d653e31e07cd71ab7101b2d1086b19ea170b7acc5472d1b stderr_sha256=c21162752e0e9447063c54afadaabcd769aa943d1dc1e9641847e990d83020b9 artifact_sha256=57eeabf40951f15e90223334c1cc380ccec58a791823f6b777ab6e312cae9f22
PASS class=tar_gz status=2 stdout_sha256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 stderr_sha256=04010bb3b0c6463e6b05a20987c8a12623f120f8318a2b6c60e3b6816b888d4a artifact_sha256=9bc6b13f4bf21eeda475a6091c419a63c5b58efd4fe438a9d3921640b7fa219d
PASS class=json_schema status=1 stdout_sha256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 stderr_sha256=e777e75e5082976321d2a3e8c4e5976ecb9975981fb04f3f7204f3a674ac97d8 artifact_sha256=751dd6a80e5fba2a18df6ff84db84c09381a8348948ceed516385770beba31b2
PASS class=oci status=1 stdout_sha256=37517e5f3dc66819f61f5a7bb8ace1921282415f10551d2defa5c3eb0985b570 stderr_sha256=584919538cff73f6f0fe28764c29b7564702a15858c94cb6a456e4852f5090c6 artifact_sha256=0000000000000000000000000000000000000000000000000000000000000000
```

La repetición también comprobó el caso de control de infraestructura: sin
acceso al socket del daemon, la clase OCI falla con “infraestructura no
disponible” y no se registra como rechazo correcto. Las salidas conservadas no
contienen secretos ni rutas temporales.
