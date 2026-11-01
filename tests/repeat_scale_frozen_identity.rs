//! **Huella de los escenarios con muchos microsatélites**, congelada para la
//! entrega 2 («escala») del plan de detección y **reescrita para la entrega 3**
//! («detección en el producto») sin recongelar lo que no tenía que cambiar.
//!
//! # Qué protege, dicho sin rodeos
//!
//! La entrega 2 cambió cómo se **recorren** los loci de repetición para que el
//! coste crezca en línea recta con su número: no podía cambiar un solo byte
//! publicado. La huella de `tests/repeat_detection_frozen_identity.rs` solo
//! cubre escenarios **sin** microsatélites; esta cubre los que tienen muchos.
//!
//! # Qué cambió con la entrega 3, y por qué los literales siguen siendo los de `dcae406`
//!
//! Desde la entrega 3 un locus declarado solo existe con la sección de
//! microsatélites (REQ-CONF-142) y tiene que coincidir con uno detectado
//! (REQ-REF-40). Los escenarios declaran ahora esa sección con **una región
//! por locus que cubre exactamente su tramo**, de modo que lo detectado es lo
//! declarado. Lo que eso cambia, y lo único:
//!
//! 1. **La copia del TOML** (`inputs/config/general.toml`), que ahora lleva la
//!    sección.
//! 2. **Lo que depende de la huella de ese TOML**: el identificador de
//!    ejecución (en la verdad, en el VCF y en el manifiesto) y la huella de la
//!    configuración que la verdad cita en cada fuente aleatoria.
//! 3. **Los campos nuevos** de la configuración efectiva (`repeat_detection` y
//!    el `origin` de cada locus) y del manifiesto (el `origin` de cada locus),
//!    y en el manifiesto la huella y el tamaño de los ficheros que cambian.
//!
//! Por eso cada escenario se juzga con **dos tablas**:
//!
//! - **La tabla normalizada, con los literales de `dcae406`.** Todos los
//!   ficheros publicados salvo la copia del TOML, con tres normalizaciones que
//!   deshacen exactamente los cambios 2 y 3: en la verdad y el VCF, el
//!   identificador de ejecución pasa a `RUN_NORMALIZED` y la huella del TOML a
//!   ceros; de la efectiva se quitan `repeat_detection` y el `origin` de cada
//!   locus; del manifiesto, el `origin` de cada locus, y el identificador de
//!   ejecución y la huella y el tamaño de la copia del TOML, la efectiva, la
//!   verdad y el VCF se normalizan igual. **Todo lo demás —los FASTQ, la
//!   matriz, la referencia, el manifiesto de evento y cada byte de la verdad,
//!   el VCF, la efectiva y el manifiesto que no sea uno de esos campos— tiene
//!   que salir idéntico al de `dcae406`.** Los literales se obtuvieron el
//!   30/10/2026 aplicando esta misma normalización (la prueba ignorada
//!   `print_normalized_tables_of_published_trees`) a los árboles que publicó el
//!   binario de `dcae406`, compilado aparte con `git archive` y su propio
//!   directorio de compilación, con los TOML de antes. Antes de normalizar se
//!   comprobó que esos árboles daban, tal cual, las huellas que `dcae406`
//!   tenía congeladas en esta misma prueba. En `dcae406` la efectiva y el
//!   manifiesto no tienen campos nuevos y la normalización solo toca lo que
//!   depende del TOML.
//! - **La tabla de los documentos que cambian a propósito** —efectiva,
//!   manifiesto y copia del TOML—, con huellas medidas con el código de la
//!   entrega 3. Congela los campos nuevos; no protege nada de lo de antes, que
//!   ya protege la primera tabla.
//!
//! **El diff de la efectiva y del manifiesto entre `dcae406` y la entrega 3**,
//! campo a campo, medido el 30/10/2026 en los escenarios 1 y 2: la efectiva
//! gana `repeat_detection` y `"origin": "declared"` en los quince loci, y
//! nada más; el manifiesto gana el mismo `origin` en los quince loci de
//! `repeat_loci.loci` y cambia el identificador de ejecución, la huella de la
//! copia del TOML, de la efectiva (en `artifacts` y en el enlace del evento
//! técnico), de la verdad y del VCF, y el tamaño de la copia del TOML y de la
//! efectiva (la verdad y el VCF miden lo mismo: el identificador nuevo tiene la
//! misma longitud). Nada más.
//!
//! # Qué cambió con la entrega 4 (01/11/2026), y por qué la tabla normalizada no se toca
//!
//! La entrega 4 publica además el **catálogo de microsatélites aparte**
//! (`repeat_catalog.tsv` y `repeat_catalog.metadata.json`, REQ-MAN-98) y lo
//! enlaza con su huella desde la efectiva. La tabla normalizada lo deshace
//! igual que deshace la entrega 3: deja fuera los dos ficheros y quita sus dos
//! artefactos de la lista y de `output_categories.minimum` del manifiesto; la
//! efectiva ya perdía `repeat_detection` entera. **Sus literales siguen siendo
//! los de `dcae406` y no se han tocado.** La tabla de los documentos que
//! cambian a propósito gana los dos ficheros del catálogo y se volvió a medir
//! con el código de la entrega 4. El diff de la efectiva y del manifiesto
//! frente a `45af26b`, comparado campo a campo sobre dos escenarios de la
//! revisión (`w-bound-0` y `v-136`) con el binario de `45af26b` compilado
//! aparte: la efectiva gana `repeat_detection.catalog` y pierde
//! `repeat_detection.retired_loci` (los retirados van al metadato); el
//! manifiesto gana los dos artefactos del catálogo en `artifacts` y en
//! `output_categories.minimum`, y cambian la huella y el tamaño de la
//! efectiva. Ningún otro fichero cambia. Los escenarios de la revisión sin la
//! sección no publican catálogo y su tabla de documentos no cambia.
//!
//! # Qué cambió con la entrega 4 bis (01/11/2026), y por qué la tabla normalizada tampoco se toca
//!
//! La entrega 4 bis saca los loci **detectados** de la efectiva y del
//! manifiesto y publica sus alelos en `repeat_alleles.tsv`, que se escribe
//! siempre con la sección (aquí, solo la cabecera: todos los loci de estos
//! escenarios son declarados). La tabla normalizada deja fuera también ese
//! fichero y su artefacto, como los del catálogo, y **sus literales siguen
//! siendo los de `dcae406`**. La tabla de los documentos gana la tabla de
//! alelos y se volvió a medir. Comparado campo a campo contra el binario de
//! `afda48b` compilado aparte, en los dos escenarios de muchos loci y en los
//! veintiuno de la revisión: en los dieciséis que llevan la sección, la
//! efectiva gana en `repeat_detection` el enlace `alleles_table`, las
//! constantes de `frequency_normalization` y, en `slippage`, las versiones,
//! el ámbito, la versión de su catálogo, las tres limitaciones y los linajes;
//! el manifiesto gana el artefacto `repeat-alleles-tsv` en `artifacts` y en
//! `output_categories.minimum`, y cambian la huella y el tamaño de la
//! efectiva (también en el enlace del evento técnico); aparece la tabla de
//! alelos de 203 bytes. Ningún otro fichero cambia: ni la verdad, ni el VCF,
//! ni los FASTQ, ni el manifiesto de evento, ni el catálogo. Los cinco
//! escenarios de la revisión sin la sección no cambian en nada.
//!
//! # Qué tiene cada escenario, y qué ejercita
//!
//! 1. **Muchos loci con tartamudeo de base y quimeras entre secuencias.**
//!    Referencia suministrada de dos secuencias con quince microsatélites
//!    plantados (nueve en `chrA` y seis en `chrB`), de periodos 1 a 4, dos
//!    parejas de loci **pegados** (el final de uno es el inicio del otro), la
//!    mayoría con alelos generados por deslizamiento y dos con alelos escritos
//!    a mano, de modo que los alelos varían entre individuos. Cuatro
//!    individuos diploides, VCF de texto y un evento técnico con tartamudeo de
//!    base encendido en todos los periodos y la mitad de las parejas
//!    quiméricas (DTC-102: el sufijo puede salir de la otra secuencia). La
//!    sección de microsatélites usa los suelos por omisión: los quince loci los
//!    alcanzan.
//! 2. **Lo mismo con PCR y tartamudeo clonal por linaje**, que sitúa los
//!    tractos por el mismo camino (`place_tracts`) desde la molécula original
//!    de la familia.
//! 3. **Los escenarios de la revisión adversarial de la fusión del VCF**
//!    (`tests/fixtures/vcf-fusion/review`, 2 de octubre de 2026). Sus TOML
//!    llevan desde la entrega 3 la sección que escribe
//!    `tests/support/detection_section.rs` —una región por locus con seis
//!    bases de margen y suelos bajos en las longitudes de motivo que declaran,
//!    porque varios declaran tractos cortos a propósito—, con un retoque a
//!    mano anotado en el propio TOML (`w-bound-0`: la región llega hasta el
//!    final del `(CA)×40` del que el locus declarado es el principio). Cada uno
//!    se congela con dos literales, el resumen de cada tabla.
//!
//! Además, en los escenarios 1 y 2 el **juez semántico del VCF**
//! (`validate_published_vcf_semantics`) tiene que aceptar lo publicado.
//!
//! # Qué discrimina
//!
//! Comprobado el 29/10/2026 para la entrega 2 (apartado 11.9 del plan): un
//! salto de línea de más al final de los JSON publicados mata las tres pruebas,
//! y los índices de la entrega 2 mutados para que pierdan lo que empieza antes
//! de la ventana matan al menos una. Para la entrega 3, comprobado el
//! 30/10/2026 mutando una copia aparte: un locus declarado con alelos escritos
//! a mano que perdiera sus alelos y heredara el deslizamiento común mata las
//! tres pruebas **en la tabla normalizada** (cambian la verdad, el VCF y los
//! FASTQ); uno que cediera su sitio al detectado que coincide mata las tres
//! antes, en el control de los orígenes o en la ejecución.
//!
//! # Qué hacer cuando esta prueba falle
//!
//! El mensaje imprime las tablas observadas. **No se copian sin más**: si
//! cambia la tabla normalizada, ha cambiado algo que la entrega 3 no podía
//! tocar.
//!
//! # Honestidad sobre el compresor
//!
//! Los `.gz` los comprime `noodles-bgzf`, que puede elegir rutina según el
//! procesador. Cada fichero comprimido lleva **dos** huellas: la del fichero y
//! la de su contenido descomprimido (entrada `#content`). Si solo fallara la
//! primera, el culpable sería el compresor y no el productor.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::results::json_document_bytes;
use genome_population_simulator::vcf_semantic_validation::validate_published_vcf_semantics;
use noodles_bgzf as bgzf;
use serde_json::Value;
use sha2::{Digest, Sha256};

const SEED: &str = "5e2a9c7b1d3f4068a2c4e6f8091b3d5f7a9c1e3f5072b4d6f8a0c2e4f6a8b0c2";
static NONCE: AtomicU64 = AtomicU64::new(0);

// ---------------------------------------------------------------------------
// Las tablas normalizadas, con los literales de `dcae406`, y las de los tres
// documentos que cambian a propósito, medidas con el código de la entrega 3.
// ---------------------------------------------------------------------------

const MANY_LOCI_BASE_STUTTER: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "ca72667bf66b0c52c1701bb6c19ded23616795b89619e007b6a0829c23f56095",
    ),
    (
        "events/reads-main/manifest.json",
        "b1694d258e89195fd398206c55c3a235769ec68cef96d01816320e6a98ff1320",
    ),
    (
        "manifest.json",
        "a3612ef43e8cbbad30e02030030886b0f406c4ac72041a2c9227c07a9f48b2d2",
    ),
    (
        "outputs/fastq/reads-main/S1_R1.fastq.gz",
        "b422901e6a74c89bc9d6fa500e5621617528f55d26547b982e5280672319b3d7",
    ),
    (
        "outputs/fastq/reads-main/S1_R1.fastq.gz#content",
        "e57f462c8f00b02ec294746c4e8687206a8873ca0a6cb35b23d68444ccdb4b9f",
    ),
    (
        "outputs/fastq/reads-main/S1_R2.fastq.gz",
        "d2e29a74e83e94453102f9f4493eb7c4b01d166e70939d614db757e00ac7cb4c",
    ),
    (
        "outputs/fastq/reads-main/S1_R2.fastq.gz#content",
        "2ab58fd86236904f32ec58540b617d9291ada0a26cab0384f34aa88cb36451bd",
    ),
    (
        "outputs/fastq/reads-main/S2_R1.fastq.gz",
        "7f67857ec426bce05bfd953d3ba669ec364e5089c2b5ba4e605c499604e88a99",
    ),
    (
        "outputs/fastq/reads-main/S2_R1.fastq.gz#content",
        "7408cb04fce5b9e7fb044d26d2369fcab37b09a24fb1fde4e90d4dea952a739b",
    ),
    (
        "outputs/fastq/reads-main/S2_R2.fastq.gz",
        "8a2b2a00b6428bd854adb3393d012372dda15fa63c1c22b5be2b6a74c6269c0e",
    ),
    (
        "outputs/fastq/reads-main/S2_R2.fastq.gz#content",
        "bbe088dae8f924deff6cd7d60ed4163516a30940b162bc1b876a1a5efe0a429c",
    ),
    (
        "outputs/fastq/reads-main/S3_R1.fastq.gz",
        "8ab8baedc2ed9bb921d1c444942cedfb3089fa7e344ec20049439f848fe93911",
    ),
    (
        "outputs/fastq/reads-main/S3_R1.fastq.gz#content",
        "73285913902cf44b3228fc85cbfd3918136d2e7c3694db5fa32d5af9ad790d69",
    ),
    (
        "outputs/fastq/reads-main/S3_R2.fastq.gz",
        "7c60fe5f249d5eee11480ae4a6b2b7439deec12f7a1e7a9744da3e199dd35b79",
    ),
    (
        "outputs/fastq/reads-main/S3_R2.fastq.gz#content",
        "7ac0fe582dd42faa27cdaca13b8927f55fe776318dca5b78346ce3e05d45285d",
    ),
    (
        "outputs/fastq/reads-main/S4_R1.fastq.gz",
        "e80d108d77f47d13ed8529558ebe0a704fb93fc3a74e3ee4810da4e3253caf13",
    ),
    (
        "outputs/fastq/reads-main/S4_R1.fastq.gz#content",
        "991f79ccfe47d22c348453a3a040b9b4a58c03b706c7f075946d9a7ac437da0e",
    ),
    (
        "outputs/fastq/reads-main/S4_R2.fastq.gz",
        "dc61707820f79b8ca7175213cac1316c435e744d218637771b255ab5ac968bb6",
    ),
    (
        "outputs/fastq/reads-main/S4_R2.fastq.gz#content",
        "f25b888ffc7471a2cd7a2b6fb29e0eca6eb95c305ca12789d2edd62a8f7e35ef",
    ),
    (
        "outputs/truth.jsonl",
        "cb56facbd9f8cc06ed9b027717b159346b2f01f60066606ad630c16cb6bcaa62",
    ),
    (
        "outputs/truth.vcf",
        "1abd7975c81243361f2c5991851a7a6e6f3dc262d8784cc2bd236a5cc37c9f0d",
    ),
    (
        "population_snvs.metadata.json",
        "408a28ebcba661a4fc9790a6fd3c18fdd7fd13de6042a87771e42b4510842a7e",
    ),
    (
        "population_snvs.tsv",
        "7fcce4afec94295fe656b24476f1915a1e2234e19d371d529f1c77d6a4c28851",
    ),
    (
        "resources/references/ref-main.fa",
        "04cd04f7f550804f32f80d6cae03895278689d34bbbb144d48486eca11a08c0a",
    ),
];

const MANY_LOCI_BASE_STUTTER_DOCUMENTS: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "5455adb3b55ce5586b1f9e1ef10a472c7e2b556243e59d44ab09eb3cb0b186c8",
    ),
    (
        "inputs/config/general.toml",
        "403cd6073aa4eea45942d9a45723745852652e1d4e6cc2bbb4444c9113ed3a3d",
    ),
    (
        "manifest.json",
        "9068907a483632a2319ce15894d50cc35ce36721e5302fd6e5cc0df03b969ba5",
    ),
    (
        "repeat_alleles.tsv",
        "c1dc0fe32075225c7d2d2cf9352ea8543c26033feeb54c2456ecc0b01b66ef6b",
    ),
    (
        "repeat_catalog.metadata.json",
        "0717a2c36e08d75b72c68b46fcf5698a93bba90ef93f796ac24dd3b3902b849b",
    ),
    (
        "repeat_catalog.tsv",
        "494ef1b9fcebf457a74da0858075a967f2ea8235ba83879f5499c6a5a9c3e56d",
    ),
];

const MANY_LOCI_PCR_CLONAL: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "3e508e8a14d5b92b1520e8db61183edca3d6ba07bca94980a1fb5e529b97d7e8",
    ),
    (
        "events/reads-main/manifest.json",
        "7cdfba62478ec28a8b184123c270a5ebec4815455aac5f7280fe3fd2f0f16b02",
    ),
    (
        "manifest.json",
        "107b44966888ffed875036275bb793a7e34c4805ce3ea02e2e9a4041e9f64b79",
    ),
    (
        "outputs/fastq/reads-main/S1_R1.fastq.gz",
        "0b0710fe8b0a9c87632f1424a3e4b2599e5d579af89bb527b881651a942371e2",
    ),
    (
        "outputs/fastq/reads-main/S1_R1.fastq.gz#content",
        "33117427107f7d1472b08059ec78a856cb57a1aeefac6ae4dbb33707dbf0edca",
    ),
    (
        "outputs/fastq/reads-main/S1_R2.fastq.gz",
        "32e4db035208399cd5ee6d62f1c212a1768a9e353434ee97775b562e4eb64bdf",
    ),
    (
        "outputs/fastq/reads-main/S1_R2.fastq.gz#content",
        "91f488ac56754ad228d2418771c83b900ef191e187e2c436e74044d61fbdd3a3",
    ),
    (
        "outputs/fastq/reads-main/S2_R1.fastq.gz",
        "1f4c8295173a5fe1cd4ab9d608bf4b8427931afe6a8be368be65ce6db5ca4cda",
    ),
    (
        "outputs/fastq/reads-main/S2_R1.fastq.gz#content",
        "1f4e646d2069bf636bc381f5a8b0b9b709fd5c393bf6c9e15cba37d1c29259c9",
    ),
    (
        "outputs/fastq/reads-main/S2_R2.fastq.gz",
        "bac6ffc672ae5f30d1ca3912b9d7cd8dc0941bd6e75b2096a407ae38612a69b8",
    ),
    (
        "outputs/fastq/reads-main/S2_R2.fastq.gz#content",
        "945980642ea4d49417b524c6e2dcb3028a8c5ceb4a49abc014fe76259ca638c2",
    ),
    (
        "outputs/fastq/reads-main/S3_R1.fastq.gz",
        "febfd6d7a694661ce96f4b8647594c34fc0f6070ff54b9ea7db222c29261c99f",
    ),
    (
        "outputs/fastq/reads-main/S3_R1.fastq.gz#content",
        "3d1d96338aa6b50c7d24cdbb12a820336ea68ab46d1df8a3d5715214e993006a",
    ),
    (
        "outputs/fastq/reads-main/S3_R2.fastq.gz",
        "70974533206f6f91a905567f33e54b5e7550da23594901574564a8749e08bfb6",
    ),
    (
        "outputs/fastq/reads-main/S3_R2.fastq.gz#content",
        "5bb5842829bec4d781150cabadda9c5eba69ce44e3a9eadbd2f548e51e816338",
    ),
    (
        "outputs/fastq/reads-main/S4_R1.fastq.gz",
        "bdfc553d6b7bf675e33cc492a2cf9a41245ddb8fb32cd9f6a3b3aef9cd62209b",
    ),
    (
        "outputs/fastq/reads-main/S4_R1.fastq.gz#content",
        "1d8152747c04cbe7a34d64a77b7ca5bc058ea85c05c6a274c30fa330a914d8b8",
    ),
    (
        "outputs/fastq/reads-main/S4_R2.fastq.gz",
        "021780cc7ebaa3bdaac7b0892801f34afc712d6aa89b8714bddd3074e33c48c2",
    ),
    (
        "outputs/fastq/reads-main/S4_R2.fastq.gz#content",
        "ed542a8864e285dbe54b570fa0efca19c9e1d99acd71a83d92cf83a6dbb7dc0c",
    ),
    (
        "outputs/truth.jsonl",
        "cb56facbd9f8cc06ed9b027717b159346b2f01f60066606ad630c16cb6bcaa62",
    ),
    (
        "outputs/truth.vcf",
        "1abd7975c81243361f2c5991851a7a6e6f3dc262d8784cc2bd236a5cc37c9f0d",
    ),
    (
        "population_snvs.metadata.json",
        "408a28ebcba661a4fc9790a6fd3c18fdd7fd13de6042a87771e42b4510842a7e",
    ),
    (
        "population_snvs.tsv",
        "7fcce4afec94295fe656b24476f1915a1e2234e19d371d529f1c77d6a4c28851",
    ),
    (
        "resources/references/ref-main.fa",
        "04cd04f7f550804f32f80d6cae03895278689d34bbbb144d48486eca11a08c0a",
    ),
];

const MANY_LOCI_PCR_CLONAL_DOCUMENTS: &[(&str, &str)] = &[
    (
        "effective-config.json",
        "deb390d8f496e75858e6ff51ba9b12891b47093ee201ad2a1ee833d79f190595",
    ),
    (
        "inputs/config/general.toml",
        "027cb66eecf74e61660fb04c37ad0d28b99f452f422f94b307305e01f1e4fd05",
    ),
    (
        "manifest.json",
        "815dd89c2c3b8c1d4ef60b80054f29568f23d8a5bb295613e4b0e3ba9d82a87f",
    ),
    (
        "repeat_alleles.tsv",
        "c1dc0fe32075225c7d2d2cf9352ea8543c26033feeb54c2456ecc0b01b66ef6b",
    ),
    (
        "repeat_catalog.metadata.json",
        "0717a2c36e08d75b72c68b46fcf5698a93bba90ef93f796ac24dd3b3902b849b",
    ),
    (
        "repeat_catalog.tsv",
        "494ef1b9fcebf457a74da0858075a967f2ea8235ba83879f5499c6a5a9c3e56d",
    ),
];

/// Resumen de la tabla normalizada de cada escenario de la revisión de la
/// fusión del VCF: SHA-256 de `ruta\tSHA-256\n` ordenada. Literales de
/// `dcae406`.
const REVIEW_FIXTURES: &[(&str, &str)] = &[
    (
        "20248",
        "7793fa1fdfaf6bcbecb2cf250d5fb0a46e3244217234f259bf1bbd02074991d4",
    ),
    (
        "40045",
        "b393aa500f566c6a7c2387238c9a3b96dd6e2886d50848d82512bcdd78b74263",
    ),
    (
        "impure-tract",
        "412a668efe858cb6b8135c020aa87420622e4fb2e8f4b9a8eeeea11a5b34b651",
    ),
    (
        "r03",
        "e9af72876afa8bad10d72b08dcd805f3d4ae62a156e0e7df56504ecdf62de4b9",
    ),
    (
        "rb-over-alt-1831",
        "9b477ba46a984344a6808858508b8aa24b90550626b20ceb809677a41c92b43a",
    ),
    (
        "v-136",
        "0f7ba17e7770cb74fe7a1b93d953d234e10c07a828a2a2d34859b7a4a6f89d6a",
    ),
    (
        "v-163",
        "57a029dffc4e7e0064f1adf50309fa7619d5e5502bdd34509b40434b16d43801",
    ),
    (
        "v-175",
        "b4bdfe2f92b3fb869867654ee706f2161fe57be3e1319fcd443b5bc9803307ee",
    ),
    (
        "w-atreg-4",
        "975a0c02fc6d60ecdb81320beb69d6200da79a5e23e22c39655e830087e7586c",
    ),
    (
        "w-bound-0",
        "1d8e0cc8fde4a0c63933546be353586befa62a6ac53bc79c953d1a8f4040623c",
    ),
    (
        "w-cap-2090",
        "c5642507fe7f8aa323d04bf39fde03ceb3af19a1064c49f18585848257bb6477",
    ),
    (
        "w-flankdel-5",
        "b019fbc52edc3020aa2be975b79a690b96ab7203c690f7528f7686676e994c43",
    ),
    (
        "w-min-a",
        "4a74a075d0c3adeb15f9332465ac584714bcaef1e23d4b189efcdc62d08cad96",
    ),
    (
        "w-plain-cross",
        "7baad90e3130e8fbae9f7f2ca1d7c31b70a89db5544ddb44088945474bb547a9",
    ),
    (
        "x-11063",
        "61a9b2a47203316422d24ba2d17469c56469b53239f72ff5051c6dd71eea4ec2",
    ),
    (
        "x-12168",
        "9f0ca1c82076b1f980784c9415b26eb8fc2456538dd22e3d99caafca9ed07210",
    ),
    (
        "x-1392",
        "11cc7436dfd041ceda2c2f13165cfa57169cee0b7bc1d8939822ce9c5f79f7ff",
    ),
    (
        "x-14771",
        "12c20bd917ed6f1b73d6b284ba1fcc07ed513ba21e7aad9defe081312206950b",
    ),
    (
        "x3338",
        "9dd80a53c12aaa8e79c0718fd49b8e90d5c160b97e36aee57af5c80ac86407da",
    ),
    (
        "y-1699",
        "489b56bb17f619fd15da77f5ee814597aacd1b240d061703bf4bff98a698a1bd",
    ),
    (
        "y-2813",
        "7466ddcc43f3f4be053de6649f2aa67cb1a1b36ea13b37dcd5ef1b998ee59785",
    ),
];

/// Resumen de los documentos que cambian a propósito en cada escenario
/// de la revisión, medido con el código de la entrega 3.
const REVIEW_FIXTURES_DOCUMENTS: &[(&str, &str)] = &[
    (
        "20248",
        "e1bdf8af95adaf2c634eb7e1948c082c4b2a1ee4041f0591e8d49b08d45415ca",
    ),
    (
        "40045",
        "05346806d2250e1c78bee5499ac0117b5cf44517e34b77a01bd25b6c0408e8e3",
    ),
    (
        "impure-tract",
        "c0994d372c5b318188a81ef81db321cb2e1636eeb5cab901d6bb0a40be4352c5",
    ),
    (
        "r03",
        "e285ca413ed522d074bb0357c6ceac520ca12ec57db81cc9644c8c7524dbe5d3",
    ),
    (
        "rb-over-alt-1831",
        "7f35a77477af258a6196580385779cf8c90d69909060703b3dec11a0b6668522",
    ),
    (
        "v-136",
        "862e3df321b1feb44e1e3ff55c789ec69ce332b65bf386458fee0c4e70af9dc1",
    ),
    (
        "v-163",
        "b9f20ec13651ed0babfcde92c5686f65be3224f94b10139c917cd3688f4110fe",
    ),
    (
        "v-175",
        "ec39fc6ace4e1355b9a7a85c37c007af883c891f5a5f12de711f17cd536262a4",
    ),
    (
        "w-atreg-4",
        "b7d8f7fdb963d3344b1cc9f08e5b0b7a4f5a62edc8ddd7ee4e7d489bca990a8f",
    ),
    (
        "w-bound-0",
        "ac2e0a590bc78c74b3f7448af49ceef5e0e1040a17dba303f3516a610afef254",
    ),
    (
        "w-cap-2090",
        "22dc19d468cc72726f984fcfa6dbafad40ce1b0a1b63b2ea04e4764591834206",
    ),
    (
        "w-flankdel-5",
        "59492afd7ea67347ad2f55e95ce67e0fd58dd4e915dd0e22f99cb18f99744ffb",
    ),
    (
        "w-min-a",
        "93a71466a646763bbe5aa5e09caadba134feb3b61b1323c3249fa5f701015b17",
    ),
    (
        "w-plain-cross",
        "11caf4f99a3b09c84ff34d7dd2900919cb15a3d89655b5f4e5f7e170f6c6136a",
    ),
    (
        "x-11063",
        "f7ae02ebf767a7874679160df5bbaf17b1132f130f836c71dfb1564d17f9e609",
    ),
    (
        "x-12168",
        "802a11c63b658cd306188dedd9f231fe5f6935d0f69deb6661aa7a1f38262078",
    ),
    (
        "x-1392",
        "fdd525a734c03853c620f22e7a733977a4279943cd64ba9207231c41dc6fff67",
    ),
    (
        "x-14771",
        "3dcb728d2ed72c046d01d61b9e045f88d75a72ce0ac08ae1dfeaab84f1cc4935",
    ),
    (
        "x3338",
        "53ecc80bb29377c343192db8418e8461c4f439a09fddf7d9988c4ebcfc605983",
    ),
    (
        "y-1699",
        "b1a4e40731ada302fb43c7100d9704b540841312f3ced526a93fc206e761773e",
    ),
    (
        "y-2813",
        "e048940ca1a97ec0a56e374a7c24f624a43e39eb900220413a94c1ddfc0f9ab0",
    ),
];

fn temporary() -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-huella-escala-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("reloj")
            .as_nanos()
    ))
}

/// Directorio temporal de un escenario: se borra al soltarse.
struct TemporaryRoot(PathBuf);

impl Drop for TemporaryRoot {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Bases pseudoaleatorias reproducibles (xorshift64), escritas aquí y sin
/// depender del programa.
fn bases(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("dos bits")]
        })
        .collect()
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// Un microsatélite plantado: `(secuencia, inicio en base 0, motivo,
/// unidades, alelos escritos a mano)`. Sin alelos escritos, los genera el
/// deslizamiento.
struct Planted {
    id: &'static str,
    sequence: &'static str,
    start: usize,
    motif: &'static str,
    units: usize,
    alleles: &'static [(usize, f64)],
}

impl Planted {
    fn end(&self) -> usize {
        self.start + self.motif.len() * self.units
    }
}

/// Quince loci en dos secuencias. `STR-A4`/`STR-A5` y `STR-B3`/`STR-B4` van
/// pegados: el primero acaba justo donde empieza el segundo.
const PLANTED: &[Planted] = &[
    Planted {
        id: "STR-A1",
        sequence: "chrA",
        start: 400,
        motif: "CA",
        units: 15,
        alleles: &[],
    },
    Planted {
        id: "STR-A2",
        sequence: "chrA",
        start: 900,
        motif: "GAT",
        units: 10,
        alleles: &[],
    },
    Planted {
        id: "STR-A3",
        sequence: "chrA",
        start: 1500,
        motif: "AAAG",
        units: 8,
        alleles: &[],
    },
    Planted {
        id: "STR-A4",
        sequence: "chrA",
        start: 2000,
        motif: "TG",
        units: 12,
        alleles: &[],
    },
    Planted {
        id: "STR-A5",
        sequence: "chrA",
        start: 2024,
        motif: "AGC",
        units: 9,
        alleles: &[],
    },
    Planted {
        id: "STR-A6",
        sequence: "chrA",
        start: 2600,
        motif: "A",
        units: 14,
        alleles: &[],
    },
    Planted {
        id: "STR-A7",
        sequence: "chrA",
        start: 3200,
        motif: "CA",
        units: 20,
        alleles: &[(20, 0.4), (18, 0.3), (23, 0.3)],
    },
    Planted {
        id: "STR-A8",
        sequence: "chrA",
        start: 3800,
        motif: "TTC",
        units: 9,
        alleles: &[],
    },
    Planted {
        id: "STR-A9",
        sequence: "chrA",
        start: 4400,
        motif: "GT",
        units: 10,
        alleles: &[],
    },
    Planted {
        id: "STR-B1",
        sequence: "chrB",
        start: 300,
        motif: "CA",
        units: 16,
        alleles: &[],
    },
    Planted {
        id: "STR-B2",
        sequence: "chrB",
        start: 800,
        motif: "ATG",
        units: 9,
        alleles: &[],
    },
    Planted {
        id: "STR-B3",
        sequence: "chrB",
        start: 1400,
        motif: "GA",
        units: 12,
        alleles: &[],
    },
    Planted {
        id: "STR-B4",
        sequence: "chrB",
        start: 1424,
        motif: "C",
        units: 12,
        alleles: &[],
    },
    Planted {
        id: "STR-B5",
        sequence: "chrB",
        start: 2200,
        motif: "AATG",
        units: 7,
        alleles: &[],
    },
    Planted {
        id: "STR-B6",
        sequence: "chrB",
        start: 2900,
        motif: "TC",
        units: 14,
        alleles: &[(14, 0.5), (11, 0.25), (16, 0.25)],
    },
];

/// Las dos secuencias con sus tractos plantados. La base anterior y la
/// siguiente de cada tracto se eligen distintas de la primera y la última del
/// motivo, para que el tracto no se prolongue por azar en los flancos.
fn reference() -> Vec<(&'static str, Vec<u8>)> {
    let mut sequences = vec![
        ("chrA", bases(0x2F6B_1C9D_83A5_47E1, 5000)),
        ("chrB", bases(0x7A3D_5E91_C2B4_6F08, 3500)),
    ];
    for planted in PLANTED {
        let bases = &mut sequences
            .iter_mut()
            .find(|(id, _)| *id == planted.sequence)
            .expect("secuencia del locus")
            .1;
        let tract = planted.motif.repeat(planted.units);
        bases[planted.start..planted.end()].copy_from_slice(tract.as_bytes());
        let motif = planted.motif.as_bytes();
        let avoid = |excluded: &[u8]| -> u8 {
            *b"ACGT"
                .iter()
                .find(|base| !excluded.contains(base))
                .expect("una base libre")
        };
        // El flanco izquierdo de un locus pegado es el tracto anterior: no se
        // toca.
        if !PLANTED
            .iter()
            .any(|other| other.sequence == planted.sequence && other.end() == planted.start)
        {
            bases[planted.start - 1] = avoid(&[motif[motif.len() - 1], motif[0]]);
        }
        if !PLANTED
            .iter()
            .any(|other| other.sequence == planted.sequence && other.start == planted.end())
        {
            bases[planted.end()] = avoid(&[motif[0], motif[motif.len() - 1]]);
        }
    }
    sequences
}

fn fasta_of(sequences: &[(&str, Vec<u8>)]) -> String {
    let mut text = String::new();
    for (id, bases) in sequences {
        text.push('>');
        text.push_str(id);
        text.push('\n');
        for chunk in bases.chunks(70) {
            text.push_str(std::str::from_utf8(chunk).expect("ASCII"));
            text.push('\n');
        }
    }
    text
}

fn loci_toml() -> String {
    let mut text = String::new();
    for planted in PLANTED {
        text.push_str(&format!(
            "[[reference.repeat_loci]]\nid = \"{}\"\nsequence = \"{}\"\nstart = {}\nend = {}\n\
             motif = \"{}\"\n\n",
            planted.id,
            planted.sequence,
            planted.start + 1,
            planted.end(),
            planted.motif
        ));
        if planted.alleles.is_empty() {
            text.push_str(&format!(
                "max_units = {}\n\n[reference.repeat_loci.slippage]\nmodel = \"stepwise\"\n\
                 mean_events = 2.0\nexpansion_bias = 0.5\n\n",
                planted.units + 10
            ));
        } else {
            for (units, frequency) in planted.alleles {
                text.push_str(&format!(
                    "[[reference.repeat_loci.alleles]]\nunits = {units}\nfrequency = {frequency}\n\n"
                ));
            }
        }
    }
    // **La sección de microsatélites** (entrega 3), con los suelos por
    // omisión y una región por locus que cubre exactamente su tramo: lo
    // detectado es lo declarado. El deslizamiento común no lo usa ningún locus,
    // porque todos declaran sus alelos o su generador.
    text.push_str(
        "[reference.repeat_detection]\n\n[reference.repeat_detection.slippage]\n\
         model = \"stepwise\"\nmean_events = 2.0\nexpansion_bias = 0.5\n\n",
    );
    for planted in PLANTED {
        text.push_str(&format!(
            "[[reference.repeat_detection.regions]]\nsequence = \"{}\"\nstart = {}\nend = {}\n\n",
            planted.sequence,
            planted.start + 1,
            planted.end()
        ));
    }
    text
}

const INDIVIDUALS: [&str; 4] = ["S1", "S2", "S3", "S4"];
const UDIS: [(&str, &str); 4] = [
    ("ACGTACGT", "TGCATGCA"),
    ("CGTACGTA", "GCATGCAT"),
    ("GTACGTAC", "CATGCATG"),
    ("TACGTACG", "ATGCATGC"),
];

const STUTTER: &str = "[technical_event.stutter]\nenabled = true\n\
     [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 1\nrate = 0.4\n\
     [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 2\nrate = 0.4\n\
     [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 3\nrate = 0.3\n\
     [[technical_event.stutter.motif_length_rates]]\nmotif_length_bases = 4\nrate = 0.3\n\
     [[technical_event.stutter.steps]]\nunits = -3\nweight = 0.05\n\
     [[technical_event.stutter.steps]]\nunits = -2\nweight = 0.15\n\
     [[technical_event.stutter.steps]]\nunits = -1\nweight = 0.6\n\
     [[technical_event.stutter.steps]]\nunits = 1\nweight = 0.2\n\n";

const PCR: &str = "library_preparation = \"pcr\"\n\n\
     [technical_event.pcr]\nduplicate_fraction = 0.5\ncycles = 10\nefficiency = 0.9\n\
     [technical_event.pcr.stutter]\n\
     [[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 1\nprobability = 0.05\n\
     [[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 2\nprobability = 0.05\n\
     [[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 3\nprobability = 0.03\n\
     [[technical_event.pcr.stutter.per_cycle_slip]]\nmotif_length_bases = 4\nprobability = 0.03\n\n";

/// Escribe el escenario 1 o 2 y devuelve la ruta de la configuración.
fn many_loci_scenario(root: &Path, preparation: &str) -> PathBuf {
    let sequences = reference();
    fs::write(root.join("reference.fa"), fasta_of(&sequences)).expect("FASTA");
    // Filas de la matriz lejos de todo tracto; la base de referencia se lee de
    // la secuencia.
    let mut matrix = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
    for (contig, position) in [("chrA", 150), ("chrA", 2300), ("chrB", 1100)] {
        let reference = sequences
            .iter()
            .find(|(id, _)| *id == contig)
            .expect("secuencia")
            .1[position - 1];
        let alternate = *b"ACGT"
            .iter()
            .find(|base| **base != reference)
            .expect("alternativa");
        matrix.push_str(&format!(
            "{contig}\t{position}\t{}\t{}\tP1\t0.5\n",
            char::from(reference),
            char::from(alternate)
        ));
    }
    fs::write(root.join("population-snvs.tsv"), matrix).expect("matriz");
    let declared: String = sequences
        .iter()
        .map(|(id, bases)| {
            format!(
                "[[reference.sequences]]\nid = \"{id}\"\nclass = \"chromosome\"\nlength = {}\n\n",
                bases.len()
            )
        })
        .collect();
    let quoted = INDIVIDUALS
        .iter()
        .map(|id| format!("\"{id}\""))
        .collect::<Vec<_>>()
        .join(", ");
    let libraries: String = INDIVIDUALS
        .iter()
        .zip(UDIS)
        .map(|(id, (i7, i5))| {
            format!(
                "[[technical_event.libraries]]\nindividual_id = \"{id}\"\ni7 = \"{i7}\"\ni5 = \"{i5}\"\n\n"
            )
        })
        .collect();
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"population-snvs.tsv\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             {declared}{}\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 4\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [{quoted}]\n\n\
             [vcf_output]\nmode = \"text\"\n\n\
             [technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [{quoted}]\n\
             target_coverage_wgs = 12\nchimera_rate = 0.5\n{preparation}{STUTTER}{libraries}",
            loci_toml()
        ),
    )
    .expect("configuración");
    config
}

/// Ejecuta `run` con el binario del proyecto, **desde el directorio del
/// escenario y con rutas relativas**, y devuelve el directorio publicado.
fn run_binary(root: &Path, config: &str) -> PathBuf {
    let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
        .current_dir(root)
        .arg("run")
        .arg(config)
        .arg("--output")
        .arg("published")
        .output()
        .expect("binario del simulador");
    assert!(
        output.status.success(),
        "la ejecución falló en {}\n{}\n{}",
        root.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    root.join("published")
}

fn walk(root: &Path, directory: &Path, found: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(directory).expect("directorio publicado") {
        let path = entry.expect("entrada").path();
        if path.is_dir() {
            walk(root, &path, found);
        } else {
            found.push(path.strip_prefix(root).expect("relativa").to_path_buf());
        }
    }
}

/// El contenido descomprimido de un fichero BGZF.
fn decompressed(path: &Path) -> Option<Vec<u8>> {
    let mut text = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).ok()?)
        .read_to_end(&mut text)
        .ok()?;
    Some(text)
}

/// Los documentos que cambian a propósito: los tres de la entrega 3, los dos
/// del catálogo de microsatélites de la entrega 4 y la tabla de alelos de la
/// entrega 4 bis.
const DOCUMENTS: [&str; 6] = [
    "effective-config.json",
    "inputs/config/general.toml",
    "manifest.json",
    "repeat_alleles.tsv",
    "repeat_catalog.metadata.json",
    "repeat_catalog.tsv",
];

/// Los dos ficheros del catálogo aparte (entrega 4) y la tabla de alelos de
/// los detectados (entrega 4 bis), que la tabla normalizada deja fuera y
/// cuyos artefactos quita del manifiesto.
const CATALOG_FILES: [&str; 3] = [
    "repeat_alleles.tsv",
    "repeat_catalog.metadata.json",
    "repeat_catalog.tsv",
];
const CATALOG_ARTIFACTS: [&str; 3] = [
    "repeat-alleles-tsv",
    "repeat-catalog-metadata",
    "repeat-catalog-tsv",
];

/// Ficheros del manifiesto cuya huella y tamaño dependen de la copia del TOML.
const TOML_DEPENDENT: [&str; 4] = [
    "effective-config.json",
    "inputs/config/general.toml",
    "outputs/truth.jsonl",
    "outputs/truth.vcf",
];

/// Quita el `origin` de cada locus de una lista de loci publicada.
fn strip_origins(loci: Option<&mut Value>) {
    for locus in loci.and_then(Value::as_array_mut).into_iter().flatten() {
        if let Some(object) = locus.as_object_mut() {
            object.remove("origin");
        }
    }
}

/// Pone a cero la huella y el tamaño de todo enlace del manifiesto a un
/// fichero que depende de la copia del TOML, esté donde esté: la efectiva, por
/// ejemplo, se cita en el catálogo y en cada evento técnico.
fn normalize_toml_dependent_links(value: &mut Value) {
    match value {
        Value::Object(object) => {
            if object
                .get("path")
                .and_then(Value::as_str)
                .is_some_and(|path| TOML_DEPENDENT.contains(&path))
                && object.contains_key("checksum")
            {
                object["checksum"]["digest"] = Value::from("0".repeat(64));
                object.insert("size_bytes".to_owned(), Value::from(0));
            }
            for child in object.values_mut() {
                normalize_toml_dependent_links(child);
            }
        }
        Value::Array(items) => items.iter_mut().for_each(normalize_toml_dependent_links),
        _ => {}
    }
}

/// **La tabla normalizada**: todos los ficheros publicados salvo la copia del
/// TOML, con la verdad y el VCF sin el identificador de ejecución ni la huella
/// del TOML, la efectiva sin sus campos nuevos y el manifiesto sin los suyos ni
/// lo que depende del TOML. Ver la cabecera.
fn normalized_table(destination: &Path) -> BTreeMap<String, String> {
    let mut paths = Vec::new();
    walk(destination, destination, &mut paths);
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join("manifest.json")).expect("manifiesto"))
            .expect("JSON");
    let execution_id = manifest["execution_id"]
        .as_str()
        .expect("identificador de ejecución")
        .to_owned();
    let toml_digest =
        digest(&fs::read(destination.join("inputs/config/general.toml")).expect("copia del TOML"));
    let mut table = BTreeMap::new();
    for relative in paths {
        let full = destination.join(&relative);
        let name = relative.to_string_lossy().into_owned();
        if name == "inputs/config/general.toml" || CATALOG_FILES.contains(&name.as_str()) {
            continue;
        }
        let bytes = fs::read(&full).expect("fichero");
        let normalized = match name.as_str() {
            "outputs/truth.jsonl" | "outputs/truth.vcf" => String::from_utf8(bytes)
                .expect("texto")
                .replace(&execution_id, "RUN_NORMALIZED")
                .replace(&toml_digest, &"0".repeat(64))
                .into_bytes(),
            "effective-config.json" => {
                let mut document: Value = serde_json::from_slice(&bytes).expect("JSON");
                if let Some(object) = document.as_object_mut() {
                    object.remove("repeat_detection");
                }
                strip_origins(document.get_mut("repeat_loci"));
                json_document_bytes(&document).expect("JSON")
            }
            "manifest.json" => {
                let mut document: Value = serde_json::from_slice(&bytes).expect("JSON");
                document["execution_id"] = Value::from("RUN_NORMALIZED");
                normalize_toml_dependent_links(&mut document);
                for pointer in ["/artifacts", "/output_categories/minimum"] {
                    if let Some(list) = document.pointer_mut(pointer).and_then(Value::as_array_mut)
                    {
                        list.retain(|entry| {
                            let id = entry.get("id").unwrap_or(entry);
                            !CATALOG_ARTIFACTS.iter().any(|artifact| id == artifact)
                        });
                    }
                }
                strip_origins(
                    document
                        .get_mut("repeat_loci")
                        .and_then(|section| section.get_mut("loci")),
                );
                json_document_bytes(&document).expect("JSON")
            }
            _ => bytes,
        };
        table.insert(name.clone(), digest(&normalized));
        if name.ends_with(".gz") {
            let content =
                decompressed(&full).unwrap_or_else(|| panic!("{name}: no es un BGZF legible"));
            table.insert(format!("{name}#content"), digest(&content));
        }
    }
    table
}

/// Huellas tal cual de los documentos que cambian a propósito. El catálogo
/// solo existe con la sección de microsatélites: sin ella no hay entrada.
fn documents_table(destination: &Path) -> BTreeMap<String, String> {
    DOCUMENTS
        .iter()
        .filter(|name| !CATALOG_FILES.contains(name) || destination.join(name).exists())
        .map(|name| {
            (
                (*name).to_owned(),
                digest(&fs::read(destination.join(name)).expect("documento")),
            )
        })
        .collect()
}

fn summary(table: &BTreeMap<String, String>) -> String {
    digest(
        table
            .iter()
            .map(|(path, sha)| format!("{path}\t{sha}\n"))
            .collect::<String>()
            .as_bytes(),
    )
}

fn expected(frozen: &[(&str, &str)]) -> BTreeMap<String, String> {
    frozen
        .iter()
        .map(|(path, sha)| ((*path).to_owned(), (*sha).to_owned()))
        .collect()
}

fn render(table: &BTreeMap<String, String>) -> String {
    table
        .iter()
        .map(|(path, sha)| format!("    (\"{path}\", \"{sha}\"),\n"))
        .collect()
}

/// El juez semántico del VCF acepta lo publicado.
fn assert_vcf_judge_accepts(destination: &Path) {
    validate_published_vcf_semantics(
        &destination.join("resources/references/ref-main.fa"),
        &destination.join("outputs/truth.jsonl"),
        &destination.join("outputs/truth.vcf"),
    )
    .unwrap_or_else(|error| panic!("el juez semántico rechaza el VCF publicado: {error:?}"));
}

fn check_many_loci(
    name: &str,
    frozen: &[(&str, &str)],
    documents: &[(&str, &str)],
    preparation: &str,
) {
    let root = temporary();
    fs::create_dir_all(&root).expect("raíz del escenario");
    let _cleanup = TemporaryRoot(root.clone());
    many_loci_scenario(&root, preparation);
    let destination = run_binary(&root, "general.toml");
    let observed = normalized_table(&destination);
    // Control positivo: el escenario publica de verdad el VCF, la verdad y los
    // FASTQ de las cuatro individuas, y el VCF tiene registros de repetición
    // de las dos secuencias.
    for needle in ["outputs/truth.vcf", "outputs/truth.jsonl", "S1_R1", "S4_R2"] {
        assert!(
            observed.keys().any(|path| path.contains(needle)),
            "{name}: el escenario tiene que publicar «{needle}»; publicó\n{}",
            render(&observed)
        );
    }
    let vcf = fs::read_to_string(destination.join("outputs/truth.vcf")).expect("VCF");
    for sequence in ["chrA", "chrB"] {
        let repeats = vcf
            .lines()
            .filter(|line| line.starts_with(&format!("{sequence}\t")) && line.contains("RN="))
            .count();
        assert!(
            repeats >= 4,
            "{name}: el VCF tiene que llevar registros de repetición en {sequence}: {repeats}"
        );
    }
    let manifest = fs::read_to_string(destination.join("events/reads-main/manifest.json"))
        .expect("manifiesto del evento");
    assert!(
        manifest.contains("stutter"),
        "{name}: el evento técnico tiene que tartamudear"
    );
    // Control positivo de la entrega 3: los quince loci se publican con su
    // origen, todos declarados, y la sección de detección está.
    let effective: Value = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("efectiva"),
    )
    .expect("JSON");
    let origins: Vec<&str> = effective["repeat_loci"]
        .as_array()
        .expect("loci")
        .iter()
        .map(|locus| locus["origin"].as_str().expect("origen"))
        .collect();
    assert_eq!(origins, vec!["declared"; PLANTED.len()], "{name}");
    assert!(effective.get("repeat_detection").is_some(), "{name}");
    assert_vcf_judge_accepts(&destination);
    assert_eq!(
        observed,
        expected(frozen),
        "{name}: la tabla normalizada no es la de dcae406. La entrega 3 solo podía cambiar la \
         copia del TOML, lo que depende de su huella y los campos nuevos. Tabla observada:\n{}",
        render(&observed)
    );
    let observed_documents = documents_table(&destination);
    assert_eq!(
        observed_documents,
        expected(documents),
        "{name}: los documentos que cambian a propósito han cambiado. Tabla observada:\n{}",
        render(&observed_documents)
    );
}

/// **Muchos loci en dos secuencias, VCF, tartamudeo de base y quimeras: los
/// bytes son los congelados.**
#[test]
fn many_loci_with_base_stutter_and_chimeras_keep_their_frozen_bytes() {
    check_many_loci(
        "tartamudeo de base",
        MANY_LOCI_BASE_STUTTER,
        MANY_LOCI_BASE_STUTTER_DOCUMENTS,
        "",
    );
}

/// **Lo mismo con PCR y tartamudeo clonal por linaje: los bytes son los
/// congelados.**
#[test]
fn many_loci_with_pcr_clonal_stutter_and_chimeras_keep_their_frozen_bytes() {
    check_many_loci(
        "PCR con tartamudeo clonal",
        MANY_LOCI_PCR_CLONAL,
        MANY_LOCI_PCR_CLONAL_DOCUMENTS,
        PCR,
    );
}

/// **Cada escenario de la revisión de la fusión del VCF publica los bytes
/// congelados.**
#[test]
fn the_vcf_fusion_review_scenarios_keep_their_frozen_bytes() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/vcf-fusion/review");
    let mut names = fs::read_dir(&fixtures)
        .expect("escenarios de la revisión")
        .map(|entry| {
            entry
                .expect("entrada")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect::<Vec<_>>();
    names.sort();
    let mut observed = BTreeMap::new();
    let mut observed_documents = BTreeMap::new();
    let mut tables = String::new();
    for name in &names {
        let root = temporary();
        fs::create_dir_all(&root).expect("raíz del escenario");
        let _cleanup = TemporaryRoot(root.clone());
        for file in ["reference.fa", "snvs.tsv", "general.toml"] {
            fs::copy(fixtures.join(name).join(file), root.join(file))
                .expect("fichero del escenario");
        }
        let destination = run_binary(&root, "general.toml");
        let table = normalized_table(&destination);
        assert!(
            table.contains_key("outputs/truth.vcf"),
            "{name}: el escenario tiene que publicar el VCF"
        );
        let documents = documents_table(&destination);
        tables.push_str(&format!(
            "{name}:\n{}{}",
            render(&table),
            render(&documents)
        ));
        observed.insert(name.clone(), summary(&table));
        observed_documents.insert(name.clone(), summary(&documents));
    }
    assert_eq!(
        observed,
        expected(REVIEW_FIXTURES),
        "las tablas normalizadas de los escenarios de la revisión no son las de dcae406. \
         Resúmenes observados:\n{}\nTablas:\n{tables}",
        render(&observed)
    );
    assert_eq!(
        observed_documents,
        expected(REVIEW_FIXTURES_DOCUMENTS),
        "los documentos que cambian a propósito en los escenarios de la revisión han cambiado. \
         Resúmenes observados:\n{}\nTablas:\n{tables}",
        render(&observed_documents)
    );
}

/// **De dónde salen los literales normalizados.** Imprime, para cada árbol
/// publicado que haya bajo `GPS_NORMALIZE_TREES`, su tabla normalizada, el
/// resumen de esa tabla y la tabla de los tres documentos. Los literales de
/// `dcae406` se obtuvieron el 30/10/2026 lanzándola sobre los árboles que
/// publicó el binario de ese commit con los TOML de antes. No forma parte de
/// la batería: se lanza a mano con `--ignored`.
#[test]
#[ignore = "herramienta para obtener literales; se lanza a mano"]
fn print_normalized_tables_of_published_trees() {
    let Ok(trees) = std::env::var("GPS_NORMALIZE_TREES") else {
        return;
    };
    let mut names = fs::read_dir(&trees)
        .expect("árboles")
        .map(|entry| entry.expect("entrada").path())
        .collect::<Vec<_>>();
    names.sort();
    for tree in names {
        let table = normalized_table(&tree);
        let documents = documents_table(&tree);
        println!(
            "== {}\nresumen\t{}\ndocumentos\t{}\n{}--\n{}",
            tree.display(),
            summary(&table),
            summary(&documents),
            render(&table),
            render(&documents)
        );
    }
}
