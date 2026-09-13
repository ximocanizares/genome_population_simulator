//! F5.5.1. Clausura de identidad: el paquete solo lee lo que F5.3 cerró.
//!
//! ## Qué se audita aquí, y por qué no es lo que audita F5.4.2
//!
//! El hito pide auditar que «el paquete publicado **solo lee** identidad,
//! copias, pertenencia y asignaciones **ya cerradas en F5.3**» y que «ningún
//! paso F5.5 puede crear, sustituir o modificar esos datos». Eso **no** es una
//! segunda validación de valores: F5.4.2 ya rederiva cada `FCP_`, `HAP_` y
//! `FHA_` desde los bytes con un oráculo independiente, y repetirlo aquí sería
//! el oráculo circular que esta base de código ya ha pagado siete veces.
//!
//! Lo que falta —y lo que este módulo produce— es la evidencia de una
//! propiedad **temporal** y de una **ausencia**:
//!
//! 1. **Precedencia, atestiguada.** Toda identidad que el paquete publica
//!    estaba cerrada en el repositorio de F5.3 **antes** de publicarse. Un
//!    fichero no puede atestiguar en qué orden se derivó lo que contiene, de
//!    modo que hace falta un testigo del propio acto de publicar. La primera
//!    versión de este hito usaba un **sello suelto** y la revisión adversarial
//!    demostró que no bastaba: sus cifras son una función determinista del
//!    repositorio sin consumir, así que un repositorio rematerializado
//!    **después** de que el paquete existiera, sellado entonces y consumido con
//!    un visitante vacío, avalaba el paquete viejo. Ahora la auditoría exige un
//!    [`F5ClosureWitness`], que solo nace de [`seal_and_publish_f5_success`] o
//!    [`seal_and_publish_f5_failure`] —sellan y **luego** publican— y que
//!    guarda el destino publicado, de modo que un testigo de otro paquete no
//!    acredita este.
//! 2. **Conservación.** El repositorio cuenta cada identidad derivada. Entre
//!    el sello —tomado **antes** de publicar— y la auditoría —tomada
//!    **después**— el registro solo puede haber crecido en exactamente una
//!    entrada por fila de `f5-founders.tsv`: las `FHA_` que la pasada única de
//!    [`F5FounderIdentities::for_each_assignment`] deriva. Ni una más, que
//!    sería una identidad creada después de F5.3; ni una menos, que sería una
//!    fila publicada sin derivar. En la mitad «fallo» la ley es la igualdad
//!    estricta: publicar un fallo **no deriva ninguna** identidad.
//! 3. **La clase fundadora también va sellada.** La Revisión F5.5 encontró la
//!    pieza que **ningún auditor cubría**: el censo **no está en el paquete**,
//!    de modo que el oráculo de F5.4.2 no puede verlo, y esta clausura lo
//!    recibía como cuarto parámetro suelto. Una segunda cadena que declarase
//!    `founder_ids` literales iguales a los `IND_` derivados **invertía** los
//!    dos recuentos que el hito publica como evidencia, y la auditoría decía
//!    `Ok`. Ahora la clase de cada subpoblación se congela en el testigo al
//!    sellar, tomada de la configuración que publica y solo si
//!    `matches_admission` la admite —`SealedF5Admission` cubre
//!    `subpopulation.founder_kind` en su material canónico—, y el censo **deja
//!    de ser un parámetro**: la sonda ya no es representable.
//! 4. **Imposibilidad de mutar.** La auditoría recibe `&F5FounderIdentities`,
//!    una referencia compartida. Derivar una identidad exige `&mut self`
//!    —`for_each_assignment` es el único camino— y el registro es privado del
//!    módulo fundador. Por tanto **ningún paso de F5.5 puede crear, sustituir
//!    ni modificar** una identidad: no es una promesa, es el comprobador de
//!    préstamos. Del lado del disco, la auditoría solo abre ficheros para
//!    leer y no escribe ninguno; la batería mide el árbol publicado antes y
//!    después y exige igualdad byte a byte.
//!
//! ## Qué no hace, y por qué
//!
//! - **No publica nada.** No añade un octavo artefacto —`artifacts` es un
//!   `prefixItems` de siete con `items:false`—, no toca el manifiesto y no
//!   reabre ningún esquema. Devuelve un informe en memoria; la evidencia vive
//!   en ese informe y en la batería que lo fija.
//! - **No rederiva identidades.** No abre BLAKE3 ni el catálogo de IDs. Si lo
//!   hiciera, estaría repitiendo F5.4.2 con menos independencia.
//! - **No juzga enlaces, tamaños ni SHA-256.** Eso es la puerta física de
//!   F5.4.3 y el oráculo de F5.4.2. Aquí un fichero alterado se detecta por
//!   la identidad que nombra, no por su huella.
//! - **No enumera el árbol.** La enumeración exhaustiva con igualdad exacta
//!   contra el conjunto declarado ya la hacen la puerta y el oráculo de fallo.
//!
//! ## Memoria
//!
//! La auditoría no retiene nada proporcional a `B*C_total`, al censo ni al
//! número de identidades: recorre `f5-founders.tsv` y `outputs/truth.jsonl`
//! por líneas, conserva una sola clave de copia y consulta las estructuras que
//! el repositorio **ya** sostenía antes de auditar. Sostiene **un registro
//! cada vez**: la línea que acaba de leer. Por eso los registros de verdad se
//! deserializan contra una estructura de cuatro campos y no contra un árbol
//! JSON —lo que descarta `segments`, `history` y `current_pieces` conforme se
//! leen—, y por eso el pico sigue al tamaño del mayor registro publicado, que
//! es el mismo término que el escritor ya sostuvo al producirlo.
//!
//! **No añade una tercera pieza fuera de `f5-resource-accounting-v9` junto al
//! oráculo y la publicación**: su pico cabe **dentro** del presupuesto de
//! memoria de trabajo que la admisión selló, de modo que no necesita ninguna
//! declaración fuera de la fórmula. La batería lo mide con el asignador
//! contador y lo acota además contra el mayor registro publicado: con `B=8`
//! pica 23 KiB frente a los 98 KiB sellados y con `B=32`, 69 KiB frente a
//! 188 KiB.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_artifacts::{F5_TRUTH_PATH, F5_TSV_CONTRACTS};
use crate::f5_founders::F5FounderIdentities;
use crate::f5_oracle::{F5_FAILURE_DOCUMENT_PATH, F5InitializationSummary, F5PackageDeclaration};
use crate::f5_preflight::{F5OtherResources, SealedF5Admission};
use crate::f5_publication::{
    F5FailureReport, F5PublicationFailure, F5PublicationRequest, publish_f5_failure,
    publish_f5_success,
};
use crate::identity::{ExecutionId, FounderCopyId};
use crate::publication::ExecutionWorkspace;

/// Posición de `f5-founders.tsv` dentro de los siete contratos TSV.
///
/// Una prueba de la batería fija que esa posición es la del contrato
/// `f5_founders`, de modo que un reordenamiento del arreglo no puede
/// convertir esta auditoría en la de otro fichero.
const FOUNDERS: usize = 5;

/// Columnas exactas de `f5-founders.tsv`, en su contrato v1.
const FOUNDER_COLUMNS: usize = 8;

/// Prefijo de clase que DTC-039 da a una copia fundadora persistente.
const FOUNDER_COPY_PREFIX: &str = "FCP_";
/// Prefijo de clase de una asignación por bloque.
const ASSIGNMENT_PREFIX: &str = "FHA_";
/// Prefijo de clase de una fila haplotípica.
const HAPLOTYPE_PREFIX: &str = "HAP_";

/// Resultado de una comprobación de clausura.
type Checked<T> = Result<T, F5ClosureRejection>;

// ---------------------------------------------------------------------------
// Rechazo tipado
// ---------------------------------------------------------------------------

/// Causa tipada de un rechazo de clausura.
///
/// Como en F5.4.2, no basta con que una corrupción falle: tiene que fallar
/// **por su causa**. La batería contrasta la causa de cada mutación en lugar
/// de limitarse a comprobar que hubo error.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum F5ClosureCause {
    /// El sello no procede del repositorio que se está auditando, o el
    /// repositorio ya no conserva lo que se selló.
    Seal,
    /// El registro de identidades no conserva el recuento: entre el sello y la
    /// auditoría se derivó algo distinto de una `FHA_` por fila publicada.
    Conservation,
    /// El paquete nombra una identidad que F5.3 no cerró.
    Unclosed,
    /// Un identificador publicado ocupa el lugar de otro de distinta clase.
    Substitution,
    /// La pertenencia publicada falta, se repite o no es la que se cerró.
    Membership,
    /// El individuo publicado no es el que el censo de F5.3 resolvió.
    Roster,
    /// Cabecera, columnas, orden o estructura de un fichero leído.
    ByteContract,
    /// Un fichero que el estado del paquete exige no puede abrirse.
    Missing,
}

impl F5ClosureCause {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Seal => "seal",
            Self::Conservation => "conservation",
            Self::Unclosed => "unclosed",
            Self::Substitution => "substitution",
            Self::Membership => "membership",
            Self::Roster => "roster",
            Self::ByteContract => "byte_contract",
            Self::Missing => "missing",
        }
    }
}

/// Rechazo de la auditoría de clausura: la causa tipada y el diagnóstico.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5ClosureRejection {
    cause: F5ClosureCause,
    diagnostic: Diagnostic,
}

impl F5ClosureRejection {
    #[must_use]
    pub const fn cause(&self) -> F5ClosureCause {
        self.cause
    }
    #[must_use]
    pub const fn diagnostic(&self) -> &Diagnostic {
        &self.diagnostic
    }
    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        self.diagnostic
    }
}

impl From<F5ClosureRejection> for Diagnostic {
    fn from(rejection: F5ClosureRejection) -> Self {
        rejection.diagnostic
    }
}

fn reject(
    cause: F5ClosureCause,
    path: &str,
    field: &str,
    detail: impl Into<String>,
) -> F5ClosureRejection {
    let diagnostic = Diagnostic::format(
        "GPS-FMT-101",
        "El paquete de perfiles F5 no satisface su formato, versión, orden o identidad declarada.",
    )
    .with("path", path)
    .with("field", field)
    .with("cause", cause.as_str())
    .with("detail", detail.into());
    F5ClosureRejection { cause, diagnostic }
}

// ---------------------------------------------------------------------------
// Sello
// ---------------------------------------------------------------------------

/// Las cuatro cifras del repositorio de identidades **antes** de publicar.
///
/// El sello por sí solo **no prueba nada**, y decirlo importa porque la
/// primera versión de este hito creyó lo contrario. Sus cuatro cifras son una
/// función determinista del repositorio sin consumir: cualquiera puede
/// rematerializar un repositorio equivalente desde la misma configuración
/// **después** de que el paquete exista, sellarlo, consumir su recorrido y
/// obtener exactamente las mismas cifras. La revisión adversarial lo
/// reprodujo y la auditoría aceptaba. Un sello suelto ya no es admisible: hay
/// que presentar un [`F5ClosureWitness`], y el único modo de obtener uno es
/// **publicar**.
///
/// Por eso el constructor es privado del módulo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5IdentitySeal {
    execution_id: String,
    copies: usize,
    units: usize,
    registered_identities: usize,
}

impl F5IdentitySeal {
    /// Sella el repositorio tal y como F5.3 lo dejó.
    fn seal(identities: &F5FounderIdentities) -> Self {
        Self {
            execution_id: identities.execution_id().as_str().to_owned(),
            copies: identities.copies().len(),
            units: identities.units().len(),
            registered_identities: identities.registered_identities(),
        }
    }

    /// La ejecución con la que se derivó cada identidad sellada.
    #[must_use]
    pub fn execution_id(&self) -> &str {
        &self.execution_id
    }
    /// Copias `FCP_` cerradas en F5.3.
    #[must_use]
    pub const fn copies(&self) -> usize {
        self.copies
    }
    /// Unidades bloque×subpoblación cerradas en F5.3.
    #[must_use]
    pub const fn units(&self) -> usize {
        self.units
    }
    /// Identidades registradas en el momento del sello: `FCP_` y `HAP_`, sin
    /// ninguna `FHA_`, que solo existe tras la pasada de publicación.
    #[must_use]
    pub const fn registered_identities(&self) -> usize {
        self.registered_identities
    }
}

// ---------------------------------------------------------------------------
// Testigo: lo que ata el sello a la publicación
// ---------------------------------------------------------------------------

/// Acredita que **este** sello se tomó antes de publicar **este** paquete.
///
/// Es la pieza que convierte la auditoría en una comprobación de precedencia.
/// Un testigo solo nace de [`seal_and_publish_f5_success`] o de
/// [`seal_and_publish_f5_failure`]: ambas sellan el repositorio y **después**
/// publican, de modo que existir ya demuestra el orden. Y como el testigo
/// guarda el destino, un testigo de otro paquete no acredita este: publicar un
/// repositorio rematerializado exige un destino nuevo —`ExecutionWorkspace`
/// rechaza un destino existente— y el testigo que devuelve nombra ese destino
/// nuevo, no el paquete que se quería avalar.
///
/// No tiene constructor público y sus campos son privados: no se puede
/// fabricar uno sin publicar.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5ClosureWitness {
    seal: F5IdentitySeal,
    destination: PathBuf,
    /// Clase fundadora de cada subpoblación, **congelada al sellar** y
    /// ordenada por `population_id`.
    ///
    /// Es la corrección de la Revisión F5.5: la clasificación
    /// aportado/derivado se leía del censo que la auditoría recibía como
    /// cuarto parámetro suelto, de modo que una segunda cadena que declarase
    /// `founder_ids` literales iguales a los `IND_` derivados **invertía** los
    /// dos recuentos que el plan nombra como evidencia del hito, y la
    /// auditoría devolvía `Ok`. Ningún auditor lo cubría: el censo **no está
    /// en el paquete**, de modo que el oráculo de F5.4.2 no puede verlo.
    ///
    /// Ahora la clase viaja en el testigo, tomada de la configuración que
    /// publica, y el censo deja de ser un parámetro: la sonda ya no es
    /// representable.
    ///
    /// `None` cuando esa configuración **no** es la que su admisión sella: el
    /// testigo existe pero no clasifica nada, y la auditoría lo rechaza.
    founder_kinds: Option<Vec<(String, bool)>>,
}

impl F5ClosureWitness {
    /// El sello que precedió a la publicación acreditada.
    #[must_use]
    pub const fn seal(&self) -> &F5IdentitySeal {
        &self.seal
    }
    /// El destino que este testigo acredita.
    #[must_use]
    pub fn destination(&self) -> &Path {
        &self.destination
    }

    /// El testigo tiene que acreditar el paquete que se está auditando, y ese
    /// paquete tiene que **existir**.
    ///
    /// Las dos rutas se canonicalizan y la canonicalización es **obligatoria**:
    /// la versión anterior caía a la comparación literal cuando alguna no
    /// existía, de modo que un destino borrado se comparaba igual a sí mismo y
    /// acreditaba un paquete inexistente. Ahora una ruta que no resuelve es
    /// `Missing`, que es lo que es.
    fn ensure_accredits(&self, payload_root: &Path) -> Checked<()> {
        let resolve = |path: &Path| -> Checked<PathBuf> {
            path.canonicalize().map_err(|error| {
                reject(
                    F5ClosureCause::Missing,
                    "witness",
                    "destination",
                    format!("{} no resuelve: {error}", path.display()),
                )
            })
        };
        let accredited = resolve(&self.destination)?;
        if accredited == resolve(payload_root)? {
            return Ok(());
        }
        Err(reject(
            F5ClosureCause::Seal,
            "witness",
            "destination",
            format!(
                "el testigo acredita {} y se audita {}",
                self.destination.display(),
                payload_root.display()
            ),
        ))
    }

    /// La clase fundadora sellada de una subpoblación: `true` si su identidad
    /// se derivó de una `founder_key`.
    fn founder_kind(&self, population_id: &str, path: &str) -> Checked<bool> {
        let kinds = self.founder_kinds.as_ref().ok_or_else(|| {
            reject(
                F5ClosureCause::Seal,
                "witness",
                "founder_kinds",
                "la configuración que publicó no es la que su admisión sella",
            )
        })?;
        kinds
            .binary_search_by(|(candidate, _)| candidate.as_str().cmp(population_id))
            .map(|position| kinds[position].1)
            .map_err(|_| {
                reject(
                    F5ClosureCause::Roster,
                    path,
                    "population_id",
                    format!("{population_id} no es ninguna de las subpoblaciones selladas"),
                )
            })
    }
}

/// Clase fundadora de cada subpoblación, tomada de la configuración que
/// publica y solo si su admisión sellada la admite.
///
/// `SealedF5Admission` cubre `subpopulation.founder_kind` en su material
/// canónico, de modo que `matches_admission` es la comprobación que ata esta
/// lectura al sello y no a un objeto en memoria cualquiera.
fn sealed_founder_kinds(
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
) -> Option<Vec<(String, bool)>> {
    if !sealed.matches_admission(config, F5OtherResources::default()) {
        return None;
    }
    let mut kinds = config
        .declared_scenario()
        .subpopulations
        .iter()
        .map(|subpopulation| {
            (
                subpopulation.population_id.value().clone(),
                subpopulation.founder_keys.is_some(),
            )
        })
        .collect::<Vec<_>>();
    kinds.sort_by(|left, right| left.0.cmp(&right.0));
    Some(kinds)
}

/// Sella el repositorio y publica el éxito, en ese orden.
///
/// Es la única vía que produce un [`F5ClosureWitness`] de éxito. No cambia
/// nada de F5.4.3: delega en [`publish_f5_success`] con sus mismos argumentos
/// y su mismo error, y se limita a tomar el sello **antes** y a atarlo al
/// destino que la publicación devuelve.
///
/// # Errors
///
/// El mismo [`F5PublicationFailure`] que [`publish_f5_success`].
pub fn seal_and_publish_f5_success<R>(
    workspace: ExecutionWorkspace,
    request: &F5PublicationRequest<'_>,
    identities: &mut F5FounderIdentities,
    recount: R,
) -> Result<F5ClosureWitness, Box<F5PublicationFailure>>
where
    R: FnOnce(&Path, &F5PackageDeclaration) -> Result<Vec<F5InitializationSummary>, Diagnostic>,
{
    let seal = F5IdentitySeal::seal(identities);
    let founder_kinds = sealed_founder_kinds(request.config(), request.sealed());
    let destination = publish_f5_success(workspace, request, identities, recount)?;
    Ok(F5ClosureWitness {
        seal,
        destination,
        founder_kinds,
    })
}

/// Sella el repositorio y publica el fallo controlado, en ese orden.
///
/// Es la única vía que produce un [`F5ClosureWitness`] de fallo. El sello se
/// toma **inmediatamente antes** de publicar, de modo que la igualdad que
/// [`audit_f5_failure_closure`] exige dice lo que tiene que decir: **publicar
/// un fallo no deriva ninguna identidad**, viniera el provisional con
/// artefactos ya escritos o sin ellos.
///
/// # Errors
///
/// El mismo diagnóstico que [`publish_f5_failure`].
pub fn seal_and_publish_f5_failure(
    workspace: ExecutionWorkspace,
    config: &EffectiveF5Config,
    sealed: &SealedF5Admission,
    execution_id: &ExecutionId,
    report: &F5FailureReport,
    identities: &F5FounderIdentities,
) -> Result<F5ClosureWitness, Diagnostic> {
    let seal = F5IdentitySeal::seal(identities);
    let founder_kinds = sealed_founder_kinds(config, sealed);
    let destination = publish_f5_failure(workspace, config, sealed, execution_id, report)?;
    Ok(F5ClosureWitness {
        seal,
        destination,
        founder_kinds,
    })
}

// ---------------------------------------------------------------------------
// Informe
// ---------------------------------------------------------------------------

/// Lo que la auditoría contó, para que el llamador pueda contrastarlo.
///
/// No es un artefacto: no se publica, no entra en el manifiesto y no tiene
/// esquema. Es la evidencia que el hito pide, en memoria.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct F5ClosureReport {
    copies: u64,
    assignments: u64,
    individuals: u64,
    supplied_individuals: u64,
    derived_individuals: u64,
}

impl F5ClosureReport {
    /// Copias `FCP_` distintas que el paquete publica y F5.3 había cerrado.
    #[must_use]
    pub const fn copies(&self) -> u64 {
        self.copies
    }
    /// Filas de `f5-founders.tsv`, que son las `FHA_` derivadas al publicar.
    #[must_use]
    pub const fn assignments(&self) -> u64 {
        self.assignments
    }
    /// Individuos distintos que el paquete publica.
    #[must_use]
    pub const fn individuals(&self) -> u64 {
        self.individuals
    }
    /// Individuos cuyo `IndividualId` lo aportó la usuaria por `founder_ids`.
    #[must_use]
    pub const fn supplied_individuals(&self) -> u64 {
        self.supplied_individuals
    }
    /// Individuos cuyo `IndividualId` se derivó de una `founder_key`.
    #[must_use]
    pub const fn derived_individuals(&self) -> u64 {
        self.derived_individuals
    }
}

// ---------------------------------------------------------------------------
// Auditoría de un paquete de éxito
// ---------------------------------------------------------------------------

/// Audita que un paquete **publicado con éxito** solo nombra identidad,
/// copias, pertenencia y asignaciones ya cerradas en F5.3.
///
/// `payload_root` es la raíz bajo la que viven `outputs/...` del destino
/// **ya renombrado**. `seal` tiene que haberse tomado antes de publicar.
///
/// Comprueba, en este orden:
///
/// 1. que el sello sea el del repositorio auditado y que este no haya perdido
///    nada de lo sellado;
/// 2. que cada fila de `f5-founders.tsv` nombre una copia, un individuo, una
///    subpoblación y una fila haplotípica que el repositorio ya contenía, sin
///    sustituir una clase de identificador por otra;
/// 3. que las copias publicadas sean **exactamente** las cerradas, cada una
///    una sola vez;
/// 4. que el registro de identidades haya crecido en exactamente una `FHA_`
///    por fila publicada;
/// 5. que la pertenencia de `outputs/truth.jsonl` sea la misma, completa y sin
///    repeticiones.
///
/// No escribe ningún byte: recibe el repositorio por referencia compartida, de
/// modo que el comprobador de préstamos impide derivar ninguna identidad, y
/// abre los dos ficheros solo para leer.
///
/// # Errors
///
/// [`F5ClosureRejection`], cuya [`F5ClosureCause`] nombra el contrato roto.
pub fn audit_f5_identity_closure(
    payload_root: &Path,
    witness: &F5ClosureWitness,
    identities: &F5FounderIdentities,
) -> Checked<F5ClosureReport> {
    witness.ensure_accredits(payload_root)?;
    let seal = &witness.seal;
    check_seal(seal, identities)?;
    check_repository_order(identities)?;
    let report = walk_founders(payload_root, witness, identities)?;
    check_conservation(seal, identities, report.assignments)?;
    walk_truth_membership(payload_root, identities, report.copies)?;
    Ok(report)
}

/// Audita que un paquete **fallido** no exponga ninguna identidad fundadora.
///
/// DTC-045 permite publicar un fallo controlado; DTC-089 le prohíbe llevar
/// artefactos biológicos. La mitad de este hito que corresponde al fallo es
/// exactamente esa: un fallo no crea, no sustituye y no modifica identidad
/// porque **no publica ninguna**. Se comprueba que ninguno de los ocho
/// ficheros biológicos exista bajo la carga útil y que el repositorio conserve
/// lo que selló.
///
/// La enumeración exhaustiva del árbol y la igualdad con el conjunto declarado
/// son de la puerta de F5.4.3 y de `audit_f5_failure_package`; aquí no se
/// repiten.
///
/// # Errors
///
/// [`F5ClosureRejection`], cuya [`F5ClosureCause`] nombra el contrato roto.
pub fn audit_f5_failure_closure(
    payload_root: &Path,
    witness: &F5ClosureWitness,
    identities: &F5FounderIdentities,
) -> Checked<F5ClosureReport> {
    witness.ensure_accredits(payload_root)?;
    // **Una afirmación positiva sobre el disco.** Sin ella la mitad «fallo»
    // era una comprobación construida entera sobre ausencias —el sello mira
    // memoria y los ocho ficheros biológicos se exigen ausentes—, y un destino
    // borrado entero, o sin su documento de fallo, la satisfacía. Se abre el
    // diagnóstico y no se lee: juzgar su contenido es de
    // `audit_f5_failure_package`, no de aquí.
    drop(open_for_reading(payload_root, F5_FAILURE_DOCUMENT_PATH)?);
    let seal = &witness.seal;
    check_seal(seal, identities)?;
    // La ley de conservación de la mitad «fallo»: publicar un fallo **no
    // deriva ninguna identidad**. `check_seal` no la cubría —solo comprobaba
    // que el registro no encogiese, y el registro nunca encoge—, de modo que
    // la mitad declarada del hito no se aplicaba.
    if identities.registered_identities() != seal.registered_identities {
        return Err(reject(
            F5ClosureCause::Conservation,
            F5_FAILURE_DOCUMENT_PATH,
            "registered_identities",
            format!(
                "el sello registró {} identidades y el repositorio registra {}: publicar un fallo no deriva ninguna",
                seal.registered_identities,
                identities.registered_identities()
            ),
        ));
    }
    for path in biological_paths() {
        // `symlink_metadata` no sigue enlaces: un enlace colgado hacia un
        // fichero biológico sigue siendo un fichero biológico publicado.
        if payload_root.join(path).symlink_metadata().is_ok() {
            return Err(reject(
                F5ClosureCause::Unclosed,
                path,
                "path",
                "un paquete fallido publica un artefacto biológico y expondría identidad fundadora",
            ));
        }
    }
    Ok(F5ClosureReport::default())
}

/// Los ocho ficheros que un éxito publica y un fallo no puede publicar.
fn biological_paths() -> impl Iterator<Item = &'static str> {
    F5_TSV_CONTRACTS
        .iter()
        .map(|contract| contract.path())
        .chain(std::iter::once(F5_TRUTH_PATH))
}

// ---------------------------------------------------------------------------
// Sello y conservación
// ---------------------------------------------------------------------------

fn check_seal(seal: &F5IdentitySeal, identities: &F5FounderIdentities) -> Checked<()> {
    if seal.execution_id != identities.execution_id().as_str() {
        return Err(reject(
            F5ClosureCause::Seal,
            "identities",
            "execution_id",
            format!(
                "el sello es de la ejecución {} y el repositorio es de {}",
                seal.execution_id,
                identities.execution_id().as_str()
            ),
        ));
    }
    if seal.copies != identities.copies().len() || seal.units != identities.units().len() {
        return Err(reject(
            F5ClosureCause::Seal,
            "identities",
            "copies,units",
            format!(
                "el sello conserva {} copias y {} unidades y el repositorio tiene {} y {}",
                seal.copies,
                seal.units,
                identities.copies().len(),
                identities.units().len()
            ),
        ));
    }
    // Aquí **no** se comprueba que el registro no haya encogido. Esa rama era
    // inalcanzable —el registro de F5.3 solo crece y el sello no tiene otro
    // constructor— y la desigualdad que sí importa la cubre la ley de
    // conservación, que exige igualdad exacta en las dos mitades.
    Ok(())
}

/// El registro solo puede haber crecido en una `FHA_` por fila publicada.
fn check_conservation(
    seal: &F5IdentitySeal,
    identities: &F5FounderIdentities,
    assignments: u64,
) -> Checked<()> {
    let derived = usize::try_from(assignments).ok().ok_or_else(|| {
        reject(
            F5ClosureCause::Conservation,
            F5_TSV_CONTRACTS[FOUNDERS].path(),
            "rows",
            "el número de filas publicadas no cabe en la máquina",
        )
    })?;
    let expected = seal
        .registered_identities
        .checked_add(derived)
        .ok_or_else(|| {
            reject(
                F5ClosureCause::Conservation,
                F5_TSV_CONTRACTS[FOUNDERS].path(),
                "rows",
                "el recuento de identidades esperado desborda la máquina",
            )
        })?;
    if identities.registered_identities() != expected {
        return Err(reject(
            F5ClosureCause::Conservation,
            F5_TSV_CONTRACTS[FOUNDERS].path(),
            "registered_identities",
            format!(
                "el sello registró {} identidades, el paquete publica {assignments} asignaciones \
                 y el repositorio registra {}: la publicación no derivó exactamente una FHA_ por fila",
                seal.registered_identities,
                identities.registered_identities()
            ),
        ));
    }
    Ok(())
}

/// La búsqueda binaria de unidades exige su orden canónico.
///
/// Es una **precondición de la propia auditoría**, no un hallazgo sobre el
/// paquete: el repositorio declara ese orden en su contrato y ninguna
/// configuración admitida lo rompe —el estudio de mutación de la batería lo
/// confirma: desactivar esta función no tumba ninguna prueba—. Se conserva
/// porque lo que protege es un **falso negativo de esta auditoría**: una
/// búsqueda binaria sobre datos desordenados devolvería «ausente» sobre
/// identidades presentes y la clausura rechazaría un paquete correcto citando
/// una causa falsa. Cuesta un recorrido lineal y no retiene nada.
///
/// La mitad que validaba el orden del censo desapareció con el censo: la
/// clase fundadora vive ahora en el testigo, ordenada al sellar.
fn check_repository_order(identities: &F5FounderIdentities) -> Checked<()> {
    if identities.units().windows(2).any(|pair| {
        (
            pair[0].population_id().as_str(),
            pair[0].block_id().as_str(),
        ) >= (
            pair[1].population_id().as_str(),
            pair[1].block_id().as_str(),
        )
    }) {
        return Err(reject(
            F5ClosureCause::Seal,
            "identities",
            "units",
            "las unidades cerradas no llegan ordenadas por población y bloque",
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// `f5-founders.tsv`
// ---------------------------------------------------------------------------

/// Clave de copia publicada, la única pieza que el recorrido retiene.
///
/// Se construye **una vez por copia**, no por fila: el recorrido compara la
/// fila en préstamo contra la clave viva y solo asigna memoria cuando la copia
/// cambia, de modo que las `B*C_total` filas no pagan cuatro cadenas cada una.
#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
struct CopyKey {
    population_id: String,
    individual_id: String,
    sequence_id: String,
    copy_index: u64,
}

impl CopyKey {
    /// Ordena la clave viva frente a los campos prestados de una fila.
    fn cmp_row(
        &self,
        population_id: &str,
        individual_id: &str,
        sequence_id: &str,
        copy_index: u64,
    ) -> std::cmp::Ordering {
        (
            self.population_id.as_str(),
            self.individual_id.as_str(),
            self.sequence_id.as_str(),
            self.copy_index,
        )
            .cmp(&(population_id, individual_id, sequence_id, copy_index))
    }

    /// `true` si la fila prestada pertenece a otro individuo.
    fn other_individual(&self, population_id: &str, individual_id: &str) -> bool {
        self.population_id != population_id || self.individual_id != individual_id
    }
}

// El recorrido es largo porque es **una sola pasada**: partir el cuerpo en
// funciones obligaría a devolver la fila prestada, la clave viva y el informe
// a través de una estructura auxiliar, y el préstamo de `line` no sobrevive a
// esa frontera sin asignar por fila. Se retira en cuanto `f5-founders.tsv`
// deje de comprobarse en una pasada única o el cuerpo pierda alguna de sus
// seis comprobaciones.
#[allow(clippy::too_many_lines)]
fn walk_founders(
    payload_root: &Path,
    witness: &F5ClosureWitness,
    identities: &F5FounderIdentities,
) -> Checked<F5ClosureReport> {
    let contract = &F5_TSV_CONTRACTS[FOUNDERS];
    let path = contract.path();
    let mut reader = open_for_reading(payload_root, path)?;
    let mut line = String::new();
    read_line(&mut reader, &mut line, path)?;
    // La cabecera exacta es también la mitad byte a byte de «ausencia de
    // migración, selección o genealogía»: el fichero que sostiene la identidad
    // fundadora no admite una columna de parental, de origen migratorio ni de
    // eficacia biológica, porque no admite ninguna columna más.
    if line.strip_suffix('\n') != Some(contract.header()) {
        return Err(reject(
            F5ClosureCause::ByteContract,
            path,
            "header",
            "la cabecera publicada no es exactamente la de su contrato",
        ));
    }

    let mut report = F5ClosureReport::default();
    let mut previous: Option<CopyKey> = None;
    let mut previous_block = String::new();
    loop {
        if read_line(&mut reader, &mut line, path)? == 0 {
            break;
        }
        let fields = split_exact(strip_lf(&line, path)?, path)?;
        let copy_index = fields[3].parse::<u64>().map_err(|_| {
            reject(
                F5ClosureCause::ByteContract,
                path,
                "copy_index",
                format!("{} no es un entero", fields[3]),
            )
        })?;
        check_class(fields[5], FOUNDER_COPY_PREFIX, path, "founder_copy_id")?;
        check_class(fields[6], ASSIGNMENT_PREFIX, path, "block_assignment_id")?;
        check_class(fields[7], HAPLOTYPE_PREFIX, path, "haplotype_id")?;

        let founder_copy_id = FounderCopyId::parse(fields[5], "f5-founders.tsv").map_err(|_| {
            reject(
                F5ClosureCause::ByteContract,
                path,
                "founder_copy_id",
                format!("{} no es una identidad tipada", fields[5]),
            )
        })?;
        let closed = identities.copy(&founder_copy_id).ok_or_else(|| {
            reject(
                F5ClosureCause::Unclosed,
                path,
                "founder_copy_id",
                format!("{} no es ninguna de las copias que F5.3 cerró", fields[5]),
            )
        })?;
        if closed.population_id().as_str() != fields[0] {
            return Err(reject(
                F5ClosureCause::Membership,
                path,
                "population_id",
                format!(
                    "la copia {} pertenece a {} y el paquete la publica en {}",
                    fields[5],
                    closed.population_id().as_str(),
                    fields[0]
                ),
            ));
        }
        if closed.individual_id().as_str() != fields[1] {
            return Err(reject(
                F5ClosureCause::Roster,
                path,
                "individual_id",
                format!(
                    "la copia {} es del fundador {} y el paquete la publica como {}",
                    fields[5],
                    closed.individual_id().as_str(),
                    fields[1]
                ),
            ));
        }
        if closed.sequence_id().as_str() != fields[2]
            || u64::from(closed.copy_index()) != copy_index
        {
            return Err(reject(
                F5ClosureCause::Unclosed,
                path,
                "sequence_id,copy_index",
                format!(
                    "la copia {} se cerró sobre {}/{} y el paquete la publica sobre {}/{}",
                    fields[5],
                    closed.sequence_id().as_str(),
                    closed.copy_index(),
                    fields[2],
                    copy_index
                ),
            ));
        }

        // El orden se juzga **después** de las comparaciones anteriores. Si se
        // juzgara antes, una identidad sustituida tumbaría el fichero por
        // «orden» en vez de por la sustitución, y la causa publicada sería
        // falsa: la prueba pasaría creyendo medir lo que no mide.
        match previous
            .as_ref()
            .map(|last| last.cmp_row(fields[0], fields[1], fields[2], copy_index))
        {
            // Misma copia: los bloques tienen que avanzar dentro de ella.
            Some(std::cmp::Ordering::Equal) => {
                if previous_block.as_str() >= fields[4] {
                    return Err(reject(
                        F5ClosureCause::ByteContract,
                        path,
                        "block_id",
                        "los bloques de una copia no avanzan",
                    ));
                }
            }
            // Copia nueva. Si además cambia el individuo, el censo cerrado dice
            // si su identidad la aportó la usuaria o se derivó de una clave.
            None | Some(std::cmp::Ordering::Less) => {
                report.copies += 1;
                if previous
                    .as_ref()
                    .is_none_or(|last| last.other_individual(fields[0], fields[1]))
                {
                    report.individuals += 1;
                    if witness.founder_kind(fields[0], path)? {
                        report.derived_individuals += 1;
                    } else {
                        report.supplied_individuals += 1;
                    }
                }
                previous = Some(CopyKey {
                    population_id: fields[0].to_owned(),
                    individual_id: fields[1].to_owned(),
                    sequence_id: fields[2].to_owned(),
                    copy_index,
                });
            }
            // Una copia ya cerrada vuelve a aparecer: su pertenencia se
            // publica dos veces, o el fichero no está en su clave.
            Some(std::cmp::Ordering::Greater) => {
                return Err(reject(
                    F5ClosureCause::Membership,
                    path,
                    "founder_copy_id",
                    format!(
                        "la copia {}/{}/{}/{copy_index} se publica fuera de su orden o más de una vez",
                        fields[0], fields[1], fields[2]
                    ),
                ));
            }
        }
        // El búfer se reutiliza: tras la primera fila no vuelve a asignar.
        previous_block.clear();
        previous_block.push_str(fields[4]);

        // La fila haplotípica asignada tiene que ser una de las que F5.3 cerró
        // para ese bloque×subpoblación.
        let unit = identities
            .units()
            .binary_search_by(|candidate| {
                (
                    candidate.population_id().as_str(),
                    candidate.block_id().as_str(),
                )
                    .cmp(&(fields[0], fields[4]))
            })
            .map(|position| &identities.units()[position])
            .map_err(|_| {
                reject(
                    F5ClosureCause::Unclosed,
                    path,
                    "block_id",
                    format!(
                        "{}×{} no es ninguna de las unidades que F5.3 cerró",
                        fields[0], fields[4]
                    ),
                )
            })?;
        if !unit
            .haplotype_ids()
            .iter()
            .any(|candidate| candidate.as_str() == fields[7])
        {
            return Err(reject(
                F5ClosureCause::Unclosed,
                path,
                "haplotype_id",
                format!(
                    "{} no es ninguna de las filas haplotípicas que F5.3 cerró para {}×{}",
                    fields[7], fields[0], fields[4]
                ),
            ));
        }

        report.assignments += 1;
    }

    let closed_copies = u64::try_from(identities.copies().len()).unwrap_or(u64::MAX);
    if report.copies != closed_copies {
        return Err(reject(
            F5ClosureCause::Membership,
            path,
            "founder_copy_id",
            format!(
                "F5.3 cerró {closed_copies} copias y el paquete publica {}",
                report.copies
            ),
        ));
    }
    // Aquí **no** se compara `report.individuals` con la suma de
    // `individual_count`. Un estudio de mutación sobre esta implementación
    // demostró que esa comparación no puede fallar: cada copia publicada ya
    // probó su individuo contra la copia cerrada, el conjunto de copias
    // publicadas ya se probó igual al cerrado, y toda copia cerrada lleva el
    // individuo que el censo resolvió. Una aserción que no puede fallar no es
    // una comprobación; se retiró en vez de decorar el módulo con ella.
    Ok(report)
}

/// Cada identificador ocupa el lugar de su propia clase, no el de otra.
fn check_class(value: &str, prefix: &str, path: &str, field: &str) -> Checked<()> {
    if value.starts_with(prefix) {
        return Ok(());
    }
    Err(reject(
        F5ClosureCause::Substitution,
        path,
        field,
        format!("{value} no es un identificador de clase {prefix}"),
    ))
}

// ---------------------------------------------------------------------------
// `outputs/truth.jsonl`
// ---------------------------------------------------------------------------

/// La pertenencia publicada es la cerrada, completa y sin repeticiones.
///
/// `truth.jsonl` es el único artefacto que publica `founder_membership`, de
/// modo que es el único sitio donde «pertenencia ausente o duplicada» puede
/// comprobarse sobre bytes.
/// Los dos recorridos exigen LF **en cada línea que leen**, como el oráculo.
///
/// La afirmación es esa y no más, y la Revisión F5.5 corrigió la anterior, que
/// decía «LF final» y era falsa: `walk_truth_membership` rompe en el primer
/// genotipo y la última línea de `truth.jsonl` es siempre un genotipo, de modo
/// que el LF **final** del documento nunca se mira. Un documento sin él se
/// acepta aquí; uno truncado dentro de la sección de copias sí se rechaza.
/// Cerrar el recorrido costaría una pasada entera por un contrato que la
/// puerta física de F5.4.3 y el oráculo de F5.4.2 ya sostienen sobre el
/// fichero completo, de modo que se rebaja la afirmación en vez de pagarla.
fn strip_lf<'a>(line: &'a str, path: &str) -> Checked<&'a str> {
    line.strip_suffix('\n').ok_or_else(|| {
        reject(
            F5ClosureCause::ByteContract,
            path,
            "eol",
            "falta el LF final",
        )
    })
}

/// Las dos claves que distinguen la cabecera v2 de un registro.
#[derive(Deserialize)]
struct TruthHeaderLine {
    record_type: String,
    /// Opcional a propósito: un registro de copia **no** la lleva, y lo que
    /// tiene que rechazarlo es la comprobación explícita de abajo, no un error
    /// de deserialización que publicaría una causa menos precisa.
    document_kind: Option<String>,
}

#[derive(Deserialize)]
struct TruthRecord {
    record_type: String,
    /// Presentes solo en `haplotype_copy`; los genotipos no los llevan.
    identity: Option<RecordIdentity>,
    founder_membership: Option<RecordMembership>,
}

#[derive(Deserialize)]
struct RecordIdentity {
    individual_id: String,
    sequence_id: String,
    copy_index: u64,
}

#[derive(Deserialize)]
struct RecordMembership {
    founder_copy_id: String,
    population_id: String,
}

fn walk_truth_membership(
    payload_root: &Path,
    identities: &F5FounderIdentities,
    published_copies: u64,
) -> Checked<()> {
    let path = F5_TRUTH_PATH;
    let mut reader = open_for_reading(payload_root, path)?;
    let mut line = String::new();
    if read_line(&mut reader, &mut line, path)? == 0 {
        return Err(reject(
            F5ClosureCause::ByteContract,
            path,
            "header",
            "documento vacío",
        ));
    }
    // La primera línea tiene que **ser** la cabecera, no solo existir. Si
    // faltara, la primera copia se consumiría en su lugar y el rechazo saldría
    // con una causa falsa: el mismo defecto que el otro recorrido ya corrigió
    // al mover el juicio de orden por detrás de las comparaciones.
    let header: TruthHeaderLine =
        serde_json::from_str(strip_lf(&line, path)?).map_err(|error| {
            reject(
                F5ClosureCause::ByteContract,
                path,
                "header",
                format!("cabecera ilegible: {error}"),
            )
        })?;
    if header.record_type != "header" || header.document_kind.as_deref() != Some("gps_truth") {
        return Err(reject(
            F5ClosureCause::ByteContract,
            path,
            "header",
            "la primera línea no es la cabecera v2 del documento de verdad",
        ));
    }
    let mut seen = 0_u64;
    let mut previous: Option<CopyKey> = None;
    loop {
        if read_line(&mut reader, &mut line, path)? == 0 {
            break;
        }
        let body = strip_lf(&line, path)?;
        // Se deserializa contra una estructura de cuatro campos, no contra un
        // árbol JSON: los `segments`, `history` y `current_pieces` del registro
        // se descartan conforme se leen y no llegan a materializarse. Es lo que
        // mantiene el pico dentro del presupuesto sellado.
        let record: TruthRecord = serde_json::from_str(body).map_err(|error| {
            reject(
                F5ClosureCause::ByteContract,
                path,
                "row",
                format!("registro ilegible: {error}"),
            )
        })?;
        if record.record_type != "haplotype_copy" {
            // La sección de copias termina en el primer genotipo.
            break;
        }
        let membership = record.founder_membership.ok_or_else(|| {
            reject(
                F5ClosureCause::Membership,
                path,
                "founder_membership",
                "una copia publicada no declara su pertenencia",
            )
        })?;
        let identity = record.identity.ok_or_else(|| {
            reject(
                F5ClosureCause::ByteContract,
                path,
                "identity",
                "una copia publicada no declara su identidad",
            )
        })?;
        let published_copy_id = membership.founder_copy_id;
        let key = CopyKey {
            population_id: membership.population_id,
            individual_id: identity.individual_id,
            sequence_id: identity.sequence_id,
            copy_index: identity.copy_index,
        };
        check_class(
            &published_copy_id,
            FOUNDER_COPY_PREFIX,
            path,
            "founder_copy_id",
        )?;
        let founder_copy_id = FounderCopyId::parse(published_copy_id.clone(), "truth.jsonl")
            .map_err(|_| {
                reject(
                    F5ClosureCause::ByteContract,
                    path,
                    "founder_copy_id",
                    format!("{published_copy_id} no es una identidad tipada"),
                )
            })?;
        let closed = identities.copy(&founder_copy_id).ok_or_else(|| {
            reject(
                F5ClosureCause::Unclosed,
                path,
                "founder_copy_id",
                format!("{published_copy_id} no es ninguna de las copias que F5.3 cerró"),
            )
        })?;
        if closed.population_id().as_str() != key.population_id {
            return Err(reject(
                F5ClosureCause::Membership,
                path,
                "founder_membership.population_id",
                format!(
                    "la copia {published_copy_id} pertenece a {} y la verdad la publica en {}",
                    closed.population_id().as_str(),
                    key.population_id
                ),
            ));
        }
        if closed.individual_id().as_str() != key.individual_id
            || closed.sequence_id().as_str() != key.sequence_id
            || u64::from(closed.copy_index()) != key.copy_index
        {
            return Err(reject(
                F5ClosureCause::Unclosed,
                path,
                "identity",
                format!(
                    "la copia {published_copy_id} se cerró sobre {}/{}/{} y la verdad la publica sobre {}/{}/{}",
                    closed.individual_id().as_str(),
                    closed.sequence_id().as_str(),
                    closed.copy_index(),
                    key.individual_id,
                    key.sequence_id,
                    key.copy_index
                ),
            ));
        }
        // Igual que en el TSV, el orden se juzga después de las
        // comparaciones de identidad, para que la causa publicada sea la real.
        if previous.as_ref().is_some_and(|last| *last >= key) {
            return Err(reject(
                F5ClosureCause::Membership,
                path,
                "founder_membership",
                format!(
                    "la pertenencia de {published_copy_id} se publica fuera de orden o más de una vez"
                ),
            ));
        }
        previous = Some(key);
        seen += 1;
    }
    if seen != published_copies {
        return Err(reject(
            F5ClosureCause::Membership,
            path,
            "founder_membership",
            format!(
                "f5-founders.tsv asigna {published_copies} copias y la verdad publica la pertenencia de {seen}"
            ),
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Lectura
// ---------------------------------------------------------------------------

fn open_for_reading(root: &Path, path: &str) -> Checked<BufReader<File>> {
    let full = root.join(path);
    File::open(&full).map(BufReader::new).map_err(|error| {
        reject(
            F5ClosureCause::Missing,
            path,
            "path",
            format!("{}: {error}", full.display()),
        )
    })
}

fn read_line(reader: &mut BufReader<File>, line: &mut String, path: &str) -> Checked<usize> {
    line.clear();
    reader.read_line(line).map_err(|error| {
        reject(
            F5ClosureCause::ByteContract,
            path,
            "encoding",
            error.to_string(),
        )
    })
}

/// Las ocho columnas de una fila fundadora, sin asignar ningún vector.
fn split_exact<'a>(body: &'a str, path: &str) -> Checked<[&'a str; FOUNDER_COLUMNS]> {
    let mut fields = [""; FOUNDER_COLUMNS];
    let mut seen = 0_usize;
    for field in body.split('\t') {
        if seen == FOUNDER_COLUMNS {
            return Err(reject(
                F5ClosureCause::ByteContract,
                path,
                "columns",
                format!("la fila tiene más de {FOUNDER_COLUMNS} columnas"),
            ));
        }
        fields[seen] = field;
        seen += 1;
    }
    if seen != FOUNDER_COLUMNS {
        return Err(reject(
            F5ClosureCause::ByteContract,
            path,
            "columns",
            format!("la fila tiene {seen} columnas y el contrato fija {FOUNDER_COLUMNS}"),
        ));
    }
    Ok(fields)
}
