//! Barrido de guardas: qué comprobaciones de los jueces independientes llegan
//! a dispararse durante la batería.
//!
//! ## Qué mide, y qué **no** mide
//!
//! Este módulo mide **una sola cosa**: si una comprobación de rechazo llega a
//! ejecutarse alguna vez mientras corre la batería. No mide si la prueba que
//! la dispara **fallaría** en caso de que la comprobación estuviera mal
//! escrita. Son cosas distintas y la segunda es la que importa: una prueba que
//! provoca un rechazo sin afirmar **cuál** rechazo es lo satisface igual que
//! una prueba que exige el motivo exacto.
//!
//! Dicho de otro modo: una cifra de cero guardas sin disparar **no demuestra
//! que los jueces estén probados**; demuestra que ninguna guarda está muerta.
//! Lo que demuestra que una guarda está probada es mutarla y ver morir la
//! prueba, y eso no lo puede hacer este barrido: lo hace quien escribe la
//! prueba, a mano, una mutación cada vez.
//!
//! Sirve, por tanto, como **suelo** y no como techo: impide que la cobertura
//! vuelva a caer sin que nadie se entere, que es exactamente lo que ocurrió
//! entre el cierre de F5 y septiembre de 2026.
//!
//! ## Cómo funciona
//!
//! Los tres jueces canalizan todos sus rechazos por dos o tres funciones
//! constructoras de error. Esas funciones llevan la anotación
//! [`track_caller`](https://doc.rust-lang.org/reference/attributes/codegen.html#the-track_caller-attribute),
//! que hace que la función sepa desde qué fichero y línea la han llamado, y
//! llaman a [`record`]. Con eso, una sola ejecución de la batería dice qué
//! guardas llegan a dispararse y cuáles no.
//!
//! El barrido está **apagado** salvo que la variable de entorno
//! `GPS_GUARD_TRACE_DIR` nombre un directorio. Apagado, cada rechazo cuesta la
//! lectura de un [`OnceLock`] ya inicializado; encendido, cuesta además un
//! bloqueo y —solo la primera vez que se ve cada línea— una escritura.
//!
//! Cada proceso escribe su propio fichero `<pid>.guards`, porque `cargo test`
//! lanza un proceso por binario de prueba y varios de las pruebas del proyecto
//! lanzan además subprocesos. El barrido junta después todos los ficheros del
//! directorio.
//!
//! ## Por qué la escritura vive aquí y no en cada juez
//!
//! Una prueba de `tests/f5_closure.rs` prohíbe que el texto de
//! `src/f5_closure.rs` contenga **cualquier vía de escritura en disco**, y ya
//! cazó la instrumentación de un revisor en cuanto la introdujo. Esa prueba
//! tiene razón y no se toca: la auditoría de identidad no escribe. Por eso el
//! único que abre un fichero es este módulo, y los jueces solo llaman a
//! [`record`].

use std::collections::BTreeSet;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::panic::Location;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

/// Variable de entorno que enciende el barrido y nombra su directorio.
pub const GUARD_TRACE_DIR_ENV: &str = "GPS_GUARD_TRACE_DIR";

/// Destino del barrido de este proceso.
///
/// `seen` evita reescribir una línea ya vista: tras el primer disparo, cada
/// rechazo cuesta una búsqueda en un conjunto ordenado y nada más.
struct Sink {
    file: File,
    seen: BTreeSet<(&'static str, u32)>,
}

/// Destino compartido, o `None` si el barrido está apagado.
///
/// Se resuelve una sola vez por proceso. Si la variable no está, si el
/// directorio no se puede crear o si el fichero no se puede abrir, el barrido
/// queda apagado en silencio: es instrumentación de medida, y no puede hacer
/// fallar una ejecución que sin ella habría ido bien.
fn sink() -> Option<&'static Mutex<Sink>> {
    static SINK: OnceLock<Option<Mutex<Sink>>> = OnceLock::new();
    SINK.get_or_init(|| {
        let directory = std::env::var_os(GUARD_TRACE_DIR_ENV)?;
        let directory = Path::new(&directory);
        std::fs::create_dir_all(directory).ok()?;
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(directory.join(format!("{}.guards", std::process::id())))
            .ok()?;
        Some(Mutex::new(Sink {
            file,
            seen: BTreeSet::new(),
        }))
    })
    .as_ref()
}

/// Anota que la guarda de quien llama se ha disparado.
///
/// La anotación `#[track_caller]` de esta función y la de las constructoras de
/// error de cada juez encadenan la posición hasta el **punto de rechazo real**,
/// no hasta la constructora. Sin ella, los ciento y pico rechazos que pasan por
/// `truth_reject` se anotarían todos en la misma línea.
#[track_caller]
pub fn record() {
    let Some(sink) = sink() else {
        return;
    };
    let location = Location::caller();
    let key = (location.file(), location.line());
    let Ok(mut sink) = sink.lock() else {
        return;
    };
    if sink.seen.insert(key) {
        // Un fallo de escritura no puede tumbar la ejecución medida.
        let _ = writeln!(sink.file, "{}:{}", key.0, key.1);
        let _ = sink.file.flush();
    }
}
