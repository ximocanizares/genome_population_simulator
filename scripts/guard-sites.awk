# Enumera los **puntos de rechazo** de un juez independiente.
#
# Un punto de rechazo es una llamada a una de las funciones constructoras de
# error del módulo. El barrido de `src/guard_trace.rs` anota la línea de esa
# llamada, de modo que esta enumeración y aquella medición hablan de lo mismo.
#
# Los rechazos de **texto suelto** —un validador que devuelve `Err("texto")` y
# deja que el llamante lo envuelva en un diagnóstico— se cuentan por su propia
# llamada a `reason(`, que anota la línea del rechazo (DT-001); el envoltorio
# que lo convierte en diagnóstico no es un punto de rechazo y no llama a una
# constructora anotada.
#
# Desde el 05/10/2026 también se cuentan las llamadas a `rejected(`, la
# constructora genérica de `src/guard_trace.rs` que envuelve el error de cada
# rechazo (`Err(rejected(Variante { .. }))`), y las definiciones de las
# constructoras `pub(crate) fn` o `pub fn` se descartan igual que las de
# columna cero con `fn` a secas.
#
# Se descartan:
#   - las **definiciones** de las constructoras y sus cuerpos, porque la
#     anotación `#[track_caller]` los atraviesa y anota al llamante;
#   - los módulos `#[cfg(test)]`, que no forman parte del juez;
#   - los comentarios, incluidos los de documentación.
#
# Uso: awk -v names="reject|resource_reject" -f guard-sites.awk fichero.rs
# Salida: una línea por punto de rechazo, con el número de línea.

BEGIN {
    in_test = 0
    in_constructor = 0
}

# Un módulo de pruebas empieza con `#[cfg(test)]` en la columna 0 y termina en
# la primera llave de cierre en la columna 0.
/^#\[cfg\(test\)\]/ { in_test = 1; next }
in_test { if ($0 == "}") { in_test = 0 } ; next }

# El cuerpo de una constructora empieza en su `fn` de columna 0 y termina
# igual. `#[track_caller]` hace que sus rechazos se anoten en el llamante, así
# que ninguna línea suya es un punto de rechazo propio.
{
    if (match($0, "^(pub(\\([a-z]+\\))? )?fn (" names ")[(<]")) { in_constructor = 1; next }
}
in_constructor { if ($0 == "}") { in_constructor = 0 } ; next }

{
    line = $0
    # Fuera comentarios: `///`, `//!` y `//` sueltos.
    index_comment = index(line, "//")
    if (index_comment > 0) { line = substr(line, 1, index_comment - 1) }
    # Una llamada es el nombre precedido de algo que no es un carácter de
    # identificador, para que `truth_reject(` no cuente como `reject(`.
    if (match(line, "(^|[^A-Za-z0-9_.])(" names ")\\(")) {
        printf "%d\n", FNR
    }
}
