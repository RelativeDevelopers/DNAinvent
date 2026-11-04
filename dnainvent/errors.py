"""Errores de DNAinvent.

Los errores no llevan frases, sino un código y unos datos. Así cada interfaz
los redacta en el idioma del usuario:

    raise InvalidSequenceError("invalid_characters", kind="dna", chars=["X"], allowed=["A", "C", "G", "T"])

    es.json -> "Caracteres no válidos para ADN: X (permitidos: A, C, G, T)."
    en.json -> "Invalid characters for DNA: X (allowed: A, C, G, T)."

Lista de códigos (cada uno debe tener su texto en todos los idiomas):
    empty_sequence              kind
    invalid_characters          kind, chars, allowed
    odd_bit_length              length, multiple
    cut_out_of_range            cut
    guide_too_short             minimum
    invalid_bit_count           (sin datos)
    negative_permutation_size   (sin datos)
    start_out_of_range          start, length
    count_too_small             minimum
    no_orf                      (sin datos)
    min_length_too_small        minimum
    invalid_accession           accession
    protein_not_found           accession
    service_unavailable         service
    invalid_email               (sin datos)
    search_failed               service
    search_timeout              service, minutes
    prediction_too_long         length, maximum
    prediction_failed           service
    gene_not_found              symbol, taxon
    sequence_not_found          accession
    variant_not_found           variant
    variant_not_supported       hgvs
    variant_outside_cds         hgvs, start, end
    variant_mismatch            hgvs, accession, position, found, expected
"""

ERROR_CODES = frozenset({
    "empty_sequence", "invalid_characters", "odd_bit_length", "cut_out_of_range",
    "guide_too_short", "invalid_bit_count", "negative_permutation_size",
    "start_out_of_range", "count_too_small", "no_orf", "min_length_too_small",
    "invalid_accession", "protein_not_found", "service_unavailable",
    "invalid_email", "search_failed", "search_timeout", "prediction_too_long",
    "prediction_failed", "gene_not_found", "sequence_not_found", "variant_not_found",
    "variant_not_supported", "variant_outside_cds", "variant_mismatch",
})


class DNAinventError(ValueError):
    """Error base: un código legible por máquinas y los datos para redactarlo."""

    def __init__(self, code: str, **params):
        if code not in ERROR_CODES:
            raise ValueError(f"Código de error desconocido: {code}")
        self.code = code
        self.params = params
        super().__init__(code)

    def __str__(self) -> str:
        # Texto técnico, solo para programadores (logs, tests). La interfaz usa el traductor.
        details = ", ".join(f"{key}={value!r}" for key, value in self.params.items())
        return f"{self.code}({details})" if details else self.code

    def to_dict(self) -> dict:
        """Forma lista para enviar como JSON (la usará la web)."""
        return {"code": self.code, "params": self.params}


class InvalidSequenceError(DNAinventError):
    """La secuencia está vacía, tiene caracteres no permitidos o una forma incorrecta."""


class InvalidParameterError(DNAinventError):
    """Un parámetro numérico (longitud de guía, número de bits...) no es válido."""


class ServiceError(DNAinventError):
    """Un servicio externo (AlphaFold DB, BLAST...) no responde o no tiene lo pedido."""
