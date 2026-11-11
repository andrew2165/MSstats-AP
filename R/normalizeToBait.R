#' Normalize label-free protein summaries to an affinity-purification bait
#'
#' Subtract the bait's summarized log2 abundance from each client in the same
#' summarized run, after fraction merging and ordinary normalization.
#'
#' @param data Output of \code{\link{dataProcess}} with \code{logTrans = 2}.
#'   Recorded scale is checked using the \code{MSstatsAbundanceScale}
#'   attribute on \code{ProteinLevelData}.
#'   Unmarked older summaries must also contain log2 abundances; the base cannot
#'   be inferred from their values.
#' @param bait One exact identifier from \code{ProteinLevelData$Protein}.
#' @param missing_bait Either \code{"error"} (default) or \code{"drop"}.
#'   Absent or non-finite bait in a run containing clients causes an error, or
#'   an explicit warning and removal of all client rows in that run. An entirely
#'   absent bait always errors. Check surviving conditions and subjects.
#' @param variance_method Either \code{"ignore"} (default), which ignores
#'   summary-level inverse-variance weights, or \code{"propagate"}, which uses
#'   client variance plus bait variance as an approximate adjusted variance.
#'
#' @details
#' Only single-label, label-free summaries are supported. Matching keys are
#' \code{RUN} and \code{LABEL}. Each RUN must map to exactly one originalRUN,
#' GROUP, SUBJECT and LABEL, and originalRUN must map to exactly one RUN.
#' Technical runs remain distinct; SUBJECT is preserved for repeated measures.
#' Duplicate protein observations and conflicting annotations are rejected.
#' Missing client abundances and zero or negative log2 ratios are retained.
#'
#' With \code{variance_method = "ignore"}, adjusted Variance is entirely NA,
#' disabling summary-level weights in \code{\link{groupComparison}}. The model
#' still estimates variation among normalized biological replicates, but does
#' not explicitly propagate individual client and bait summary uncertainty.
#'
#' With \code{variance_method = "propagate"}, the approximation is
#' \code{Var(client - bait) = Var(client) + Var(bait)}, assuming zero covariance
#' between summary errors. In general the formula also subtracts
#' \code{2 * Cov(client, bait)}. This approximation is neither exact nor always
#' conservative. Finite nonnegative client and bait variances are required for
#' every retained finite client abundance. Missing clients need no variance.
#' Missing-bait policy is applied first; invalid variances never drop samples
#' or trigger fallback. Existing model weights remain
#' \code{1 / pmax(Variance, 0.001)}: zero and very small adjusted variances
#' receive weight 1000. Both policies give identical abundances.
#'
#' Use the original summary for \code{\link{dataProcessPlots}} and bait QC,
#' and the adjusted summary for groupComparison. FeatureLevelData stays on its
#' original scale. Finite imputed bait summaries are used as supplied. When
#' present, NumImputedFeature and feature-level predicted/censored information
#' can help inspect imputation; older summaries may lack this information.
#' The function cannot reconstruct unavailable imputation provenance.
#'
#' Comparisons describe client recovery relative to recovered bait, not
#' absolute molecular stoichiometry. The method assumes prey recovery scales
#' with bait recovery. Saturation, inconsistent bait peptide coverage, and
#' noisy bait measurements can undermine this assumption. Inspect original
#' bait abundance and peptide coverage before interpreting comparisons.
#' This personal feature operates after summarization; published workflows may
#' normalize earlier. The variance policies are explicit software choices.
#'
#' saveRDS/readRDS preserve normalization and scale attributes. Rewrapping the
#' tables keeps the protein marker but loses the outer audit; rebuilding tables
#' may lose both markers and defeat repeat-normalization and plotting checks.
#'
#' @return The usual summarization list, with bait rows removed and adjusted
#'   ProteinLevelData in input row order. Unused Protein and GROUP factor levels
#'   are removed for downstream comparisons. The \code{BaitNormalization}
#'   attribute contains bait, matching_keys, missing_bait, variance_method,
#'   original_abundance and original_variance in input row order, retained_rows
#'   (input indices), bait_abundance (original bait rows including available QC
#'   and variances), dropped_samples (annotations, reason, bait abundance and
#'   ClientRows), coverage_before and coverage_after (samples, condition counts
#'   and condition/subject counts), lost_conditions, adjusted_variance, and
#'   covariance_assumption. Original metadata and FeatureLevelData are preserved.
#' @seealso \code{\link{dataProcess}}, \code{\link{groupComparison}}
#' @export
#' @examples
#' \dontrun{
#' summarized <- dataProcess(converted, logTrans = 2)
#' adjusted <- normalizeToBait(summarized, bait = "exact protein identifier")
#' audit <- attr(adjusted, "BaitNormalization")
#' result <- groupComparison(contrast.matrix = contrasts, data = adjusted)
#' }
normalizeToBait = function(data, bait, missing_bait = c("error", "drop"),
                           variance_method = c("ignore", "propagate")) {
    missing_bait = match.arg(missing_bait)
    variance_method = match.arg(variance_method)
    if (!is.list(data) || !is.data.frame(data[["FeatureLevelData"]]) ||
        !is.data.frame(data[["ProteinLevelData"]])) {
        stop("Supply a dataProcess() result with FeatureLevelData and ProteinLevelData.",
             call. = FALSE)
    }
    if (.hasBaitNormalization(data)) {
        stop("Bait normalization has already been applied. Use the original protein summary.",
             call. = FALSE)
    }
    if (!is.character(bait) || length(bait) != 1L || is.na(bait) ||
        !nzchar(trimws(bait))) {
        stop("Supply one nonempty bait protein identifier.", call. = FALSE)
    }
    scale = attr(data$ProteinLevelData, "MSstatsAbundanceScale", exact = TRUE)
    if (!is.null(scale) && !identical(scale, "log2")) {
        stop("Bait normalization requires log2 input. Rerun dataProcess(logTrans = 2).",
             call. = FALSE)
    }
    # A private copy is essential: callers may supply data.tables.
    proteins = data.table::copy(data$ProteinLevelData)
    input = as.data.frame(proteins)
    keys = c("RUN", "LABEL")
    annotations = c("RUN", "originalRUN", "GROUP", "SUBJECT", "LABEL")
    required = c(annotations, "Protein", "LogIntensities")
    absent = setdiff(required, names(input))
    if (length(absent)) {
        stop("Missing ProteinLevelData columns: ", paste(absent, collapse = ", "),
             call. = FALSE)
    }
    for (column in c(annotations, "Protein")) {
        value = as.character(input[[column]])
        if (anyNA(value) || any(!nzchar(trimws(value)))) {
            stop("Missing or blank values in ProteinLevelData$", column, call. = FALSE)
        }
    }
    if (!is.numeric(input$LogIntensities) || !is.null(dim(input$LogIntensities))) {
        stop("LogIntensities must be a numeric vector of log2 abundances.", call. = FALSE)
    }
    labels = unique(as.character(input$LABEL))
    feature_labels = unique(as.character(data$FeatureLevelData$LABEL))
    if (length(labels) != 1L || length(feature_labels) > 1L ||
        (length(feature_labels) && !identical(labels, feature_labels))) {
        stop("Only single-label label-free summaries are supported; multiple or inconsistent labels are unsupported.",
             call. = FALSE)
    }
    # Exact integer codes avoid delimiter collisions and data.table numeric rounding.
    codes = as.data.frame(lapply(input[c(annotations, "Protein")],
                                function(x) match(x, unique(x))))
    sample_rows = which(!duplicated(codes$RUN))
    sample_id = match(codes$RUN, codes$RUN[sample_rows])
    if (any(vapply(annotations, function(x)
        any(codes[[x]] != codes[[x]][sample_rows][sample_id]), logical(1))) ||
        anyDuplicated(codes$originalRUN[sample_rows])) {
        stop("Conflicting RUN, originalRUN, GROUP, SUBJECT or LABEL sample annotations.",
             call. = FALSE)
    }
    duplicate = anyDuplicated(codes[c(keys, "Protein")])
    if (duplicate) {
        stop("Duplicate protein observation at input row ", duplicate,
             " (Protein / RUN / LABEL).", call. = FALSE)
    }
    is_bait = as.character(input$Protein) == bait
    if (!any(is_bait)) {
        stop("Bait protein ", sQuote(bait), " is absent from ProteinLevelData. Check the exact Protein identifier.",
             call. = FALSE)
    }
    bait_rows = which(is_bait)
    bait_match = bait_rows[match(seq_along(sample_rows), sample_id[bait_rows])]
    sample_bait = input$LogIntensities[bait_match]
    clients = which(!is_bait)
    counts = tabulate(sample_id[clients], length(sample_rows))
    dropped = which(!is.finite(sample_bait) & counts > 0L)
    dropped_samples = input[sample_rows[dropped], annotations, drop = FALSE]
    dropped_samples$BaitAbundance = sample_bait[dropped]
    dropped_samples$Reason = ifelse(is.na(bait_match[dropped]), "absent", "non-finite")
    dropped_samples$ClientRows = counts[dropped]
    if (length(dropped)) {
        details = apply(utils::head(dropped_samples, 10L), 1L, paste, collapse = " / ")
        problem = paste0("Bait ", sQuote(bait), " is absent or non-finite in ",
                         length(dropped), " sample(s) (RUN / originalRUN / GROUP / SUBJECT / LABEL / bait / reason / client rows): ",
                         paste(details, collapse = "; "), ".")
        if (missing_bait == "error") stop(problem, call. = FALSE)
    }
    retained = clients[!sample_id[clients] %in% dropped]
    adjusted = as.double(input$LogIntensities[retained]) - sample_bait[sample_id[retained]]
    if (any(is.infinite(adjusted))) {
        stop("Infinite adjusted client abundance. Check client and bait summaries.", call. = FALSE)
    }
    if (!any(is.finite(adjusted))) {
        stop("No usable client measurements remain after bait normalization.", call. = FALSE)
    }
    original_variance = input$Variance
    adjusted_variance = rep(NA_real_, length(retained))
    if (variance_method == "propagate") {
        if (!is.numeric(original_variance) || !is.null(dim(original_variance))) {
            stop("Propagation requires numeric client and bait Variance; use variance_method = \"ignore\".", call. = FALSE)
        }
        client_var = original_variance[retained]
        bait_var = original_variance[bait_match[sample_id[retained]]]
        finite_client = is.finite(adjusted)
        invalid = finite_client & (!is.finite(client_var) | client_var < 0 |
                                  !is.finite(bait_var) | bait_var < 0 |
                                  !is.finite(client_var + bait_var))
        if (any(invalid)) {
            details = apply(utils::head(input[retained[invalid], c("Protein", annotations), drop = FALSE], 10L),
                            1L, paste, collapse = " / ")
            stop("Invalid or unavailable client/bait Variance for retained observations (Protein / RUN / originalRUN / GROUP / SUBJECT / LABEL): ",
                 paste(details, collapse = "; "),
                 ". Require finite nonnegative variances; use variance_method = \"ignore\".", call. = FALSE)
        }
        adjusted_variance[finite_client] = client_var[finite_client] + bait_var[finite_client]
    }
    before = .baitCoverage(input, clients, sample_id, sample_rows)
    after = .baitCoverage(input, retained, sample_id, sample_rows)
    lost = setdiff(as.character(before$conditions$GROUP), as.character(after$conditions$GROUP))
    if (length(dropped)) {
        warning(problem, " Dropped ", sum(counts[dropped]), " client row(s). ",
                if (length(lost)) paste0("Lost entire condition(s): ", paste(lost, collapse = ", "), ". ") else "",
                "Check condition and biological-replicate coverage in the BaitNormalization audit before constructing contrasts.",
                call. = FALSE)
    }
    output = proteins[retained, , drop = FALSE]
    output$LogIntensities = adjusted
    output$Variance = adjusted_variance
    for (column in c("Protein", "GROUP")) {
        if (is.factor(output[[column]])) output[[column]] = droplevels(output[[column]])
    }
    attr(output, "MSstatsAbundanceScale") = "log2_client_bait"
    data$ProteinLevelData = output
    attr(data, "BaitNormalization") = list(
        bait = bait, matching_keys = keys, missing_bait = missing_bait,
        variance_method = variance_method, original_abundance = input$LogIntensities,
        original_variance = original_variance, retained_rows = retained,
        bait_abundance = input[is_bait, , drop = FALSE],
        dropped_samples = dropped_samples, coverage_before = before,
        coverage_after = after, lost_conditions = lost,
        adjusted_variance = if (variance_method == "propagate") adjusted_variance else NULL,
        covariance_assumption = if (variance_method == "propagate")
            "Zero covariance between client and bait summary errors; Var(client - bait) approximated by Var(client) + Var(bait)." else
            "Summary-level variance weights ignored; client/bait summary uncertainty is not explicitly propagated."
    )
    data
}

.hasBaitNormalization = function(data) {
    !is.null(attr(data, "BaitNormalization", exact = TRUE)) ||
        (is.list(data) && identical(attr(data[["ProteinLevelData"]],
                                        "MSstatsAbundanceScale", exact = TRUE),
                                   "log2_client_bait"))
}

.baitCoverage = function(input, rows, sample_id, sample_rows) {
    ids = unique(sample_id[rows])
    samples = input[sample_rows[ids], c("RUN", "originalRUN", "GROUP", "SUBJECT", "LABEL"), drop = FALSE]
    samples$ClientRows = tabulate(match(sample_id[rows], ids), length(ids))
    samples$FiniteClientRows = tabulate(match(sample_id[rows[is.finite(input$LogIntensities[rows])]], ids), length(ids))
    groups = unique(as.character(samples$GROUP))
    conditions = data.frame(GROUP = groups,
        Samples = vapply(groups, function(g) sum(samples$GROUP == g), integer(1)),
        BiologicalReplicates = vapply(groups, function(g) length(unique(samples$SUBJECT[samples$GROUP == g])), integer(1)),
        UsableSamples = vapply(groups, function(g) sum(samples$GROUP == g & samples$FiniteClientRows > 0L), integer(1)),
        row.names = NULL)
    subjects = unique(samples[c("GROUP", "SUBJECT")])
    subjects$Samples = vapply(seq_len(nrow(subjects)), function(i)
        sum(samples$GROUP == subjects$GROUP[i] & samples$SUBJECT == subjects$SUBJECT[i]), integer(1))
    list(samples = samples, conditions = conditions, biological_replicates = subjects)
}
