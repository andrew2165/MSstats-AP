# Synthetic run summaries: variable bait, balanced client deviations within group.
bait_example = function() {
    samples = data.frame(RUN = factor(1:6), originalRUN = paste0("sample", 1:6),
        GROUP = factor(rep(c("Control", "Treatment"), each = 3)),
        SUBJECT = paste0("subject", 1:6), LABEL = "L")
    p = samples[rep(1:6, each = 3), ]
    p$Protein = factor(rep(c("bait|id", "tracks", "changes"), 6))
    p$LogIntensities = rep(c(20, 22, 19, 23, 21, 24), each = 3) +
        rep(c(0, -2, 0), 6) + rep(rep(c(-0.2, 0, 0.2), 2), each = 3) *
        rep(c(0, 1, 1), 6) + rep(c(rep(0, 9), rep(c(0, 0, 1), 3)), 1)
    p$Variance = rep(c(0.09, 0.04, 0.04), 6)
    p$more50missing = FALSE
    p$NumMeasuredFeature = 3L
    p$TotalGroupMeasurements = 9L
    list(FeatureLevelData = data.frame(LABEL = factor("L"), INTENSITY = 100),
         ProteinLevelData = p, SummaryMethod = "linear", custom = "preserved")
}
q = bait_example()
before = q
x = normalizeToBait(q, "bait|id")
y = normalizeToBait(q, "bait|id", variance_method = "ignore")
z = normalizeToBait(q, "bait|id", variance_method = "propagate")
keep = which(q$ProteinLevelData$Protein != "bait|id")
expected = c(-2.2, -0.2, -2, 0, -1.8, 0.2, -2.2, 0.8, -2, 1, -1.8, 1.2)
expect_equal(x$ProteinLevelData$LogIntensities, expected)
expect_identical(x, y)
expect_identical(x$ProteinLevelData$LogIntensities, z$ProteinLevelData$LogIntensities)
expect_true(all(is.na(x$ProteinLevelData$Variance)))
expect_equal(z$ProteinLevelData$Variance, rep(0.13, 12))
expect_identical(q, before)
expect_identical(x$FeatureLevelData, q$FeatureLevelData)
expect_identical(x$custom, "preserved")
audit = attr(z, "BaitNormalization")
expect_identical(audit$matching_keys, c("RUN", "LABEL"))
expect_identical(audit$retained_rows, keep)
expect_identical(audit$original_variance, q$ProteinLevelData$Variance)
expect_equal(audit$bait_abundance$Variance, rep(0.09, 6))
expect_equal(audit$coverage_after$conditions$BiologicalReplicates, c(3L, 3L))
expect_true(grepl("Zero covariance", audit$covariance_assumption))
expect_equal(nrow(audit$dropped_samples), 0)

# Exercise actual public comparison and inspect the fitted model, not just Variance.
contrast = matrix(c(-1, 1), 1, dimnames = list("Treatment-Control", c("Control", "Treatment")))
for (adjusted in list(x, z)) {
    fit = groupComparison(contrast, adjusted, use_log_file = FALSE, verbose = FALSE)
    expect_equal(fit$ComparisonResult$log2FC[match(c("tracks", "changes"), fit$ComparisonResult$Protein)], c(0, 1), tolerance = 1e-10)
    expect_equal(length(fit$FittedModel), 2L)
    for (public_model in fit$FittedModel) {
        expect_true(inherits(public_model, "lm"))
        if (identical(attr(adjusted, "BaitNormalization")$variance_method, "ignore")) {
            expect_true(is.null(public_model$weights))
        } else expect_equal(public_model$weights, rep(1/0.13, 6))
    }
    prepared = MSstatsPrepareForGroupComparison(adjusted)
    single = MSstatsAP:::.prepareSingleProteinForGC(prepared[["tracks"]])
    model = MSstatsAP:::.fitModelForGroupComparison(single, FALSE, FALSE, FALSE)$full_fit
    if (identical(attr(adjusted, "BaitNormalization")$variance_method, "ignore")) {
        expect_true(is.null(model$weights))
    } else {
        expect_equal(model$weights, rep(1/0.13, 6))
    }
}

# Private copies, ordering, markers and saveRDS/readRDS.
for (as_dt in c(FALSE, TRUE)) {
    q = bait_example()
    if (as_dt) q$ProteinLevelData = data.table::as.data.table(q$ProteinLevelData)
    q$ProteinLevelData = q$ProteinLevelData[18:1, , drop = FALSE]
    before = data.table::copy(q)
    adjusted = normalizeToBait(q, "bait|id")
    expect_equal(adjusted$ProteinLevelData$LogIntensities, rev(expected))
    expect_identical(q, before)
    if (as_dt) {
        expect_true(data.table::is.data.table(adjusted$ProteinLevelData))
        data.table::set(adjusted$ProteinLevelData, i = 1L, j = "LogIntensities", value = -10)
        expect_identical(q, before)
    }
    path = tempfile(fileext = ".rds")
    saveRDS(adjusted, path)
    restored = readRDS(path)
    unlink(path)
    expect_identical(attr(restored, "BaitNormalization"), attr(adjusted, "BaitNormalization"))
    for (obj in list(restored, restored[c("FeatureLevelData", "ProteinLevelData")])) {
        expect_error(normalizeToBait(obj, "tracks"), "already been applied")
        for (plot in c("ProfilePlot", "QCPlot", "ConditionPlot")) {
            expect_error(dataProcessPlots(obj, plot), "original protein summary")
        }
    }
}
q = bait_example()
for (bait in list("absent", "bait", NA_character_, character(), c("bait|id", "tracks"), " ")) {
    expect_error(normalizeToBait(q, bait))
}
expect_error(normalizeToBait(q[-1], "bait|id"), "FeatureLevelData")
p = q$ProteinLevelData
q$ProteinLevelData = p[c(1:18, 1), ]
expect_error(normalizeToBait(q, "bait|id"), "Duplicate protein")
for (column in c("originalRUN", "GROUP", "SUBJECT", "LABEL")) {
    q = bait_example()
    q$ProteinLevelData[[column]] = as.character(q$ProteinLevelData[[column]])
    q$ProteinLevelData[[column]][1] = "conflict"
    expect_error(normalizeToBait(q, "bait|id"), "Conflicting|labels")
}
q = bait_example()
q$ProteinLevelData$originalRUN[4:6] = "sample1"
expect_error(normalizeToBait(q, "bait|id"), "Conflicting")
for (column in c("RUN", "originalRUN", "GROUP", "SUBJECT", "LABEL", "Protein")) {
    q = bait_example()
    q$ProteinLevelData[[column]][1] = NA
    expect_error(normalizeToBait(q, "bait|id"), "Missing or blank")
}
q = bait_example()
attr(q$ProteinLevelData, "MSstatsAbundanceScale") = "log2"
expect_equal(normalizeToBait(q, "bait|id")$ProteinLevelData$LogIntensities, expected)
attr(q$ProteinLevelData, "MSstatsAbundanceScale") = "log10"
expect_error(normalizeToBait(q, "bait|id"), "log2")
# Only the protein-table scale attribute is recorded, including after rewrapping.
for (log_base in c(2, 10)) {
    q = dataProcess(DDARawData, logTrans = log_base, use_log_file = FALSE,
                    verbose = FALSE)
    expect_true(is.null(attr(q, "MSstatsLogBase", exact = TRUE)))
    repacked = unserialize(serialize(q[c("FeatureLevelData", "ProteinLevelData")], NULL))
    expect_identical(attr(repacked$ProteinLevelData, "MSstatsAbundanceScale", exact = TRUE),
                     paste0("log", log_base))
    if (log_base == 10) {
        expect_error(normalizeToBait(repacked, as.character(q$ProteinLevelData$Protein[1])), "log2")
    }
}
q = bait_example()
q$FeatureLevelData$LABEL = factor("H")
expect_error(normalizeToBait(q, "bait|id"), "labels")

# Missing-bait policy precedes variance validation; drops whole runs only.
for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    q = bait_example()
    q$ProteinLevelData$LogIntensities[1] = bad
    q$ProteinLevelData$Variance[1:3] = -1
    expect_error(normalizeToBait(q, "bait|id"), "sample1")
    expect_warning(d <- normalizeToBait(q, "bait|id", "drop", "propagate"), "Dropped 2 client row")
    a = attr(d, "BaitNormalization")
    expect_equal(a$dropped_samples$ClientRows, 2L)
    expect_equal(a$coverage_after$conditions$BiologicalReplicates, c(2L, 3L))
}
q = bait_example()
q$ProteinLevelData = q$ProteinLevelData[-1, ]
expect_error(normalizeToBait(q, "bait|id"), "absent")
expect_warning(d <- normalizeToBait(q, "bait|id", "drop"), "Dropped 2")
expect_identical(attr(d, "BaitNormalization")$dropped_samples$Reason, "absent")
q = bait_example()
q$ProteinLevelData$LogIntensities[c(10, 13, 16)] = NA
expect_warning(d <- normalizeToBait(q, "bait|id", "drop"), "Lost entire condition.*Treatment")
expect_identical(attr(d, "BaitNormalization")$lost_conditions, "Treatment")
expect_error(groupComparison(contrast, d, use_log_file = FALSE, verbose = FALSE))
q$ProteinLevelData$LogIntensities[c(1, 4, 7)] = NA
expect_error(normalizeToBait(q, "bait|id", "drop"), "No usable client")
q = bait_example()
q$ProteinLevelData$LogIntensities[keep] = NA
expect_error(normalizeToBait(q, "bait|id"), "No usable client")
q$ProteinLevelData$LogIntensities[2] = Inf
expect_error(normalizeToBait(q, "bait|id"), "Infinite adjusted")

# Variance validation, missing clients and numerical weight floor.
for (bad in c(NA_real_, NaN, -1, Inf, -Inf)) {
    for (row in c(1, 2)) {
        q = bait_example()
        q$ProteinLevelData$Variance[row] = bad
        expect_true(all(is.na(normalizeToBait(q, "bait|id")$ProteinLevelData$Variance)))
        expect_error(normalizeToBait(q, "bait|id", variance_method = "propagate"), "ignore")
    }
}
q = bait_example()
q$ProteinLevelData$Variance = NULL
expect_true(all(is.na(normalizeToBait(q, "bait|id")$ProteinLevelData$Variance)))
expect_error(normalizeToBait(q, "bait|id", variance_method = "propagate"), "ignore")
q = bait_example()
q$ProteinLevelData$LogIntensities[2] = NA
q$ProteinLevelData$Variance[2] = NA
expect_true(is.na(normalizeToBait(q, "bait|id", variance_method = "propagate")$ProteinLevelData$LogIntensities[1]))
q = bait_example()
for (variance in c(0, 1e-8)) {
    q$ProteinLevelData$Variance = variance
    d = normalizeToBait(q, "bait|id", variance_method = "propagate")
    single = MSstatsAP:::.prepareSingleProteinForGC(MSstatsPrepareForGroupComparison(d)[["tracks"]])
    fit = MSstatsAP:::.fitModelForGroupComparison(single, FALSE, FALSE, FALSE)$full_fit
    expect_equal(fit$weights, rep(1000, 6))
}

# Repeated subjects and technical runs remain separate and reach mixed models.
q = bait_example()
q$ProteinLevelData$SUBJECT = rep(rep(c("s1", "s2", "s3"), 2), each = 3)
q$ProteinLevelData$LogIntensities[10:18] = q$ProteinLevelData$LogIntensities[10:18] +
    rep(c(-0.1, 0.04, 0.06), each = 3) * rep(c(0, 1, 1), 3)
d = normalizeToBait(q, "bait|id")
expect_true(checkRepeatedDesign(d))
expect_identical(d$ProteinLevelData$SUBJECT, q$ProteinLevelData$SUBJECT[keep])
expect_equal(nrow(groupComparison(contrast, d, use_log_file = FALSE, verbose = FALSE)$ComparisonResult), 2L)
p2 = q$ProteinLevelData
p2$RUN = paste0(p2$RUN, "tech2")
p2$originalRUN = paste0(p2$originalRUN, "tech2")
p2$LogIntensities = p2$LogIntensities + rep(c(1, 1.1, 0.9), 6) +
    rep(c(-0.04, 0.02, 0.06, 0.01, -0.03, 0.04), each = 3) * rep(c(0, 1, 1), 6)
q$ProteinLevelData$RUN = as.character(q$ProteinLevelData$RUN)
q$ProteinLevelData = rbind(q$ProteinLevelData, p2)
d = normalizeToBait(q, "bait|id", variance_method = "propagate")
expect_equal(nrow(d$ProteinLevelData), 24L)
expect_equal(attr(d, "BaitNormalization")$coverage_after$conditions$Samples, c(6L, 6L))
expect_equal(nrow(groupComparison(contrast, d, use_log_file = FALSE, verbose = FALSE)$ComparisonResult), 2L)

mixed = groupComparison(contrast, d, use_log_file = FALSE, verbose = FALSE)
expect_true(all(is.finite(mixed$ComparisonResult$log2FC)))
for (model in mixed$FittedModel) {
    expect_true(inherits(model, "merMod"))
    expect_equal(unname(stats::weights(model)), rep(1/0.13, 12))
}
# Exact sample identifiers are independent of data.table's global rounding.
q = bait_example()
q$ProteinLevelData$RUN = rep(c(1, 1 + 1e-12, 2, 3, 4, 5), each = 3)
previous_rounding = data.table::getNumericRounding()
data.table::setNumericRounding(2L)
tryCatch({
    d = normalizeToBait(q, "bait|id")
    expect_equal(d$ProteinLevelData$LogIntensities, expected)
}, finally = data.table::setNumericRounding(previous_rounding))
q = bait_example()
q$ProteinLevelData$Variance = "invalid"
expect_true(all(is.na(normalizeToBait(q, "bait|id")$ProteinLevelData$Variance)))
expect_error(normalizeToBait(q, "bait|id", variance_method = "propagate"), "ignore")
q = bait_example()
q$ProteinLevelData$LogIntensities = as.character(q$ProteinLevelData$LogIntensities)
expect_error(normalizeToBait(q, "bait|id"), "numeric vector")

# Explicit contrasts must be rebuilt after an entire third condition is dropped.
q = bait_example()
p = q$ProteinLevelData
p$RUN = paste0(p$RUN, "lost")
p$originalRUN = paste0(p$originalRUN, "lost")
p$SUBJECT = paste0(p$SUBJECT, "lost")
p$GROUP = "Lost"
p$LogIntensities[p$Protein == "bait|id"] = NA_real_
q$ProteinLevelData$RUN = as.character(q$ProteinLevelData$RUN)
q$ProteinLevelData$GROUP = as.character(q$ProteinLevelData$GROUP)
q$ProteinLevelData = rbind(q$ProteinLevelData, p)
expect_warning(d <- normalizeToBait(q, "bait|id", "drop"), "Lost entire condition.*Lost")
expect_equal(sum(attr(d, "BaitNormalization")$dropped_samples$ClientRows), 12L)
expect_identical(attr(d, "BaitNormalization")$lost_conditions, "Lost")
expect_error(groupComparison(cbind(contrast, Lost = 0), d, use_log_file = FALSE, verbose = FALSE))
r = groupComparison(contrast, d, use_log_file = FALSE, verbose = FALSE)
expect_equal(r$ComparisonResult$log2FC[match(c("tracks", "changes"), r$ComparisonResult$Protein)], c(0, 1), tolerance = 1e-10)
