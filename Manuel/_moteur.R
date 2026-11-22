# ===========================================================================
# _moteur.R — copie conforme des fonctions de calcul de l'application simEFA
# (app.R). Le manuel les utilise pour que ses figures et ses nombres soient
# exactement ceux que l'étudiant obtient à l'écran. À recopier si app.R change.
# ===========================================================================

generate_efa_simulation <- function(n, num_factors, vars_per_factor,
                                    loading_strength, cross_loading,
                                    cross_proportion, error_multiplier,
                                    correlated, factor_cor, seed,
                                    loading_pattern = NULL) {
  factor_cov <- if (correlated) {
    tmp <- matrix(factor_cor, nrow = num_factors, ncol = num_factors)
    diag(tmp) <- 1
    tmp
  } else {
    diag(num_factors)
  }
  set.seed(seed)
  total_vars <- num_factors * vars_per_factor
  if (is.null(loading_pattern)) {
    loading_pattern <- matrix(0, nrow = total_vars, ncol = num_factors)
    for (f in seq_len(num_factors)) {
      variable_indices <- ((f - 1) * vars_per_factor + 1):(f * vars_per_factor)
      loading_pattern[variable_indices, f] <- loading_strength
      if (num_factors > 1 && cross_loading > 0 && cross_proportion > 0) {
        ambiguous <- variable_indices[runif(length(variable_indices)) < cross_proportion]
        for (var_idx in ambiguous) {
          # Attention : sample(x, 1) avec un x de longueur 1 tire dans 1:x.
          # On choisit donc explicitement parmi les AUTRES facteurs, sans
          # jamais écraser la charge principale de l'item.
          autres <- setdiff(seq_len(num_factors), f)
          other <- if (length(autres) == 1) autres else sample(autres, size = 1)
          loading_pattern[var_idx, other] <- cross_loading
        }
      }
    }
  }
  if (!is.matrix(loading_pattern) || nrow(loading_pattern) != total_vars ||
      ncol(loading_pattern) != num_factors) {
    stop("Le patron de charges simulées est incompatible avec les réglages demandés.")
  }
  common_variance <- diag(loading_pattern %*% factor_cov %*% t(loading_pattern))
  if (any(!is.finite(common_variance)) || any(common_variance >= 0.98)) {
    stop("Cette combinaison de charges et de corrélation est trop forte.")
  }

  if (correlated) {
    factors <- matrix(rnorm(n * num_factors), nrow = n, ncol = num_factors) %*%
      chol(factor_cov)
  } else {
    factors <- matrix(rnorm(n * num_factors), nrow = n, ncol = num_factors)
  }
  factors <- scale(factors)

  data <- matrix(0, nrow = n, ncol = total_vars)
  for (var_idx in seq_len(total_vars)) {
    signal <- factors %*% loading_pattern[var_idx, ]
    noise_sd <- sqrt(1 - common_variance[var_idx]) * error_multiplier
    data[, var_idx] <- scale(signal + rnorm(n, 0, noise_sd))
  }
  data <- as.data.frame(data)
  names(data) <- paste0("X", seq_len(total_vars))
  attr(data, "common_variance") <- common_variance
  attr(data, "loading_pattern") <- loading_pattern
  data
}

extract_standardized_efa_loadings <- function(efa_fit, item_names = NULL) {
  fit <- if (inherits(efa_fit, "efaList")) efa_fit[[1]] else efa_fit
  std <- tryCatch(lavaan::standardizedSolution(fit), error = function(e) NULL)
  if (is.null(std) || !"est.std" %in% names(std)) return(NULL)

  std <- std[std$op == "=~" & !is.na(std$est.std), , drop = FALSE]
  if (nrow(std) == 0) return(NULL)

  if (is.null(item_names)) item_names <- unique(std$rhs)
  std <- std[std$rhs %in% item_names, , drop = FALSE]
  factor_names <- unique(std$lhs)
  if (length(factor_names) == 0) return(NULL)

  loadings <- matrix(0, nrow = length(item_names), ncol = length(factor_names),
                     dimnames = list(item_names, paste0("F", seq_along(factor_names))))
  for (i in seq_len(nrow(std))) {
    row <- match(std$rhs[i], item_names)
    col <- match(std$lhs[i], factor_names)
    if (!is.na(row) && !is.na(col)) loadings[row, col] <- std$est.std[i]
  }
  residuals <- tryCatch(lavaan::standardizedSolution(fit), error = function(e) NULL)
  if (!is.null(residuals)) {
    residuals <- residuals[residuals$op == "~~" & residuals$lhs == residuals$rhs &
                             residuals$lhs %in% item_names, , drop = FALSE]
    communalities <- rep(NA_real_, length(item_names))
    names(communalities) <- item_names
    communalities[match(residuals$lhs, item_names)] <- 1 - residuals$est.std
    attr(loadings, "communalities") <- communalities
  }
  loadings
}

extract_factor_correlation <- function(efa_fit, nfactors) {
  fit <- if (inherits(efa_fit, "efaList")) efa_fit[[1]] else efa_fit
  psi <- tryCatch(lavaan::inspect(fit, "est")$psi, error = function(e) NULL)
  if (is.null(psi) || !is.matrix(psi) || nrow(psi) != nfactors) {
    return(diag(nfactors))
  }
  out <- tryCatch(stats::cov2cor(psi), error = function(e) diag(nfactors))
  diag(out) <- 1
  colnames(out) <- rownames(out) <- paste0("F", seq_len(nfactors))
  out
}

summarise_simple_structure <- function(loadings, communalities, threshold = 0.30,
                                       cross_gap = 0.15) {
  q <- ncol(loadings)
  ordre <- t(apply(abs(loadings), 1, order, decreasing = TRUE))
  if (q == 1) ordre <- matrix(1, nrow = nrow(loadings))
  idx1 <- ordre[, 1]
  prim <- loadings[cbind(seq_len(nrow(loadings)), idx1)]
  sec <- if (q > 1) loadings[cbind(seq_len(nrow(loadings)), ordre[, 2])] else rep(0, nrow(loadings))
  gap <- abs(prim) - abs(sec)
  statut <- ifelse(abs(prim) < threshold, "charge faible",
            ifelse(abs(sec) >= threshold & gap < cross_gap, "charge croisée",
            ifelse(abs(sec) >= threshold, "seconde charge notable", "structure simple")))
  data.frame(
    Item = rownames(loadings),
    Facteur_principal = paste0("F", idx1),
    Charge_principale = round(prim, 3),
    Seconde_charge = round(sec, 3),
    Écart = round(gap, 3),
    Communalité = round(communalities, 3),
    Statut = statut,
    stringsAsFactors = FALSE
  )
}

all_factor_permutations <- function(x) {
  if (length(x) <= 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_factor_permutations(x[-i]))
  }))
}

align_bootstrap_loadings <- function(reference, candidate) {
  if (!identical(dim(reference), dim(candidate))) return(candidate)
  q <- ncol(reference)
  permutations <- all_factor_permutations(seq_len(q))
  scores <- apply(permutations, 1, function(p) {
    sum(abs(colSums(reference * candidate[, p, drop = FALSE], na.rm = TRUE)))
  })
  aligned <- candidate[, permutations[which.max(scores), ], drop = FALSE]
  for (j in seq_len(q)) {
    if (sum(reference[, j] * aligned[, j], na.rm = TRUE) < 0) {
      aligned[, j] <- -aligned[, j]
    }
  }
  colnames(aligned) <- colnames(reference)
  aligned
}

ajuster_afe <- function(data_complete, nfactors, rotation, method,
                        n_total = nrow(data_complete), maxit = 1000) {
  if (ncol(data_complete) <= nfactors) {
    stop("Il faut plus d'items que de facteurs demandés.")
  }
  modele <- lavaan::efa(
    data = data_complete, nfactors = nfactors, rotation = rotation,
    estimator = method,
    # Optimiseur par défaut de lavaan (nlminb) : nettement plus rapide et plus
    # stable que BFGS lorsque l'on demande un facteur de trop.
    control = list(iter.max = maxit)
  )
  fit <- if (inherits(modele, "efaList")) modele[[1]] else modele
  L <- extract_standardized_efa_loadings(modele, names(data_complete))
  if (is.null(L)) stop("Charges standardisées indisponibles : la solution n'a pas convergé.")
  h2 <- attr(L, "communalities"); attr(L, "communalities") <- NULL
  if (is.null(h2) || all(is.na(h2))) h2 <- rowSums(L^2)
  names(h2) <- rownames(L)
  phi <- extract_factor_correlation(modele, nfactors)
  ss <- colSums(L^2)
  list(
    efa = modele, loadings = L, communality = h2, uniquenesses = 1 - h2,
    method = method, rotation = rotation, nfactors = nfactors,
    # Somme des carrés des charges (patrons) par facteur. Avec une rotation
    # oblique, ces quantités se chevauchent : leur total n'est pas la variance
    # commune. La variance commune réelle est la moyenne des communalités.
    Vaccounted = data.frame(SS_loadings = ss, Proportion = ss / nrow(L),
                            Cumulative = cumsum(ss / nrow(L))),
    variance_commune = mean(h2),
    fit_indices = tryCatch(lavaan::fitMeasures(fit), error = function(e) NULL),
    converged = isTRUE(tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)),
    factor_correlation = phi, structure_loadings = L %*% phi,
    data_complete = data_complete, variables = names(data_complete),
    n_total = n_total, n_used = nrow(data_complete)
  )
}

population_simulee <- function(loading_pattern, factor_cor, error_multiplier) {
  q <- ncol(loading_pattern)
  phi <- if (q == 1) matrix(1) else factor_cor
  commune <- diag(loading_pattern %*% phi %*% t(loading_pattern))
  var_item <- commune + (1 - commune) * error_multiplier^2
  L <- loading_pattern / sqrt(var_item)
  dimnames(L) <- list(paste0("X", seq_len(nrow(L))), paste0("F", seq_len(q)))
  list(loadings = L, communality = setNames(commune / var_item, rownames(L)), phi = phi)
}

omega_sous_echelle <- function(d) {
  m <- tryCatch(suppressWarnings(psych::fa(d, nfactors = 1, fm = "ml")), error = function(e) NULL)
  if (is.null(m)) return(list(omega = NA_real_, lambda = rep(NA_real_, ncol(d))))
  l <- as.numeric(m$loadings[, 1]); psi <- 1 - l^2
  list(omega = sum(l)^2 / (sum(l)^2 + sum(psi)), lambda = setNames(l, colnames(d)))
}

congruence_tucker <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

suggestion_parallele <- function(observees, simulees) {
  sous <- which(observees <= simulees)
  if (length(sous) == 0) length(observees) else sous[1] - 1
}

convert_efa_to_cfa <- function(efa_result, loading_threshold = 0.3,
                               crossloading_rule = c("primary", "all")) {
  crossloading_rule <- match.arg(crossloading_rule)
  # Extraire les charges factorielles 
  loadings_matrix <- efa_result$loadings
  
  # Identifier les facteurs
  factors <- colnames(loadings_matrix)
  
  # Initialiser la syntaxe CFA
  cfa_syntax <- c()
  
  primary_factor <- apply(abs(loadings_matrix), 1, which.max)
  primary_value <- apply(abs(loadings_matrix), 1, max)

  # Par défaut, la CFA candidate retient seulement l'appartenance principale.
  # Une charge croisée peut être conservée explicitement, mais elle ne doit pas
  # s'ajouter silencieusement au modèle sous prétexte qu'elle dépasse un seuil.
  for (factor_idx in 1:length(factors)) {
    factor <- factors[factor_idx]
    factor_items <- if (crossloading_rule == "all") {
      rownames(loadings_matrix)[abs(loadings_matrix[, factor_idx]) >= loading_threshold]
    } else {
      rownames(loadings_matrix)[primary_factor == factor_idx & primary_value >= loading_threshold]
    }
    
    if (length(factor_items) > 0) {
      # Construire la syntaxe du modèle CFA pour ce facteur
      cfa_syntax <- c(cfa_syntax, paste(factor, "=~", paste(factor_items, collapse=" + ")))
    }
  }
  
  # Combiner en une seule chaîne
  cfa_model <- paste(cfa_syntax, collapse = "\n")
  
  return(cfa_model)
}
