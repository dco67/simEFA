library(shiny)
library(bslib)
library(EFAfactors)
library(DT)
library(ggplot2)
library(reshape2)
library(corrplot)
library(EFA.dimensions)
library(utils)
library(plotly)
library(nortest)
library(base64enc)

# Fonction pour exécuter les tests de dimensionnalité avec RAWPAR
run_dimtests <- function(data, use_rawpar = TRUE, randtype = "generated", 
                         extraction = "PAF", ndatasets = 100, 
                         percentile = 95, corkind = "pearson", 
                         corkindRAND = "pearson") {
  
  # Calculer les valeurs propres de la matrice de corrélation
  cor_matrix <- cor(data, method = corkind)
  eigenvalues <- eigen(cor_matrix)$values
  
  # Résultats à retourner
  results <- list(
    eigenvalues = eigenvalues
  )
  
  # Test parallèle avec RAWPAR si disponible
  if (use_rawpar && requireNamespace("EFA.dimensions", quietly = TRUE)) {
    tryCatch({
      # Exécuter RAWPAR
      rawpar_results <- EFA.dimensions::RAWPAR(
        data = data,
        randtype = randtype,
        extraction = extraction,
        Ndatasets = ndatasets,
        percentile = percentile,
        corkind = corkind,
        corkindRAND = corkindRAND
      )
      
      # Vérifier que rawpar_results est une liste avant d'accéder à ses éléments
      if (is.list(rawpar_results) && !is.null(rawpar_results$eigenvalues)) {
        # Stocker les résultats du test parallèle
        results$parallel_data <- list(
          Mean = rawpar_results$eigenvalues$Random.Mean,
          Percentile = rawpar_results$eigenvalues$Random.Percentile,
          NfactorsPA = rawpar_results$NfactorsPA
        )
      } else {
        warning("Format de résultat RAWPAR inattendu")
      }
    }, error = function(e) {
      warning("Erreur lors de l'exécution de RAWPAR : ", e$message)
    })
  } else {
    # Méthode alternative si RAWPAR n'est pas disponible
    # (code existant pour l'analyse parallèle)
  }
  
  return(results)
}

# Modification de la fonction compare_factor_models
compare_factor_models <- function(data, min_factors = 1, max_factors = 5, 
                                  extraction = "paf", rotation = "promax", 
                                  corkind = "pearson") {
  
  results <- list()
  
  # Stocker les résultats pour chaque nombre de facteurs
  for (nfactors in min_factors:max_factors) {
    tryCatch({
      model <- EFA.dimensions::EFA(
        data = data,
        extraction = extraction,
        corkind = corkind,
        Nfactors = nfactors,
        rotation = rotation,
        verbose = FALSE
      )
      
      # Si model$loadingsNOROT n'existe pas mais qu'on a model$loadings,
      # créer une copie pour assurer la compatibilité
      if (!exists("loadingsNOROT", where = model) && exists("loadings", where = model)) {
        model$loadingsNOROT <- model$loadings
      }
      
      # Extraire les indices d'ajustement
      fit_indices <- c()
      if (!is.null(model$fit_coefs)) {
        fit_coefs <- model$fit_coefs
        # Indices d'ajustement basiques
        if (!is.null(fit_coefs$RMSR)) fit_indices["RMSR"] <- fit_coefs$RMSR
        if (!is.null(fit_coefs$RMSEA)) fit_indices["RMSEA"] <- fit_coefs$RMSEA
        if (!is.null(fit_coefs$TLI)) fit_indices["TLI"] <- fit_coefs$TLI
        if (!is.null(fit_coefs$CFI)) fit_indices["CFI"] <- fit_coefs$CFI
        if (!is.null(fit_coefs$GFI)) fit_indices["GFI"] <- fit_coefs$GFI
        if (!is.null(fit_coefs$CAF)) fit_indices["CAF"] <- fit_coefs$CAF
        
        # Critères d'information
        if (!is.null(fit_coefs$AIC)) fit_indices["AIC"] <- fit_coefs$AIC
        if (!is.null(fit_coefs$BIC)) fit_indices["BIC"] <- fit_coefs$BIC
        if (!is.null(fit_coefs$CAIC)) fit_indices["CAIC"] <- fit_coefs$CAIC
        if (!is.null(fit_coefs$SABIC)) fit_indices["SABIC"] <- fit_coefs$SABIC
      }
      
      # Calculer la variance expliquée
      variance_explained <- NA
      if (!is.null(model$varexplROT)) {
        if (is.matrix(model$varexplROT)) {
          # Si c'est une matrice (format complet avec colonnes)
          variance_explained <- model$varexplROT[nrow(model$varexplROT), 3] * 100
        } else if (is.vector(model$varexplROT)) {
          # Si c'est un vecteur (sommes des carrés par facteur)
          variance_explained <- sum(model$varexplROT) * 100
        }
      } else if (!is.null(model$varexplNOROT2) && is.matrix(model$varexplNOROT2)) {
        variance_explained <- model$varexplNOROT2[nrow(model$varexplNOROT2), 3] * 100
      } else if (!is.null(model$varexplNOROT1) && is.matrix(model$varexplNOROT1)) {
        variance_explained <- sum(model$varexplNOROT1[1:model$nfactors, 2]) * 100
      } else {
        variance_explained <- NA
      }
      
      # Stocker les résultats
      results[[paste0("factors_", nfactors)]] <- list(
        nfactors = nfactors,
        fit_indices = fit_indices,
        variance_explained = variance_explained,
        model = model  # Stocker le modèle complet pour référence
      )
      
    }, error = function(e) {
      message(paste("Erreur pour", nfactors, "facteurs:", e$message))
      results[[paste0("factors_", nfactors)]] <- list(
        nfactors = nfactors,
        error = e$message
      )
    })
  }
  
  # Créer un dataframe des résultats pour faciliter la comparaison
  comparison_df <- data.frame(
    Factors = min_factors:max_factors,
    stringsAsFactors = FALSE
  )
  
  # Remplir le dataframe avec les indices d'ajustement
  all_metrics <- unique(unlist(lapply(results, function(x) {
    if (!is.null(x$fit_indices)) names(x$fit_indices) else NULL
  })))
  
  for (metric in all_metrics) {
    comparison_df[[metric]] <- NA
    for (i in 1:nrow(comparison_df)) {
      nfactors <- comparison_df$Factors[i]
      result_key <- paste0("factors_", nfactors)
      if (!is.null(results[[result_key]]) && 
          !is.null(results[[result_key]]$fit_indices) && 
          !is.null(results[[result_key]]$fit_indices[metric])) {
        comparison_df[i, metric] <- results[[result_key]]$fit_indices[metric]
      }
    }
  }
  
  # Ajouter la variance expliquée
  comparison_df$Variance_Explained <- NA
  for (i in 1:nrow(comparison_df)) {
    nfactors <- comparison_df$Factors[i]
    result_key <- paste0("factors_", nfactors)
    if (!is.null(results[[result_key]]) && 
        !is.null(results[[result_key]]$variance_explained)) {
      comparison_df[i, "Variance_Explained"] <- results[[result_key]]$variance_explained
    }
  }
  
  return(list(
    results = results,
    comparison_df = comparison_df
  ))
}

ui <- page_sidebar(
  title = "Simulation de données pour l'analyse factorielle exploratoire",
  sidebar = sidebar(
    width = 500,
    
    # Préréglages avec design amélioré en haut
    card(
      card_header(
        div(
          span(icon("sliders"), style = "margin-right: 5px;"),
          "Préréglages de simulation"
        )
      ),
      card_body(
        p("Appliquez rapidement des configurations prédéfinies :"),
        
        div(
          style = "display: grid; grid-template-columns: 1fr 1fr; gap: 10px;",
          
          # Préréglage Personnalité
          div(
            class = "preset-card",
            style = "border: 1px solid #e0e0e0; border-radius: 8px; padding: 12px; cursor: pointer; background-color: #f8f9fa; transition: all 0.2s;",
            onclick = "$(this).css('background-color', '#e2f0fd').delay(200).queue(function(next){ $(this).css('background-color', '#f8f9fa'); next(); }); $('#preset_personality').click();",
            
            div(
              style = "display: flex; align-items: center;",
              div(
                style = "background-color: #17a2b8; color: white; width: 40px; height: 40px; border-radius: 50%; display: flex; align-items: center; justify-content: center; margin-right: 10px;",
                icon("user")
              ),
              div(
                style = "flex-grow: 1;",
                div(style = "font-weight: bold; color: #333;", "Personnalité"),
                div(style = "font-size: 0.8rem; color: #666;", "5 facteurs - BFI")
              )
            ),
            actionButton("preset_personality", "", style = "display: none;")
          ),
          
          # Préréglage Cognitif
          div(
            class = "preset-card",
            style = "border: 1px solid #e0e0e0; border-radius: 8px; padding: 12px; cursor: pointer; background-color: #f8f9fa; transition: all 0.2s;",
            onclick = "$(this).css('background-color', '#e2f0fd').delay(200).queue(function(next){ $(this).css('background-color', '#f8f9fa'); next(); }); $('#preset_cognitive').click();",
            
            div(
              style = "display: flex; align-items: center;",
              div(
                style = "background-color: #28a745; color: white; width: 40px; height: 40px; border-radius: 50%; display: flex; align-items: center; justify-content: center; margin-right: 10px;",
                icon("brain")
              ),
              div(
                style = "flex-grow: 1;",
                div(style = "font-weight: bold; color: #333;", "Tests cognitifs"),
                div(style = "font-size: 0.8rem; color: #666;", "3 facteurs - Aptitudes")
              )
            ),
            actionButton("preset_cognitive", "", style = "display: none;")
          ),
          
          # Préréglage Attitudes
          div(
            class = "preset-card",
            style = "border: 1px solid #e0e0e0; border-radius: 8px; padding: 12px; cursor: pointer; background-color: #f8f9fa; transition: all 0.2s;",
            onclick = "$(this).css('background-color', '#e2f0fd').delay(200).queue(function(next){ $(this).css('background-color', '#f8f9fa'); next(); }); $('#preset_social').click();",
            
            div(
              style = "display: flex; align-items: center;",
              div(
                style = "background-color: #fd7e14; color: white; width: 40px; height: 40px; border-radius: 50%; display: flex; align-items: center; justify-content: center; margin-right: 10px;",
                icon("comments")
              ),
              div(
                style = "flex-grow: 1;",
                div(style = "font-weight: bold; color: #333;", "Échelles d'attitudes"),
                div(style = "font-size: 0.8rem; color: #666;", "4 facteurs - Opinion")
              )
            ),
            actionButton("preset_social", "", style = "display: none;")
          ),
          
          # Préréglage Performance
          div(
            class = "preset-card",
            style = "border: 1px solid #e0e0e0; border-radius: 8px; padding: 12px; cursor: pointer; background-color: #f8f9fa; transition: all 0.2s;",
            onclick = "$(this).css('background-color', '#e2f0fd').delay(200).queue(function(next){ $(this).css('background-color', '#f8f9fa'); next(); }); $('#preset_academic').click();",
            
            div(
              style = "display: flex; align-items: center;",
              div(
                style = "background-color: #6f42c1; color: white; width: 40px; height: 40px; border-radius: 50%; display: flex; align-items: center; justify-content: center; margin-right: 10px;",
                icon("graduation-cap")
              ),
              div(
                style = "flex-grow: 1;",
                div(style = "font-weight: bold; color: #333;", "Performance académique"),
                div(style = "font-size: 0.8rem; color: #666;", "2 facteurs - Education")
              )
            ),
            actionButton("preset_academic", "", style = "display: none;")
          )
        ),
        
        tags$script(HTML("
          $(document).ready(function() {
            $('.preset-card').hover(
              function() { $(this).css('transform', 'translateY(-3px)').css('box-shadow', '0 4px 8px rgba(0,0,0,0.1)'); },
              function() { $(this).css('transform', 'translateY(0)').css('box-shadow', 'none'); }
            );
          });
        ")),
        
        tags$hr(style = "margin-top: 0.5rem; margin-bottom: 0.5rem;"),
        actionButton("generate", "Générer les données", 
                     icon = icon("play"),
                     class = "btn-primary", 
                     width = "100%")
      )
    ),
    
    # Paramètres de simulation dans un accordéon
    accordion(
      open = FALSE,
      
      # Structure des facteurs
      accordion_panel(
        "Structure des facteurs",
        layout_columns(
          col_widths = c(6, 6),
          numericInput("nfact", "Nombre de facteurs:", 
                       value = 3, min = 1, max = 10, step = 1),
          selectInput("vpf_type", "Variables par facteur:",
                      choices = c("Fixe" = "fixed", 
                                  "Aléa petit (5-10)" = "S",
                                  "Aléa moyen (5-15)" = "M", 
                                  "Aléa grand (5-20)" = "L",
                                  "Aléa personnalisé" = "R"),
                      selected = "M")
        ),
        
        conditionalPanel(
          condition = "input.vpf_type == 'fixed'",
          numericInput("vpf_fixed", "Nombre de variables par facteur:", 
                       value = 5, min = 3, max = 20, step = 1)
        ),
        
        conditionalPanel(
          condition = "input.vpf_type == 'R'",
          layout_columns(
            col_widths = c(6, 6),
            numericInput("low_vpf", "Minimum variables:", 
                         value = 5, min = 3, max = 19, step = 1),
            numericInput("up_vpf", "Maximum variables:", 
                         value = 15, min = 4, max = 20, step = 1)
          )
        )
      ),
      
      # Distribution et échantillon
      accordion_panel(
        "Distribution et échantillon",
        numericInput("N", "Nombre d'observations:", 
                     value = 500, min = 100, max = 5000, step = 100),
        
        selectInput("distri", "Distribution des données:",
                    choices = c("Normale" = "normal", "Beta" = "beta"),
                    selected = "normal"),
        
        conditionalPanel(
          condition = "input.distri == 'beta'",
          selectInput("beta_params", "Paramètres Beta:",
                      choices = c("Aléatoires" = "random", "Personnalisés" = "custom"),
                      selected = "random"),
          conditionalPanel(
            condition = "input.beta_params == 'custom'",
            layout_columns(
              col_widths = c(6, 6),
              numericInput("a_param", "Paramètre a:", 
                           value = 2, min = 1, max = 10, step = 1),
              numericInput("b_param", "Paramètre b:", 
                           value = 2, min = 1, max = 10, step = 1)
            )
          )
        )
      ),
      
      # Corrélations et charges
      accordion_panel(
        "Corrélations et charges",
        
        # Corrélations entre facteurs
        selectInput("fc_type", "Corrélation entre facteurs:",
                    choices = c("Fixe" = "fixed", "Aléatoire (0.0-0.5)" = "R"),
                    selected = "R"),
        
        conditionalPanel(
          condition = "input.fc_type == 'fixed'",
          sliderInput("fc_fixed", "Valeur de corrélation:", 
                      min = 0, max = 0.75, value = 0.3, step = 0.05)
        ),
        
        # Charges primaires
        selectInput("pl_type", "Charges factorielles primaires:",
                    choices = c("Fixe" = "fixed", 
                                "Faibles (0.35-0.50)" = "L",
                                "Moyennes (0.50-0.65)" = "M",
                                "Élevées (0.65-0.80)" = "H",
                                "Aléatoires (0.35-0.80)" = "R"),
                    selected = "R"),
        
        conditionalPanel(
          condition = "input.pl_type == 'fixed'",
          sliderInput("pl_fixed", "Valeur de charge primaire:", 
                      min = 0, max = 1, value = 0.6, step = 0.05)
        ),
        
        # Charges secondaires
        selectInput("cl_type", "Charges factorielles secondaires:",
                    choices = c("Fixe" = "fixed", 
                                "Faibles (-0.1-0.1)" = "L",
                                "Élevées (-0.2--0.1 ∪ 0.1-0.2)" = "H",
                                "Aucune (0)" = "None",
                                "Aléatoires (-0.2-0.2)" = "R"),
                    selected = "R"),
        
        conditionalPanel(
          condition = "input.cl_type == 'fixed'",
          sliderInput("cl_fixed", "Valeur de charge secondaire:", 
                      min = 0, max = 0.5, value = 0.1, step = 0.05)
        )
      )
    ),
    
    br(),
    
    # Zone d'information déroulante en bas
    accordion(
      open = FALSE,
      accordion_panel(
        "Aide",
        tags$div(
          style = "font-size: 0.9em;",
          HTML("
            <p><strong>Préréglages :</strong> Sélectionnez un type de simulation prédéfini puis cliquez sur Générer.</p>
            <p><strong>Structure des facteurs :</strong> Configurez le nombre de facteurs et de variables.</p>
            <p><strong>Distribution :</strong> Choisissez entre distribution normale ou asymétrique (Beta).</p>
            <p><strong>Corrélations :</strong> Définissez la force des relations entre facteurs latents.</p>
            <p><strong>Charges :</strong> Contrôlez la force des liens entre variables et facteurs.</p>
          ")
        )
      )
    )
  ),
  
  # Contenu principal
  navset_card_tab(
    nav_panel(
      title = "Introduction",
      icon = icon("info-circle"),
      
      # Introduction générale à l'analyse factorielle exploratoire
      card(
        card_header(
          div(
            class = "d-flex align-items-center",
            span(icon("book"), style = "margin-right: 10px;"),
            h3("Introduction à l'analyse factorielle exploratoire", class = "m-0")
          )
        ),
        card_body(
          p("L'analyse factorielle exploratoire (AFE) est une technique statistique puissante qui permet de réduire un ensemble de variables observées à un nombre réduit de facteurs latents. Elle est largement utilisée en sciences sociales, psychologie, marketing et dans de nombreuses autres disciplines pour explorer la structure sous-jacente des données et identifier des construits non directement observables."),
          
          h4("Concepts fondamentaux"),
          tags$ul(
            tags$li(strong("Facteurs latents"), " - Variables non directement observées qui influencent les variables mesurées. Ces facteurs représentent des construits théoriques sous-jacents."),
            tags$li(strong("Charges factorielles"), " - Coefficients qui indiquent la force de la relation entre les variables observées et les facteurs latents."),
            tags$li(strong("Communalité"), " - Proportion de la variance d'une variable qui est expliquée par les facteurs communs."),
            tags$li(strong("Unicité"), " - Proportion de la variance d'une variable qui n'est pas expliquée par les facteurs communs (spécificité + erreur de mesure)."),
            tags$li(strong("Rotation"), " - Procédure mathématique visant à obtenir une solution factorielle plus interprétable tout en préservant la variance expliquée.")
          ),
          
          h4("Étapes de l'analyse factorielle exploratoire"),
          tags$ol(
            tags$li(strong("Préparation des données"), " - Examen des statistiques descriptives, gestion des valeurs manquantes et des valeurs aberrantes."),
            tags$li(strong("Évaluation de la factorabilité"), " - Vérification que les données sont appropriées pour l'AFE (test de sphéricité de Bartlett, mesure KMO)."),
            tags$li(strong("Détermination du nombre de facteurs"), " - Utilisation de différentes méthodes (valeurs propres > 1, test parallèle, méthode Hull) pour identifier le nombre optimal de facteurs."),
            tags$li(strong("Extraction des facteurs"), " - Application d'une méthode d'extraction comme l'analyse en composantes principales, la factorisation en axes principaux, ou le maximum de vraisemblance."),
            tags$li(strong("Rotation des facteurs"), " - Application d'une rotation orthogonale (varimax, quartimax) ou oblique (promax, oblimin) pour simplifier l'interprétation."),
            tags$li(strong("Interprétation des résultats"), " - Analyse des charges factorielles pour donner du sens aux facteurs identifiés."),
            tags$li(strong("Validation de la solution"), " - Évaluation de la cohérence interne, calcul des scores factoriels, vérification de la stabilité de la solution.")
          )
        )
      ),
      
      # Comment utiliser cette application
      card(
        card_header(
          div(
            class = "d-flex align-items-center",
            span(icon("hand-point-right"), style = "margin-right: 10px;"),
            h3("Utilisation de cette application", class = "m-0")
          )
        ),
        card_body(
          p("Cette application vous permet de simuler des données pour explorer et comprendre les concepts de l'analyse factorielle exploratoire. Elle offre une approche pratique et interactive pour se familiariser avec les différentes étapes et méthodes."),
          
          h4("Comment procéder"),
          tags$ol(
            tags$li(strong("Sélectionnez un préréglage"), " ou configurez manuellement les paramètres de simulation dans le panneau latéral."),
            tags$li(strong("Générez les données"), " en cliquant sur le bouton correspondant."),
            tags$li(strong("Explorez les données simulées"), " dans l'onglet 'Données et statistiques descriptives'."),
            tags$li(strong("Effectuez des tests préliminaires"), " pour évaluer la factorabilité et déterminer le nombre optimal de facteurs."),
            tags$li(strong("Réalisez l'analyse factorielle"), " en configurant les paramètres d'extraction et de rotation."),
            tags$li(strong("Interprétez les résultats"), " à l'aide des différentes visualisations et indicateurs fournis.")
          ),
          
          h4("Points forts de l'application"),
          tags$ul(
            tags$li("Simulation de données avec une structure factorielle connue pour faciliter l'apprentissage"),
            tags$li("Implémentation de multiples méthodes de détermination du nombre de facteurs"),
            tags$li("Variété d'options pour l'extraction et la rotation des facteurs"),
            tags$li("Visualisations interactives pour faciliter l'interprétation"),
            tags$li("Évaluation de la qualité des solutions factorielles")
          )
        )
      ),

      # Références bibliographiques
      card(
        card_header(
          div(
            class = "d-flex align-items-center",
            span(icon("book"), style = "margin-right: 10px;"),
            h3("Références et ressources", class = "m-0")
          )
        ),
        card_body(
          h4("Références bibliographiques"),
          tags$ul(
            tags$li("Bartholomew, D. J., Knott, M., & Moustaki, I. (2011). ", tags$em("Latent variable models and factor analysis: A unified approach"), ". John Wiley & Sons."),
            tags$li("Fabrigar, L. R., & Wegener, D. T. (2011). ", tags$em("Exploratory factor analysis"), ". Oxford University Press."),
            tags$li("Lorenzo-Seva, U., & Ferrando, P. J. (2006). FACTOR: A computer program to fit the exploratory factor analysis model. ", tags$em("Behavior Research Methods"), ", 38(1), 88-91."),
            tags$li("Tabachnick, B. G., & Fidell, L. S. (2019). ", tags$em("Using multivariate statistics"), " (7th ed.). Pearson."),
            tags$li("Thompson, B. (2004). ", tags$em("Exploratory and confirmatory factor analysis: Understanding concepts and applications"), ". American Psychological Association.")
          ),
          
          h4("Ressources en ligne"),
          tags$ul(
            tags$li(tags$a(href = "https://psych.wisc.edu/henriques/factor.html", target = "_blank", "An Introduction to Factor Analysis - University of Wisconsin-Madison")),
            tags$li(tags$a(href = "https://stats.idre.ucla.edu/spss/seminars/introduction-to-factor-analysis/", target = "_blank", "Introduction to Factor Analysis - UCLA")),
            tags$li(tags$a(href = "https://personality-project.org/r/psych/", target = "_blank", "psych Package for R - Documentation and Tutorials")),
            tags$li(tags$a(href = "https://cran.r-project.org/web/packages/EFA.dimensions/index.html", target = "_blank", "EFA.dimensions Package for R - Documentation")),
            tags$li(tags$a(href = "https://www.statology.org/exploratory-factor-analysis-r/", target = "_blank", "Guide to Exploratory Factor Analysis in R - Statology"))
          ),
          
          h4("Conseils pratiques"),
          tags$ul(
            tags$li("Utilisez plusieurs méthodes pour déterminer le nombre de facteurs et comparez les résultats."),
            tags$li("Tenez compte à la fois des résultats statistiques et de la pertinence théorique lors de l'interprétation."),
            tags$li("Pour l'interprétation des facteurs, concentrez-vous sur les variables ayant des charges supérieures à 0.3 ou 0.4 en valeur absolue."),
            tags$li("Soyez attentif aux variables à faible communalité ou à celles qui chargent fortement sur plusieurs facteurs."),
            tags$li("Le choix entre rotation orthogonale et oblique dépend de vos hypothèses concernant la corrélation entre les facteurs.")
          )
        )
      )
    ),
    
    # Section "Données et statistiques descriptives" uniquement

    nav_panel(
      title = "Données et statistiques descriptives",
      icon = icon("table"),

      # Navset principal pour les onglets
      navset_card_tab(
        # 1. Onglet Données simulées
        nav_panel(
          title = "Données simulées",
          icon = icon("table"),
          card(
            full_screen = TRUE,
            card_header("Aperçu des données de réponse simulées"),
            card_body(
              div(
                class = "d-flex justify-content-end mb-3",
                downloadButton("download_data", "Télécharger les données", class = "btn-sm")
              ),
              DTOutput("simulated_data")
            ),
            card_footer("Ces données représentent les réponses simulées aux items du questionnaire.")
          )
        ),
        
        # 2. Onglet Statistiques univariées
        nav_panel(
          title = "Statistiques univariées",
          icon = icon("chart-simple"),
          
          navset_pill(
            # A. Sous-onglet Statistiques descriptives
            nav_panel(
              title = "Statistiques descriptives",
              card(
                full_screen = TRUE,
                height = "800px",
                card_header("Statistiques descriptives complètes"),
                card_body(
                  DTOutput("descriptive_stats_table")
                )
              )
            ),
            
            # B. Sous-onglet Distributions
            nav_panel(
              title = "Distributions",
              layout_columns(
                col_widths = c(10, 2),
                card(
                  full_screen = TRUE,
                  card_header("Visualisation des distributions"),
                  card_body(
                    plotOutput("distribution_plot", height = "500px")
                  )
                ),
                card(
                  card_header("Options de visualisation"),
                  card_body(
                    selectInput("distribution_plot_type", "Type de graphique:",
                                choices = c("Histogrammes" = "hist",
                                            "Boîtes à moustaches" = "boxplot",
                                            "Courbes de densité" = "density",
                                            "QQ-plots" = "qqplot"),
                                selected = "hist"),
                    uiOutput("distrib_var_selector"),
                    hr(),
                    downloadButton("download_distrib_plot", "Télécharger le graphique", class = "btn-sm btn-secondary")
                  )
                )
              )
            ),
            
            # C. Sous-onglet Tests de normalité
            nav_panel(
              title = "Tests de normalité",
              card(
                full_screen = TRUE,
                card_header("Résultats des tests de normalité"),
                card_body(
                  p("Ce tableau présente les résultats de tous les tests de normalité disponibles pour chaque variable."),
                  p("La conclusion est basée sur l'ensemble des tests : une variable est considérée comme normale si la majorité des tests l'indique."),
                  br(),
                  DTOutput("normality_results_table")
                )
              )
            )
          )
        ),
        
        # 3. Onglet Statistiques bivariées
        nav_panel(
          title = "Statistiques bivariées",
          icon = icon("chart-line"),
          
          navset_pill(
            # A. Sous-onglet Corrélations - avec groupCheckbox pour la sélection des variables
            nav_panel(
              title = "Corrélations",
              layout_columns(
                col_widths = c(9, 3),
                card(
                  full_screen = TRUE,
                  card_header("Matrice de corrélation"),
                  card_body(
                    plotOutput("item_corr_plot", height = "600px")
                  )
                ),
                card(
                  card_header("Options"),
                  card_body(
                    selectInput("corr_method", "Méthode de corrélation:",
                                choices = c("Pearson" = "pearson", 
                                            "Spearman" = "spearman", 
                                            "Kendall" = "kendall"),
                                selected = "pearson"),
                    
                    # Interface pour sélectionner les variables (similaire aux distributions)
                    p("Variables à inclure:"),
                    uiOutput("corr_var_selector"),
                    hr(),
                    
                    # Options de visualisation améliorées
                    checkboxInput("show_corr_values", "Afficher les valeurs numériques", TRUE),
                    selectInput("corr_plot_method", "Style de visualisation:",
                                choices = c("Couleurs" = "color", 
                                            "Cercles" = "circle", 
                                            "Carrés" = "square",
                                            "Ellipses" = "ellipse",
                                            "Ombrages" = "shade"),
                                selected = "color"),
                    selectInput("corr_plot_order", "Ordre des variables:",
                                choices = c("Original" = "original", 
                                            "Par groupes" = "hclust", 
                                            "Par valeurs" = "FPC"),
                                selected = "original"),
                    uiOutput("clust_method_ui"),
                    selectInput("corr_plot_type", "Orientation de la matrice:",
                                choices = c("Triangulaire supérieure" = "upper", 
                                            "Triangulaire inférieure" = "lower", 
                                            "Matrice complète" = "full"),
                                selected = "upper"),
                    hr(),
                    downloadButton("download_corr_plot", "Télécharger le graphique", class = "btn-sm btn-secondary"),
                    downloadButton("download_corr_matrix", "Télécharger la matrice", class = "btn-sm btn-secondary", style = "margin-top: 10px;")
                  )
                )
              )
            ),
            
            # B. Sous-onglet Nuages de points - également avec groupCheckbox
            nav_panel(
              title = "Nuages de points",
              layout_columns(
                col_widths = c(9, 3),
                card(
                  full_screen = TRUE,
                  card_header("Matrice de nuages de points"),
                  card_body(
                    plotOutput("pairs_plot", height = "600px")
                  )
                ),
                card(
                  card_header("Options"),
                  card_body(
                    p("Variables à inclure:"),
                    uiOutput("pairs_var_selector"),

                    # Options de visualisation améliorées
                    checkboxInput("show_ellipses", "Afficher les ellipses de confiance", TRUE),
                    checkboxInput("show_smooth", "Ajouter courbe de lissage", TRUE),
                    selectInput("pairs_color_scheme", "Palette de couleurs:",
                                choices = c("Bleu" = "blue", 
                                            "Rouge" = "red", 
                                            "Vert" = "green",
                                            "Violet" = "purple"),
                                selected = "blue"),
                    sliderInput("pairs_alpha", "Transparence:", 
                                min = 0.1, max = 1, value = 0.5, step = 0.1),
                    hr(),
                    downloadButton("download_pairs_plot", "Télécharger le graphique", class = "btn-sm btn-secondary")
                  )
                )
              )
            )
          )
        )
      )
    ),
    
    nav_panel(
      title = "Tests préliminaires",
      icon = icon("check-double"),
      navset_card_tab(
        # Premier sous-onglet: Tests factoriels
        nav_panel(
          title = "Tests factoriels",
          layout_columns(
            col_widths = c(3, 9),
            # Card des inputs à gauche
            card(
              card_header("Paramètres des tests"),
              card_body(
                checkboxGroupInput("dimtests", tags$h5("Tests à effectuer:"),
                                   choices = c("Coordonnées optimales" = "CD",
                                               "Critère d'information de Kaiser" = "EMPKC",
                                               "Méthode Hull" = "HULL",
                                               "Minimum Average Partial" = "MAP",
                                               "Valeurs propres > 1" = "NEVALSGT1",
                                               "Scree test" = "SESCREE",
                                               "Test de dimensionnalité séquentielle" = "SMT",
                                               "Test Parallèle" = "PARALLEL"),
                                   
                                   selected = c("EMPKC", "HULL", "PARALLEL")),
                
                selectInput("dimtests_cor", "Type de corrélation:", 
                            choices = c("Pearson" = "pearson", 
                                        "Kendall" = "kendall", 
                                        "Spearman" = "spearman",
                                        "Gamma" = "gamma",
                                        "Polychorique" = "polychoric"),
                            selected = "pearson",
                            width = "150px"),
                
                conditionalPanel(
                  condition = "input.dimtests.includes('HULL')",
                  tags$h5("Paramètres de la méthode HULL"),
                  radioButtons("dimtests_hull_method", "Méthode d'estimation:", 
                               choices = c("PAF" = "PAF", 
                                           "ULS" = "ULS", 
                                           "ML" = "ML"),
                               selected = "PAF",
                               inline = TRUE),
                  radioButtons("dimtests_hull_gof", "Indice d'ajustement:", 
                               choices = c("CAF" = "CAF", 
                                           "CFI" = "CFI", 
                                           "RMSEA" = "RMSEA"),
                               selected = "CAF",
                               inline = TRUE)
                ),
                
                conditionalPanel(
                  condition = "input.dimtests.includes('PARALLEL')",
                  tags$h5("Paramètres du test parallèle"),
                  
                  radioButtons("rawpar_randtype", "Type de données aléatoires",
                               choices = c("Générées" = "generated", "Permutées" = "permuted"),
                               selected = "generated",
                               inline = TRUE),
                  
                  radioButtons("rawpar_extraction", "Méthode d'extraction",
                               choices = c("PAF", 
                                           "PCA",
                                           "image"),
                               selected = "PAF",
                               inline = TRUE),
                  
                  selectInput("rawpar_corkindRAND", "Type de corrélation",
                               choices = c("Pearson" = "pearson", 
                                           "Kendall" = "kendall", 
                                           "Spearman" = "spearman",
                                           "Gamma" = "gamma",
                                           "Polychorique" = "polychoric"),
                               selected = "pearson",
                              width = "240px"),
                  fluidRow(
                    column(6,
                           numericInput("rawpar_ndatasets", "# Itérations",
                                         value = 100, min = 10, max = 1000, step = 10,
                                         width = "100px")
                  ),
                  column(6,
                         numericInput("rawpar_percentile", "Percentile",
                                       value = 95, min = 50, max = 99, step = 1,
                                       width = "100px")
                  )
                )
                ),
                
                actionButton("run_dimtests", "Exécuter les tests", 
                             icon = icon("calculator"),
                             class = "btn-primary", 
                             style = "margin-top: 15px; width: 100%;")
              )
            
            ),
            
            # Zone de résultats à droite
            card(
              card_header("Résultats des tests dimensionnels"),
              navset_card_tab(
                nav_panel(
                  title = "Résumé",
                  htmlOutput("dimtests_summary")
                ),
                nav_panel(
                  title = "Scree Plot",
                  plotlyOutput("dimtests_plot", height = "600px", width = "100%")
                ),
                nav_panel(
                  title = "Détails techniques",
                  verbatimTextOutput("dimtests_results")
                )
              )
            )
          )
        ),
        
        # Deuxième sous-onglet: Factorabilité
        nav_panel(
          title = "Factorabilité",
          icon = icon("check-circle"),
          layout_columns(
            col_widths = c(2, 10),
            # Card des inputs à gauche
            card(
              card_header("Paramètres de factorabilité"),
              card_body(
                radioButtons("factor_cor", "Type de corrélation:", 
                             choices = c("Pearson" = "pearson", 
                                         "Kendall" = "kendall", 
                                         "Spearman" = "spearman"),
                             selected = "pearson",
                             width = "150px"),
#                actionButton("run_factorability", "Exécuter l'analyse", 
#                             icon = icon("calculator"),
#                             class = "btn-primary", 
#                             style = "margin-top: 15px; width: 100%;")
              )
            ),
            
            # Zone de résultats à droite
            card(
              card_header("Résultats de l'analyse de factorabilité"),
              verbatimTextOutput("factorability_raw_results")
            )
          )
        )
      )
    
  ),
    
nav_panel(
  title = "Analyse factorielle",
  icon = icon("project-diagram"),
  layout_columns(
    col_widths = c(3, 9),
    # Card des inputs à gauche
    card(
      card_header("Paramètres de l'analyse factorielle"),
      card_body(
        radioButtons("analysis_type", "Type d'analyse:", 
                     choices = c("EFA" = "efa",
                                 "PCA" = "pca"),
                     selected = "efa",
                     inline = TRUE),
        selectInput("extraction_method", "Méthode d'extraction:", 
                    choices = c("PAF" = "paf",
                                "ML" = "ml",
                                "Image" = "image",
                                "MinRes" = "minres",
                                "ULS" = "uls",
                                "OLS" = "ols",
                                "WLS" = "wls",
                                "GLS" = "gls",
                                "Alpha" = "alpha",
                                "FullInfo" = "fullinfo"),
                    selected = "paf"),
        
        selectInput("factor_method", "Choix du nombre de facteurs:",
                    choices = c("Automatique (EMPKC)" = "auto",
                                "Spécifier manuellement" = "manual"),
                    selected = "auto"),
        
        conditionalPanel(
          condition = "input.factor_method == 'manual'",
          numericInput("num_factors", "Nombre de facteurs:", 
                       value = 3, min = 1, max = 10, step = 1)
        ),
        
        # Options de rotation
        radioButtons("rotation_type", "Type de rotation:", 
                     choices = c("Orthogonal" = "orthogonal", 
                                 "Oblique" = "oblique",
                                 "Aucune" = "none"),
                     selected = "oblique",
                     inline = TRUE),
        
        conditionalPanel(
          condition = "input.rotation_type == 'orthogonal'",
          selectInput("rotation_method", "Méthode de rotation orthogonale:", 
                      choices = c("Varimax" = "varimax",
                                  "Quartimax" = "quartimax",
                                  "BentlerT" = "bentlerT",
                                  "Equamax" = "equamax",
                                  "GeominT" = "geominT",
                                  "BifactorT" = "bifactorT",
                                  "Entropy" = "entropy"),
                      selected = "varimax")
        ),
        
        conditionalPanel(
          condition = "input.rotation_type == 'oblique'",
          selectInput("rotation_method", "Méthode de rotation oblique:", 
                      choices = c("Promax" = "promax",
                                  "Quartimin" = "quartimin",
                                  "Oblimin" = "oblimin",
                                  "Oblimax" = "oblimax",
                                  "Simplimax" = "simplimax",
                                  "BentlerQ" = "bentlerQ",
                                  "GeominQ" = "geominQ",
                                  "BifactorQ" = "bifactorQ"),
                      selected = "promax")
        ),
        
        # Type de corrélation
        selectInput("efa_cor", "Type de corrélation:", 
                    choices = c("Pearson" = "pearson", 
                                "Kendall" = "kendall", 
                                "Spearman" = "spearman",
                                "Gamma" = "gamma",
                                "Polychorique" = "polychoric"),
                    selected = "pearson"),
        sliderInput("loading_threshold", "Seuil minimal des saturations à afficher:", 
                    min = 0.0, max = 0.7, value = 0.3, step = 0.05),
        
        actionButton("run_efa", "Exécuter l'analyse factorielle", 
                     icon = icon("calculator"),
                     class = "btn-primary", 
                     style = "margin-top: 15px; width: 100%")
      )
    ),
    
    # Zone de résultats à droite
    card(
      card_header(
        uiOutput("analysis_result_title") 
      ),
      navset_card_tab(
        nav_panel(
          title = "Résultats",
          verbatimTextOutput("efa_results")
        ),
        nav_panel(
          title = "Biplot",
          fluidRow(
            column(12,
                   uiOutput("factor_selector_ui"),
                   plotlyOutput("efa_biplot")
            )
          )
        ),
        nav_panel(
          title = "Complexité",
          verbatimTextOutput("complexity_results")
        ),
        nav_panel(
          title = "Cohérence interne",
          icon = icon("sync"),
          layout_columns(
            col_widths = c(4, 8),
            # Card des inputs à gauche
            card(
              card_header("Paramètres de cohérence interne"),
              card_body(
                selectInput("ic_extraction", "Méthode d'extraction:", 
                            choices = c("ML", 
                                        "PAF"),
                            selected = "ML"),
                actionButton("run_internal_consistency", "Exécuter l'analyse", 
                             icon = icon("calculator"),
                             class = "btn-primary", 
                             style = "margin-top: 15px; width: 100%;")
              )
            ),
            
            # Zone de résultats à droite
            card(
              card_header("Résultats de l'analyse de cohérence interne"),
              verbatimTextOutput("internal_consistency_results")
            )
          )
        ),
        # Remplacez la définition de l'onglet Analyse des résidus par celui-ci:
        nav_panel(
          title = "Analyse des résidus",
          value = "residuals_tab",
          fluidRow(
            column(
              width = 12,
              card(
                card_header("Matrice des corrélations résiduelles"),
                card_body(
                  p("Cette section présente l'analyse des résidus du modèle factoriel, permettant d'évaluer la qualité d'ajustement du modèle aux données."),
                  p("La matrice des résidus est calculée comme la différence entre la matrice de corrélations observée et la matrice de corrélations reproduite par le modèle factoriel."),
                  
                  # Message de diagnostic visible
                  verbatimTextOutput("residuals_debug_info"),
                  
                  # Indicateur de statut du calcul
                  textOutput("residuals_status"),
                  
                  fluidRow(
                    column(
                      width = 6,
                      h5("Indicateurs d'ajustement"),
                      uiOutput("residual_stats")
                    ),
                    column(
                      width = 6,
                      h5("Interprétation"),
                      htmlOutput("residual_interpretation")
                    )
                  ),
                  
                  hr(),
                  
                  h5("Matrice des résidus"),
                  p("Les cellules sont colorées selon l'amplitude des résidus. Les résidus importants (>0.1) sont surlignés."),
                  DTOutput("residual_matrix_table"),
                  
                  tags$small(
                    "Note : Des résidus élevés indiquent des relations entre variables que le modèle ne capture pas bien. 
            Un modèle bien ajusté devrait avoir peu de résidus supérieurs à 0.05 en valeur absolue."
                  )
                )
              )
            )
          )
        ),
        nav_panel(
          title = "Scores factoriels",
          layout_columns(
            col_widths = c(4, 8),
            # Card des inputs à gauche
            card(
              card_header("Paramètres des scores factoriels"),
              card_body(
                selectInput("factor_scores_method", "Méthode de calcul:", 
                            choices = c("Thurstone (regression)" = "Thurstone", 
                                        "Harman (idealized)" = "Harman",
                                        "Bartlett" = "Bartlett",
                                        "ten Berge" = "tenBerge",
                                        "Anderson & Rubin" = "Anderson",
                                        "PCA (non rotationnels)" = "PCA"),
                            selected = "Thurstone"),
                actionButton("run_factor_scores", "Calculer les scores factoriels", 
                             icon = icon("calculator"),
                             class = "btn-primary", 
                             style = "margin-top: 15px; width: 100%;")
              )
            ),
            
            # Zone de résultats à droite
            card(
              card_header("Résultats des scores factoriels"),
              navset_card_tab(
                nav_panel(
                  title = "Coefficients et statistiques",
                  verbatimTextOutput("factor_scores_stats")
                ),
                nav_panel(
                  title = "Scores factoriels",
                  downloadButton("download_factor_scores", 
                                 "Télécharger les scores", 
                                 style = "margin-bottom: 10px;"),
                  
                  DTOutput("factor_scores_table")
                )
              )
            )
          )
        ),
        nav_panel(
          title = "Comparaison de modèles",
          icon = icon("balance-scale"),
          layout_columns(
            col_widths = c(3, 9),
            # Card des inputs à gauche
            card(
              card_header("Paramètres de comparaison"),
              card_body(
                numericInput("compare_min_factors", "Nombre minimum de facteurs:", 
                             value = 1, min = 1, max = 10, step = 1),
                numericInput("compare_max_factors", "Nombre maximum de facteurs:", 
                             value = 5, min = 2, max = 10, step = 1),
                selectInput("compare_extraction", "Méthode d'extraction:", 
                            choices = c("PAF" = "paf",
                                        "ML" = "ml",
                                        "MinRes" = "minres",
                                        "ULS" = "uls"),
                            selected = "paf"),
                selectInput("compare_rotation", "Méthode de rotation:", 
                            choices = c("Promax" = "promax",
                                        "Varimax" = "varimax",
                                        "Oblimin" = "oblimin",
                                        "Aucune" = "none"),
                            selected = "promax"),
                selectInput("compare_corkind", "Type de corrélation:", 
                            choices = c("Pearson" = "pearson", 
                                        "Kendall" = "kendall", 
                                        "Spearman" = "spearman"),
                            selected = "pearson"),
                hr(),
                checkboxGroupInput("compare_metrics", "Indices à afficher:",
                                   choices = c("AIC" = "AIC", 
                                               "BIC" = "BIC",
                                               "RMSEA" = "RMSEA",
                                               "CFI" = "CFI", 
                                               "TLI" = "TLI",
                                               "RMSR" = "RMSR",
                                               "Variance expliquée" = "Variance_Explained"),
                                   selected = c("AIC", "BIC", "CFI", "RMSEA")),
                
                actionButton("run_compare", "Comparer les modèles", 
                             icon = icon("balance-scale"),
                             class = "btn-primary", 
                             style = "margin-top: 15px; width: 100%;")
              )
            ),
            
            # Zone de résultats à droite
            card(
              card_header("Résultats de la comparaison de modèles"),
              navset_card_tab(
                nav_panel(
                  title = "Tableau comparatif",
                  DTOutput("compare_table")
                ),
                nav_panel(
                  title = "Graphiques",
                  plotlyOutput("compare_plots", height = "600px")
                ),
                nav_panel(
                  title = "Interprétation",
                  verbatimTextOutput("compare_interpretation")
                )
              )
            )
          )
        )
      )
    )
  )
),
nav_panel(
  title = "Rapport Final",
  icon = icon("file-alt"),
  layout_columns(
    col_widths = c(3, 9),
    # Card des inputs à gauche
    card(
      card_header("Configuration du rapport"),
      card_body(
        p("Générez un rapport synthétisant l'ensemble des analyses effectuées."),
        selectInput("rapport_nfactors", "Nombre de facteurs retenus:", 
                    choices = 1:10, selected = 3),
        selectInput("rapport_structure", "Structure factorielle:", 
                    choices = c("Simple" = "simple", 
                                "Complexe" = "complex"),
                    selected = "simple"),
        checkboxInput("rapport_include_communalities", "Inclure l'analyse des communalités", TRUE),
        checkboxInput("rapport_include_reliability", "Inclure l'analyse de fiabilité", TRUE),
        checkboxInput("rapport_include_methodology", "Inclure la méthodologie", TRUE),
        actionButton("generate_report", "Générer le rapport", 
                     icon = icon("file-alt"),
                     class = "btn-primary", 
                     style = "margin-top: 15px; width: 100%;")
      )
    ),
    
    # Zone de rapport à droite
    card(
      full_screen = TRUE,
      card_header(
        div(
          style = "display: flex; justify-content: space-between; align-items: center;",
          h3("Rapport d'analyse factorielle", class = "m-0"),
          downloadButton("download_report", "Télécharger le rapport", class = "btn-sm btn-secondary")
        )
      ),
      card_body(
        uiOutput("rapport_contenu")
      )
    )
  )
)
)
)



server <- function(input, output, session) {
  # Réactive pour stocker les données simulées
  sim_data <- reactiveVal(NULL)
  
  # Réactive pour stocker les résultats des tests dimensionnels
  dimtests_results <- reactiveVal(NULL)
  
  # Réactive pour stocker les résultats de factorabilité
  factorability_results <- reactiveVal(NULL)
  
  # Réactive pour stocker les résultats de cohérence interne
  internal_consistency_results <- reactiveVal(NULL)
  
  # Réactive pour stocker les résultats de l'analyse factorielle
  efa_result <- reactiveVal(NULL)
  
  # Réactive pour stocker les scores factoriels
  factor_scores_result <- reactiveVal(NULL)
  
  # Réactive pour stocker les résultats de comparaison
  compare_results <- reactiveVal(NULL)

  # Fonction réactive pour calculer les statistiques descriptives (à placer dans la partie serveur)
  descriptive_stats_reactive <- reactive({
    req(sim_data())
    
    # S'assurer que le package moments est disponible
    if (!requireNamespace("moments", quietly = TRUE)) {
      warning("Le package 'moments' n'est pas installé. Les mesures d'asymétrie et kurtosis ne seront pas disponibles.")
      has_moments <- FALSE
    } else {
      has_moments <- TRUE
      library(moments)
    }
    
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Obtenir toutes les variables
    all_vars <- colnames(response_data)
    
    # Créer le dataframe qui contiendra les résultats
    results <- data.frame(
      Variable = all_vars,
      N = sapply(all_vars, function(var) sum(!is.na(response_data[, var]))),
      Min = sapply(all_vars, function(var) min(response_data[, var], na.rm = TRUE)),
      Q1 = sapply(all_vars, function(var) quantile(response_data[, var], 0.25, na.rm = TRUE)),
      Médiane = sapply(all_vars, function(var) median(response_data[, var], na.rm = TRUE)),
      Moyenne = sapply(all_vars, function(var) mean(response_data[, var], na.rm = TRUE)),
      Q3 = sapply(all_vars, function(var) quantile(response_data[, var], 0.75, na.rm = TRUE)),
      Max = sapply(all_vars, function(var) max(response_data[, var], na.rm = TRUE)),
      Écart_type = sapply(all_vars, function(var) sd(response_data[, var], na.rm = TRUE)),
      CV = sapply(all_vars, function(var) {
        m <- mean(response_data[, var], na.rm = TRUE)
        if (m == 0) return(NA) else return(sd(response_data[, var], na.rm = TRUE) / m)
      }),
      stringsAsFactors = FALSE
    )
    
    # Ajouter les statistiques avancées si le package moments est disponible
    if (has_moments) {
      results$Asymétrie <- sapply(all_vars, function(var) moments::skewness(response_data[, var], na.rm = TRUE))
      results$Kurtosis <- sapply(all_vars, function(var) moments::kurtosis(response_data[, var], na.rm = TRUE))
      results$Kurtosis_excess <- results$Kurtosis - 3  # Kurtosis excédentaire (0 pour une loi normale)
    }
    
    # Ajouter l'intervalle de confiance à 95% pour la moyenne
    results$IC_95_inf <- sapply(all_vars, function(var) {
      x <- response_data[, var]
      se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      mean(x, na.rm = TRUE) - qt(0.975, sum(!is.na(x)) - 1) * se
    })
    
    results$IC_95_sup <- sapply(all_vars, function(var) {
      x <- response_data[, var]
      se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
      mean(x, na.rm = TRUE) + qt(0.975, sum(!is.na(x)) - 1) * se
    })
    
    # Calculer l'étendue interquartile (IQR)
    results$IQR <- results$Q3 - results$Q1
    
    # Ajouter le mode (valeur la plus fréquente)
    results$Mode <- sapply(all_vars, function(var) {
      x <- response_data[, var]
      # Arrondir pour gérer les variables continues
      x_rounded <- round(x, 3)
      tab <- table(x_rounded)
      as.numeric(names(tab)[which.max(tab)])
    })
    
    return(results)
  })
  
  # Rendu du tableau des statistiques descriptives complet et élégant
  output$descriptive_stats_table <- renderDT({
    req(descriptive_stats_reactive())
    
    results <- descriptive_stats_reactive()
    
    # Renommer certaines colonnes pour une meilleure présentation
    col_names <- colnames(results)
    col_names[col_names == "IC_95_inf"] <- "IC 95% (inf)"
    col_names[col_names == "IC_95_sup"] <- "IC 95% (sup)"
    col_names[col_names == "Kurtosis_excess"] <- "Kurtosis excess"
    colnames(results) <- col_names
    
    # Création d'un tableau plus élégant avec formatage personnalisé
    dt <- datatable(
      results,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'tip',
        columnDefs = list(
          list(className = 'dt-center', targets = 1:ncol(results)-1) # Centrer toutes les colonnes sauf la première
        )
      ),
      caption = "Statistiques descriptives complètes"
    )
    
    # Formater les différentes colonnes selon leur type
    # Nombre entier
    dt <- formatRound(dt, columns = "N", digits = 0)
    
    # Nombre avec 3 décimales pour la plupart des statistiques
    decimal_cols <- which(colnames(results) %in% c("Min", "Q1", "Médiane", "Moyenne", 
                                                   "Q3", "Max", "Écart_type", "IQR", "Mode"))
    if (length(decimal_cols) > 0) {
      dt <- formatRound(dt, columns = decimal_cols, digits = 3)
    }
    
    # Nombre avec 4 décimales pour les statistiques plus précises
    precise_cols <- which(colnames(results) %in% c("CV", "Asymétrie", "Kurtosis", 
                                                   "Kurtosis excess", "IC 95% (inf)", "IC 95% (sup)"))
    if (length(precise_cols) > 0) {
      dt <- formatRound(dt, columns = precise_cols, digits = 4)
    }
    
    # Appliquer un style conditionnel pour l'asymétrie
    if ("Asymétrie" %in% colnames(results)) {
      dt <- formatStyle(
        dt,
        'Asymétrie',
        backgroundColor = styleInterval(
          c(-0.5, 0.5),
          c('rgba(173, 216, 230, 0.5)', 'rgba(255, 255, 255, 0)', 'rgba(255, 182, 193, 0.5)')
        ),
        color = styleInterval(
          c(-2, -1, 1, 2),
          c('darkred', 'darkblue', 'black', 'darkblue', 'darkred')
        )
      )
    }
    
    # Appliquer un style conditionnel pour le kurtosis excess
    if ("Kurtosis excess" %in% colnames(results)) {
      dt <- formatStyle(
        dt,
        'Kurtosis excess',
        backgroundColor = styleInterval(
          c(-1, 1),
          c('rgba(173, 216, 230, 0.5)', 'rgba(255, 255, 255, 0)', 'rgba(255, 182, 193, 0.5)')
        ),
        color = styleInterval(
          c(-3, -2, 2, 3),
          c('darkred', 'darkblue', 'black', 'darkblue', 'darkred')
        )
      )
    }
    
    return(dt)
  })
  
  # Variable réactive pour le calcul des tests de normalité (à placer dans la section serveur)
  normality_results_reactive <- reactive({
    req(sim_data())
    
    # S'assurer que le package nortest est disponible pour les tests d'Anderson-Darling
    if (!requireNamespace("nortest", quietly = TRUE)) {
      warning("Le package 'nortest' n'est pas installé. Le test Anderson-Darling ne sera pas disponible.")
      has_nortest <- FALSE
    } else {
      has_nortest <- TRUE
    }
    
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Obtenir toutes les variables
    all_vars <- colnames(response_data)
    
    # Créer le dataframe avec uniquement le nom des variables
    results <- data.frame(
      Variable = all_vars,
      stringsAsFactors = FALSE
    )
    
    # Ajouter les colonnes pour les tests et statistiques
    results$Shapiro_Stat <- NA
    results$Shapiro_p <- NA
    if (has_nortest) {
      results$Anderson_Stat <- NA
      results$Anderson_p <- NA
    }
    results$KS_Stat <- NA
    results$KS_p <- NA
    results$Conclusion <- NA
    
    # Calculer les tests pour chaque variable
    for (i in 1:length(all_vars)) {
      var_name <- all_vars[i]
      var_data <- response_data[, var_name]
      
      # Test de Shapiro-Wilk
      shapiro_test <- tryCatch({
        test <- shapiro.test(var_data)
        results$Shapiro_Stat[i] <- test$statistic
        results$Shapiro_p[i] <- test$p.value
      }, error = function(e) {
        results$Shapiro_Stat[i] <- NA
        results$Shapiro_p[i] <- NA
      })
      
      # Test d'Anderson-Darling si disponible
      if (has_nortest) {
        ad_test <- tryCatch({
          test <- nortest::ad.test(var_data)
          results$Anderson_Stat[i] <- test$statistic
          results$Anderson_p[i] <- test$p.value
        }, error = function(e) {
          results$Anderson_Stat[i] <- NA
          results$Anderson_p[i] <- NA
        })
      }
      
      # Test de Kolmogorov-Smirnov
      ks_test <- tryCatch({
        test <- ks.test(var_data, "pnorm", mean = mean(var_data), sd = sd(var_data))
        results$KS_Stat[i] <- test$statistic
        results$KS_p[i] <- test$p.value
      }, error = function(e) {
        results$KS_Stat[i] <- NA
        results$KS_p[i] <- NA
      })
    }
    
    # Déterminer la conclusion basée sur tous les tests disponibles
    results$Conclusion <- apply(results, 1, function(row) {
      # Récupérer les p-values des tests de normalité
      p_values <- c()
      
      # Shapiro-Wilk
      shapiro_idx <- which(names(row) == "Shapiro_p")
      if (!is.na(row[shapiro_idx])) {
        p_values <- c(p_values, row[shapiro_idx])
      }
      
      # Anderson-Darling (si disponible)
      anderson_idx <- which(names(row) == "Anderson_p")
      if (length(anderson_idx) > 0 && !is.na(row[anderson_idx])) {
        p_values <- c(p_values, row[anderson_idx])
      }
      
      # Kolmogorov-Smirnov
      ks_idx <- which(names(row) == "KS_p")
      if (!is.na(row[ks_idx])) {
        p_values <- c(p_values, row[ks_idx])
      }
      
      # Calculer la conclusion basée sur les p-values disponibles
      if (length(p_values) == 0) {
        return("Indéterminé")
      }
      
      # Compter les tests qui suggèrent la normalité (p > 0.05)
      normal_count <- sum(p_values > 0.05)
      non_normal_count <- sum(p_values <= 0.05)
      
      if (normal_count > non_normal_count) {
        return("Normale")
      } else if (non_normal_count > normal_count) {
        return("Non-normale")
      } else {
        return("Mitigé")
      }
    })
    
    return(results)
  })
  
  # Rendu du tableau des résultats - sans les statistiques descriptives et sans la boîte de recherche
  output$normality_results_table <- renderDT({
    req(normality_results_reactive())
    
    results <- normality_results_reactive()
    
    # Renommer les colonnes pour une meilleure présentation
    col_names <- colnames(results)
    col_names[col_names == "Shapiro_Stat"] <- "Statistique W (Shapiro)"
    col_names[col_names == "Shapiro_p"] <- "p-value (Shapiro)"
    col_names[col_names == "Anderson_Stat"] <- "Statistique A² (Anderson)"
    col_names[col_names == "Anderson_p"] <- "p-value (Anderson)"
    col_names[col_names == "KS_Stat"] <- "Statistique D (KS)"
    col_names[col_names == "KS_p"] <- "p-value (KS)"
    colnames(results) <- col_names
    
    # Créer le tableau avec formatage - pas de boîte de recherche (f retiré du dom)
    dt <- datatable(
      results,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'tip'  # Retiré le 'f' pour supprimer la boîte de recherche
      ),
      caption = "Résultats complets des tests de normalité"
    )
    
    # Formater les statistiques et p-values
    stat_cols <- grep("Statistique", colnames(results))
    if (length(stat_cols) > 0) {
      dt <- formatRound(dt, columns = stat_cols, digits = 4)
    }
    
    p_value_cols <- grep("p-value", colnames(results))
    if (length(p_value_cols) > 0) {
      dt <- formatSignif(dt, columns = p_value_cols, digits = 3)
    }
    
    # Colorer la conclusion
    dt <- formatStyle(
      dt,
      'Conclusion',
      backgroundColor = styleEqual(
        c("Normale", "Non-normale", "Mitigé", "Indéterminé"),
        c("rgba(0, 255, 0, 0.2)", "rgba(255, 0, 0, 0.2)", 
          "rgba(255, 165, 0, 0.2)", "rgba(200, 200, 200, 0.2)")
      )
    )
    
    return(dt)
  })
  

  
  normality_results_data <- reactiveVal(NULL)
  
  # Fonction pour comparer différents nombres de facteurs
  compare_factor_models_reactive <- function() {
    req(sim_data())
    
    # Récupérer les données
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Configurer les paramètres de comparaison
    min_factors <- input$compare_min_factors
    max_factors <- input$compare_max_factors
    
    # Utiliser les valeurs des inputs
    extraction <- input$compare_extraction  
    rotation <- input$compare_rotation      
    corkind <- input$compare_corkind        
    
    # Exécuter la comparaison
    compare_results <- compare_factor_models(
      data = response_data,
      min_factors = min_factors,
      max_factors = max_factors,
      extraction = extraction,
      rotation = rotation,
      corkind = corkind
    )
    
    return(compare_results)
  }
  
  # Observer pour l'analyse factorielle exploratoire
  observeEvent(input$run_efa, {
    req(sim_data())
    
    # Afficher un indicateur de chargement
    withProgress(message = 'Analyse en cours...', {
      # Récupérer les données
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Configurer les paramètres
      extraction <- input$extraction_method
      rotation <- if(input$rotation_type == "none") {
        "none"
      } else {
        input$rotation_method 
      }

      corkind <- input$efa_cor
      
      # Déterminer le nombre de facteurs et le sauvegarder dans une variable
      determined_nfactors <- if(input$factor_method == "manual") {
        input$num_factors
      } else {
        EMPKC(response_data, input$efa_cor, verbose = FALSE)$NfactorsEMPKC
      }
      
      # Exécuter l'analyse avec gestion d'erreur
      tryCatch({
        # Choix entre EFA et PCA
        if(input$analysis_type == "efa") {
          # Effectuer l'analyse factorielle
          result <- EFA.dimensions::EFA(
            data = response_data,
            extraction = extraction,
            corkind = corkind,
            Nfactors = determined_nfactors,
            rotation = rotation,
            verbose = FALSE
          )
          analysisName <- "factorielle exploratoire"
        } else {
          # Effectuer l'analyse en composantes principales
          result <- PCA(
            data = response_data,
            corkind = corkind,
            Nfactors = determined_nfactors,
            rotation = rotation,
            verbose = FALSE
          )
          analysisName <- "en composantes principales"
        }
        
        # Stocker les résultats pour affichage
        result$nfactors <- determined_nfactors
        efa_result(result)
        
        # Notification
        showNotification(paste0("Analyse ", analysisName, " complétée !"), type = "message")
        
      }, error = function(e) {
        # En cas d'erreur
        showNotification(paste("Erreur lors de l'analyse:", e$message), 
                         type = "error", 
                         duration = 10)
      })
    })
  })
  
  # Observateurs pour les boutons Sélectionner/Désélectionner tout
  observeEvent(input$select_all_distrib_vars, {
    req(sim_data())
    data_obj <- sim_data()
    var_names <- colnames(data_obj$response)
    max_choices <- min(20, length(var_names))
    var_choices <- var_names[1:max_choices]
    updateCheckboxGroupInput(session, "selected_distrib_vars", selected = var_choices)
  })
  
  observeEvent(input$deselect_all_distrib_vars, {
    updateCheckboxGroupInput(session, "selected_distrib_vars", selected = character(0))
  })
  
  # Observateurs pour les boutons Sélectionner/Désélectionner tout pour nuages de points
  observeEvent(input$select_all_pairs_vars, {
    req(sim_data())
    data_obj <- sim_data()
    var_names <- colnames(data_obj$response)
    max_choices <- min(15, length(var_names))
    var_choices <- var_names[1:max_choices]
    updateCheckboxGroupInput(session, "selected_pairs_vars", selected = var_choices)
  })
  
  observeEvent(input$deselect_all_pairs_vars, {
    updateCheckboxGroupInput(session, "selected_pairs_vars", selected = character(0))
  })
  
  # Observer pour valider les entrées de vpf personnalisé
  observe({
    if (input$vpf_type == "R") {
      if (input$low_vpf >= input$up_vpf) {
        updateNumericInput(session, "low_vpf", value = input$up_vpf - 1)
      }
    }
  })
  
  # Fonction pour générer les données
  observeEvent(input$generate, {
    # Afficher un indicateur de chargement
    withProgress(message = 'Génération des données...', {
      
      # Préparer les paramètres pour EFAsim.data
      vpf_param <- if(input$vpf_type == "fixed") input$vpf_fixed else input$vpf_type
      fc_param <- if(input$fc_type == "fixed") input$fc_fixed else input$fc_type
      pl_param <- if(input$pl_type == "fixed") input$pl_fixed else input$pl_type
      cl_param <- if(input$cl_type == "fixed") input$cl_fixed else input$cl_type
      
      # Paramètres beta
      a_param <- NULL
      b_param <- NULL
      if(input$distri == "beta" && input$beta_params == "custom") {
        a_param <- input$a_param
        b_param <- input$b_param
      }
      
      # Générer les données
      tryCatch({
        data_obj <- EFAsim.data(
          nfact = input$nfact, 
          vpf = vpf_param, 
          N = input$N, 
          distri = input$distri, 
          fc = fc_param, 
          pl = pl_param, 
          cl = cl_param,
          low.vpf = input$low_vpf, 
          up.vpf = input$up_vpf, 
          a = a_param, 
          b = b_param, 
          vis = FALSE
        )
        
        # Stocker les données dans la réactive
        sim_data(data_obj)
        
        # Notification de succès
        showNotification("Données générées avec succès!", type = "message")
        
      }, error = function(e) {
        # En cas d'erreur
        showNotification(paste("Erreur lors de la génération des données:", e$message), 
                         type = "error", 
                         duration = 10)
      })
    })
  })
  
  # Préréglage pour la structure de personnalité
  observeEvent(input$preset_personality, {
    # Mettre à jour les paramètres pour simuler des données de personnalité
    updateNumericInput(session, "nfact", value = 5)  # Les 5 grands facteurs
    updateSelectInput(session, "vpf_type", selected = "fixed")
    updateNumericInput(session, "vpf_fixed", value = 5)  # 5 items par facteur
    updateNumericInput(session, "N", value = 500)
    updateSelectInput(session, "distri", selected = "normal")
    updateSelectInput(session, "fc_type", selected = "fixed")
    updateSliderInput(session, "fc_fixed", value = 0.3)  # Facteurs modérément corrélés
    updateSelectInput(session, "pl_type", selected = "M")  # Charges moyennes
    updateSelectInput(session, "cl_type", selected = "L")  # Faibles charges secondaires
    
    # Notification confirmant le préréglage
    showNotification(
      HTML("<b>Préréglage Personnalité appliqué</b><br>5 facteurs, 5 items par facteur<br>Simuler des données semblables aux inventaires de personnalité"),
      type = "message", 
      duration = 5
    )
  })
  
  # Préréglage pour les tests cognitifs
  observeEvent(input$preset_cognitive, {
    # Mettre à jour les paramètres pour simuler des données de tests cognitifs
    updateNumericInput(session, "nfact", value = 3)  # Les habilités cognitives sont souvent 2-3 facteurs
    updateSelectInput(session, "vpf_type", selected = "fixed")
    updateNumericInput(session, "vpf_fixed", value = 6)  # Généralement plus d'items par facteur
    updateNumericInput(session, "N", value = 400)
    updateSelectInput(session, "distri", selected = "normal")
    updateSelectInput(session, "fc_type", selected = "fixed")
    updateSliderInput(session, "fc_fixed", value = 0.4)  # Facteurs plus corrélés pour cognitif
    updateSelectInput(session, "pl_type", selected = "H")  # Charges plus élevées
    updateSelectInput(session, "cl_type", selected = "None")  # Typiquement pas de charges croisées
    
    # Notification confirmant le préréglage
    showNotification(
      HTML("<b>Préréglage Tests Cognitifs appliqué</b><br>3 facteurs, 6 items par facteur<br>Simuler des données semblables aux tests d'aptitudes cognitives"),
      type = "message", 
      duration = 5
    )
  })
  
  # Préréglage pour les échelles d'attitudes
  observeEvent(input$preset_social, {
    # Mettre à jour les paramètres pour simuler des données d'échelles d'attitudes
    updateNumericInput(session, "nfact", value = 4)  # Souvent 3-5 dimensions d'attitudes
    updateSelectInput(session, "vpf_type", selected = "M")  # Nombre variable d'items
    updateNumericInput(session, "N", value = 600)  # Échantillon plus grand
    updateSelectInput(session, "distri", selected = "beta")  # Distribution souvent asymétrique 
    updateSelectInput(session, "beta_params", selected = "custom")
    updateNumericInput(session, "a_param", value = 2)
    updateNumericInput(session, "b_param", value = 3)  # Légère asymétrie
    updateSelectInput(session, "fc_type", selected = "fixed")
    updateSliderInput(session, "fc_fixed", value = 0.25)  # Corrélations modérées
    updateSelectInput(session, "pl_type", selected = "M")  # Charges moyennes
    updateSelectInput(session, "cl_type", selected = "H")  # Plus de charges croisées
    
    # Notification confirmant le préréglage
    showNotification(
      HTML("<b>Préréglage Échelles d'Attitudes appliqué</b><br>4 facteurs, distribution Beta<br>Simuler des données semblables aux échelles d'attitudes sociales"),
      type = "message", 
      duration = 5
    )
  })
  
  # Préréglage pour la performance académique
  observeEvent(input$preset_academic, {
    # Mettre à jour les paramètres pour simuler des données de performance
    updateNumericInput(session, "nfact", value = 2)  # Souvent 2 facteurs (compétence/performance)
    updateSelectInput(session, "vpf_type", selected = "fixed")
    updateNumericInput(session, "vpf_fixed", value = 8)  # Plus d'items par facteur
    updateNumericInput(session, "N", value = 300)
    updateSelectInput(session, "distri", selected = "normal")
    updateSelectInput(session, "fc_type", selected = "fixed")
    updateSliderInput(session, "fc_fixed", value = 0.5)  # Corrélations plus élevées
    updateSelectInput(session, "pl_type", selected = "H")  # Charges élevées
    updateSelectInput(session, "cl_type", selected = "L")  # Charges croisées faibles
    
    # Notification confirmant le préréglage
    showNotification(
      HTML("<b>Préréglage Performance Académique appliqué</b><br>2 facteurs, 8 items par facteur<br>Simuler des données semblables aux évaluations de performance académique"),
      type = "message", 
      duration = 5
    )
  })
  
  # Observer pour la comparaison de modèles factoriels
  observeEvent(input$run_compare, {
    # Afficher un indicateur de chargement
    withProgress(message = 'Comparaison des modèles factoriels...', {
      # Exécuter la comparaison avec gestion d'erreur
      tryCatch({
        results <- compare_factor_models_reactive()
        compare_results(results)
        showNotification("Comparaison des modèles complétée !", type = "message")
      }, error = function(e) {
        showNotification(paste("Erreur lors de la comparaison:", e$message), 
                         type = "error", duration = 10)
      })
    })
  })
  
  # Tableau comparatif
  output$compare_table <- renderDT({
    req(compare_results())
    
    results <- compare_results()
    comparison_df <- results$comparison_df
    
    # Formater les noms des colonnes pour l'affichage
    col_names <- names(comparison_df)
    display_names <- col_names
    display_names[col_names == "Factors"] <- "Facteurs"
    display_names[col_names == "Variance_Explained"] <- "Variance expliquée (%)"
    
    # Mettre en forme le tableau avec DT
    datatable(
      comparison_df,
      colnames = display_names,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 't'
      ),
      caption = "Comparaison des indices d'ajustement par nombre de facteurs"
    ) %>%
      formatRound(columns = 2:ncol(comparison_df), digits = 3)
  })
  
  # Interprétation des résultats de comparaison
  output$compare_interpretation <- renderPrint({
    req(compare_results())
    
    results <- compare_results()
    comparison_df <- results$comparison_df
    
    cat("INTERPRÉTATION DES RÉSULTATS DE COMPARAISON DE MODÈLES\n")
    cat("===================================================\n\n")
    
    cat("1. INDICES D'AJUSTEMENT PAR NOMBRE DE FACTEURS\n")
    cat("------------------------------------------\n")
    
    # Afficher les indices d'ajustement dans un format lisible
    print(comparison_df)
    
    cat("\n\n2. ANALYSE DES CRITÈRES D'INFORMATION\n")
    cat("---------------------------------\n")
    
    # Analyser les critères d'information (AIC, BIC, etc.)
    info_criteria <- c("AIC", "BIC", "CAIC", "SABIC")
    available_criteria <- intersect(info_criteria, names(comparison_df))
    
    if (length(available_criteria) > 0) {
      cat("Critères d'information disponibles :\n")
      for (criterion in available_criteria) {
        best_model <- which.min(comparison_df[[criterion]])
        best_n_factors <- comparison_df$Factors[best_model]
        cat(sprintf("- %s : Favorise le modèle à %d facteurs (valeur = %.2f)\n", 
                    criterion, best_n_factors, comparison_df[[criterion]][best_model]))
      }
    } else {
      cat("Aucun critère d'information (AIC, BIC, etc.) n'est disponible dans les résultats.\n")
    }
    
    cat("\n3. ANALYSE DES INDICES D'AJUSTEMENT\n")
    cat("--------------------------------\n")
    
    # Analyser les indices d'ajustement (CFI, TLI, RMSEA, etc.)
    fit_indices <- c("CFI", "TLI", "RMSEA", "RMSR", "GFI", "CAF")
    available_indices <- intersect(fit_indices, names(comparison_df))
    
    if (length(available_indices) > 0) {
      cat("Indices d'ajustement disponibles :\n")
      for (index in available_indices) {
        if (index %in% c("RMSEA", "RMSR")) {
          # Pour ces indices, une valeur plus basse est meilleure
          best_model <- which.min(comparison_df[[index]])
          best_value <- comparison_df[[index]][best_model]
          
          # Interprétation pour RMSEA
          if (index == "RMSEA") {
            quality <- ifelse(best_value <= 0.05, "Excellent",
                              ifelse(best_value <= 0.08, "Bon",
                                     ifelse(best_value <= 0.10, "Acceptable", "Médiocre")))
          } 
          # Interprétation pour RMSR
          else if (index == "RMSR") {
            quality <- ifelse(best_value <= 0.05, "Excellent",
                              ifelse(best_value <= 0.08, "Bon", "Médiocre"))
          }
        } else {
          # Pour ces indices, une valeur plus élevée est meilleure
          best_model <- which.max(comparison_df[[index]])
          best_value <- comparison_df[[index]][best_model]
          
          # Interprétation pour CFI et TLI
          if (index %in% c("CFI", "TLI")) {
            quality <- ifelse(best_value >= 0.95, "Excellent",
                              ifelse(best_value >= 0.90, "Acceptable", "Médiocre"))
          } 
          # Interprétation pour GFI
          else if (index == "GFI") {
            quality <- ifelse(best_value >= 0.95, "Excellent",
                              ifelse(best_value >= 0.90, "Bon", "Médiocre"))
          }
          # Interprétation par défaut
          else {
            quality <- "Non disponible"
          }
        }
        
        best_n_factors <- comparison_df$Factors[best_model]
        cat(sprintf("- %s : Favorise le modèle à %d facteurs (valeur = %.2f, qualité = %s)\n", 
                    index, best_n_factors, best_value, quality))
      }
    } else {
      cat("Aucun indice d'ajustement standard (CFI, TLI, RMSEA, etc.) n'est disponible dans les résultats.\n")
    }
    
    cat("\n4. ANALYSE DE LA VARIANCE EXPLIQUÉE\n")
    cat("------------------------------\n")
    
    # Analyser la variance expliquée
    if ("Variance_Explained" %in% names(comparison_df)) {
      cat("Variance expliquée par nombre de facteurs :\n")
      for (i in 1:nrow(comparison_df)) {
        n_factors <- comparison_df$Factors[i]
        var_explained <- comparison_df$Variance_Explained[i]
        
        # Calculer l'augmentation de variance par rapport au modèle précédent
        increase <- if (i > 1) {
          var_explained - comparison_df$Variance_Explained[i-1]
        } else {
          var_explained
        }
        
        cat(sprintf("- Modèle à %d facteur(s) : %.2f%% de variance expliquée", 
                    n_factors, var_explained))
        
        if (i > 1) {
          cat(sprintf(" (augmentation de %.2f%%)\n", increase))
        } else {
          cat("\n")
        }
      }
      
      # Règle du "coude" sur la variance expliquée
      cat("\nAnalyse du 'coude' sur la variance expliquée :\n")
      increases <- numeric(nrow(comparison_df)-1)
      for (i in 2:nrow(comparison_df)) {
        increases[i-1] <- comparison_df$Variance_Explained[i] - comparison_df$Variance_Explained[i-1]
      }
      
      # Trouver où se produit la diminution la plus importante dans les augmentations
      if (length(increases) >= 2) {
        diff_increases <- diff(increases)
        elbow_idx <- which.min(diff_increases) + 1
        
        if (length(elbow_idx) > 0 && elbow_idx <= length(increases)) {
          elbow_n_factors <- comparison_df$Factors[elbow_idx]
          cat(sprintf("Le 'coude' dans la variance expliquée suggère %d facteurs.\n", elbow_n_factors))
          cat(sprintf("(L'ajout de facteurs au-delà de %d apporte des gains de variance moins importants)\n", elbow_n_factors))
        }
      }
    } else {
      cat("Les données sur la variance expliquée ne sont pas disponibles.\n")
    }
    
    cat("\n5. SYNTHÈSE ET RECOMMANDATION\n")
    cat("---------------------------\n")
    
    # Compiler les recommandations de tous les critères
    recommended_factors <- c()
    
    # Ajouter les recommandations basées sur les critères d'information
    for (criterion in available_criteria) {
      best_model <- which.min(comparison_df[[criterion]])
      recommended_factors <- c(recommended_factors, comparison_df$Factors[best_model])
    }
    
    # Ajouter les recommandations basées sur les indices d'ajustement
    for (index in available_indices) {
      if (index %in% c("RMSEA", "RMSR")) {
        best_model <- which.min(comparison_df[[index]])
      } else {
        best_model <- which.max(comparison_df[[index]])
      }
      recommended_factors <- c(recommended_factors, comparison_df$Factors[best_model])
    }
    
    # Ajouter la recommandation basée sur le "coude" de variance expliquée
    if (exists("elbow_n_factors")) {
      recommended_factors <- c(recommended_factors, elbow_n_factors)
    }
    
    # Calculer le mode (valeur la plus fréquente)
    if (length(recommended_factors) > 0) {
      recommendation_table <- table(recommended_factors)
      mode_factors <- as.numeric(names(recommendation_table)[which.max(recommendation_table)])
      mode_count <- max(recommendation_table)
      mode_percent <- round(mode_count / length(recommended_factors) * 100)
      
      cat(sprintf("Sur la base de l'analyse des critères ci-dessus, la solution à %d facteurs est recommandée.\n", mode_factors))
      cat(sprintf("Cette solution est favorisée par %d critères sur %d (%d%%).\n", 
                  mode_count, length(recommended_factors), mode_percent))
      
      # Ajouter un commentaire sur la cohérence des recommandations
      if (mode_percent >= 75) {
        cat("Les critères montrent une forte convergence vers cette solution.\n")
      } else if (mode_percent >= 50) {
        cat("Les critères montrent une convergence modérée vers cette solution.\n")
      } else {
        cat("Les critères ne montrent pas une forte convergence. Considérez d'examiner plusieurs solutions factorielles.\n")
      }
      
      # Recommandation finale
      cat("\nRECOMMANDATION FINALE :\n")
      cat(sprintf("Examinez en détail la solution à %d facteurs.\n", mode_factors))
      
      # Si plusieurs solutions sont proches, les mentionner également
      runner_up <- sort(recommendation_table, decreasing = TRUE)
      if (length(runner_up) > 1 && runner_up[2] >= runner_up[1] * 0.7) {
        runner_up_factors <- as.numeric(names(runner_up)[2])
        cat(sprintf("Considérez également la solution à %d facteurs comme alternative viable.\n", runner_up_factors))
      }
    } else {
      cat("Insuffisance de critères pour formuler une recommandation solide.\n")
      cat("Basez votre décision sur des considérations théoriques ou d'autres analyses.\n")
    }
  })
  
  # Afficher la matrice de corrélation des facteurs
  output$factor_corr_plot <- renderPlot({
    req(sim_data())
    data_obj <- sim_data()
    
    corrplot(data_obj$cor.factors, method = "color", 
             type = "upper", 
             addCoef.col = "black",
             tl.col = "black", tl.srt = 45,
             diag = FALSE,
             title = "Corrélation entre facteurs",
             mar = c(0, 0, 2, 0))
  })
  
  # UI pour la matrice de corrélation avec méthode de clustering corrigée
  output$item_corr_plot <- renderPlot({
    req(sim_data(), input$selected_corr_vars)
    
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Vérifier qu'au moins 2 variables sont sélectionnées
    if (length(input$selected_corr_vars) < 2) {
      plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
      text(0, 0, "Sélectionnez au moins 2 variables pour afficher la matrice de corrélation", cex = 1.5)
      return()
    }
    
    # Sélectionner les variables à inclure
    var_subset <- response_data[, input$selected_corr_vars, drop = FALSE]
    
    # Calculer la matrice de corrélation
    cor_matrix <- cor(var_subset, method = input$corr_method)
    
    # Déterminer le titre
    n_vars_total <- ncol(response_data)
    n_vars_shown <- ncol(var_subset)
    title <- if (n_vars_shown < n_vars_total) {
      paste0("Matrice de corrélation (", n_vars_shown, " variables sur ", n_vars_total, ")")
    } else {
      "Matrice de corrélation complète"
    }
    
    # Important: permettre le rafraîchissement quand la méthode de clustering change
    if (input$corr_plot_order == "hclust") {
      req(input$corr_clust_method)
    }
    
    # Solution 1: Essayons avec un bloc tryCatch pour gérer l'erreur
    tryCatch({
      # Visualiser avec corrplot et les options améliorées
      corrplot(
        cor_matrix, 
        method = input$corr_plot_method,
        type = input$corr_plot_type,
        order = input$corr_plot_order,
        addCoef.col = if(input$show_corr_values) "black" else NULL,
        number.cex = 0.7,
        tl.col = "black",
        tl.srt = 45,
        tl.cex = 0.8,
        diag = input$corr_plot_diag,
        title = title,
        mar = c(0, 0, 0.5, 0),
        cl.ratio = 0.2,
        col = COL2("RdBu", 200),
        hclust.method = if(input$corr_plot_order == "hclust") input$corr_clust_method else "complete"
      )
    }, error = function(e) {
      # En cas d'erreur, essayer une version simplifiée
      corrplot(
        cor_matrix, 
        method = input$corr_plot_method,
        type = input$corr_plot_type,
        order = "original",  # Défaut à original en cas d'erreur avec hclust
        addCoef.col = if(input$show_corr_values) "black" else NULL,
        number.cex = 0.7,
        tl.col = "black",
        tl.srt = 45,
        tl.cex = 0.8,
        diag = input$corr_plot_diag,
        title = paste0(title, "\n(Mode d'ordonnancement simplifié en raison d'une erreur)"),
        mar = c(0, 0, 1, 0),
        cl.ratio = 0.2,
        col = COL2("RdBu", 200)
      )
    })
  })
  
  # Téléchargement de la matrice de corrélation - avec la même correction
  output$download_corr_plot <- downloadHandler(
    filename = function() {
      paste0("correlation_matrix_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Créer un fichier PNG plus grand pour une meilleure qualité
      png(file, width = 1200, height = 1000, res = 120)
      
      # Reproduire le même graphique
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Vérifier qu'au moins 2 variables sont sélectionnées
      if (length(input$selected_corr_vars) < 2) {
        plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
        text(0, 0, "Sélectionnez au moins 2 variables pour afficher la matrice de corrélation", cex = 1.5)
      } else {
        # Sélectionner les variables à inclure
        var_subset <- response_data[, input$selected_corr_vars, drop = FALSE]
        
        # Calculer la matrice de corrélation
        cor_matrix <- cor(var_subset, method = input$corr_method)
        
        # Déterminer le titre
        n_vars_total <- ncol(response_data)
        n_vars_shown <- ncol(var_subset)
        title <- if (n_vars_shown < n_vars_total) {
          paste0("Matrice de corrélation (", n_vars_shown, " variables sur ", n_vars_total, ")")
        } else {
          "Matrice de corrélation complète"
        }
        
        # Solution équivalente pour le téléchargement
        tryCatch({
          corrplot(
            cor_matrix, 
            method = input$corr_plot_method,
            type = input$corr_plot_type,
            order = input$corr_plot_order,
            addCoef.col = if(input$show_corr_values) "black" else NULL,
            number.cex = 0.7,
            tl.col = "black",
            tl.srt = 45,
            tl.cex = 0.8,
            diag = input$corr_plot_diag,
            title = title,
            mar = c(0, 0, 0.5, 0),
            cl.ratio = 0.2,
            col = COL2("RdBu", 200),
            hclust.method = if(input$corr_plot_order == "hclust") input$corr_clust_method else "complete"
          )
        }, error = function(e) {
          corrplot(
            cor_matrix, 
            method = input$corr_plot_method,
            type = input$corr_plot_type,
            order = "original",
            addCoef.col = if(input$show_corr_values) "black" else NULL,
            number.cex = 0.7,
            tl.col = "black",
            tl.srt = 45,
            tl.cex = 0.8,
            diag = input$corr_plot_diag,
            title = paste0(title, "\n(Mode d'ordonnancement simplifié)"),
            mar = c(0, 0, 1, 0),
            cl.ratio = 0.2,
            col = COL2("RdBu", 200)
          )
        })
      }
      
      dev.off()
    }
  )
  

  
  # UI pour sélectionner manuellement les variables pour le nuage de points
  output$pairs_var_selector <- renderUI({
    req(sim_data())
    data_obj <- sim_data()
    var_names <- colnames(data_obj$response)
    
    # Limiter à 15 variables maximum pour éviter une interface trop chargée
    max_choices <- min(15, length(var_names))
    var_choices <- var_names[1:max_choices]
    
    # Sélectionner par défaut les 5 premières variables ou moins
    default_selected <- var_choices[1:min(5, length(var_choices))]
    
    # Déterminer le nombre de colonnes en fonction du nombre de variables
    n_columns <- if(length(var_choices) > 15) 3 else if(length(var_choices) > 8) 2 else 1
    
    tagList(
      checkboxGroupInput(
        "selected_pairs_vars",
        NULL, # Pas de label car déjà indiqué au-dessus
        choices = var_choices,
        selected = default_selected,
        inline = FALSE
      ) %>% tagAppendAttributes(class = paste0("multicol-", n_columns)),  # Ajouter une classe pour le style multi-colonnes
      
      # Ajouter du CSS pour le style multi-colonnes
      tags$style(HTML(paste0("
      .multicol-2 .shiny-options-group {
        -webkit-column-count: 2;
        -moz-column-count: 2;
        column-count: 2;
        -webkit-column-gap: 10px;
        -moz-column-gap: 10px;
        column-gap: 10px;
      }
      .multicol-3 .shiny-options-group {
        -webkit-column-count: 3;
        -moz-column-count: 3;
        column-count: 3;
        -webkit-column-gap: 10px;
        -moz-column-gap: 10px;
        column-gap: 10px;
      }
    "))),
      
      div(style = "display: flex; justify-content: space-between; margin-top: 10px;",
          actionLink("select_all_pairs_vars", "Tout sélectionner"),
          actionLink("deselect_all_pairs_vars", "Tout désélectionner")
      )
    )
  })
  
  
  
  # Matrice de nuages de points
  output$pairs_plot <- renderPlot({
    req(sim_data(), input$selected_pairs_vars)
    
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Vérifier qu'au moins 2 variables sont sélectionnées
    if (length(input$selected_pairs_vars) < 2) {
      plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
      text(0, 0, "Sélectionnez au moins 2 variables pour afficher les nuages de points", cex = 1.5)
      return()
    }
    
    # Sélectionner les variables
    selected_vars <- input$selected_pairs_vars
    
    # Créer la matrice de nuages de points
    pairs_data <- response_data[, selected_vars, drop = FALSE]
    
    # Ajuster les options du graphique
    pairs(
      pairs_data,
      pch = 19,
      col = adjustcolor(input$pairs_color_scheme, alpha.f = input$pairs_alpha),
      cex = 0.7,
      main = paste0("Matrice de nuages de points (", length(selected_vars), " variables)"),
      gap = 0.5,
      lower.panel = function(x, y, ...) {
        points(x, y, pch = 19, col = adjustcolor(input$pairs_color_scheme, alpha.f = input$pairs_alpha), cex = 0.7)
        if (input$show_ellipses) {
          # Vérifier que car::dataEllipse est disponible
          if (requireNamespace("car", quietly = TRUE)) {
            car::dataEllipse(x, y, add = TRUE, plot.points = FALSE, 
                             levels = 0.95, col = "red", lwd = 1.5)
          }
        }
        if (input$show_smooth) {
          # Ajouter la ligne de régression avec gestion d'erreur
          tryCatch({
            abline(lm(y ~ x), col = "darkgreen", lwd = 1.5)
          }, error = function(e) {
            # Ne rien faire si la régression échoue
          })
        }
      },
      upper.panel = function(x, y, ...) {
        # Calculer la corrélation avec gestion d'erreur
        r <- tryCatch({
          cor(x, y, method = input$corr_method, use = "pairwise.complete.obs")
        }, error = function(e) {
          NA
        })
        
        # Vérifier que r n'est pas NA avant de l'utiliser dans les conditions
        if (!is.na(r)) {
          # Colorier selon la force de la corrélation
          text_col <- if (abs(r) > 0.7) {
            "darkred"
          } else if (abs(r) > 0.4) {
            "darkorange"
          } else {
            "darkblue"
          }
          text(mean(range(x)), mean(range(y)), sprintf("r = %.2f", r), 
               cex = 1.25 + abs(r), col = text_col, font = 2)
        } else {
          # Afficher un message d'erreur si la corrélation n'a pas pu être calculée
          text(mean(range(x)), mean(range(y)), "r = NA", 
               cex = 1, col = "gray", font = 2)
        }
      }
    )
  })

  # Téléchargement de la matrice de nuages de points
  output$download_pairs_plot <- downloadHandler(
    filename = function() {
      paste0("scatterplot_matrix_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      # Créer un fichier PNG plus grand pour une meilleure qualité
      png(file, width = 1200, height = 1000, res = 120)
      
      # Reproduire le même graphique
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Vérifier qu'au moins 2 variables sont sélectionnées
      if (length(input$selected_pairs_vars) < 2) {
        plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
        text(0, 0, "Sélectionnez au moins 2 variables pour afficher les nuages de points", cex = 1.5)
      } else {
        # Sélectionner les variables
        selected_vars <- input$selected_pairs_vars
        
        # Créer la matrice de nuages de points
        pairs_data <- response_data[, selected_vars, drop = FALSE]
        
        # Ajuster les options du graphique
        pairs(
          pairs_data,
          pch = 19,
          col = adjustcolor(input$pairs_color_scheme, alpha.f = input$pairs_alpha),
          cex = 0.7,
          main = paste0("Matrice de nuages de points (", length(selected_vars), " variables)"),
          gap = 0.5,
          lower.panel = function(x, y, ...) {
            points(x, y, pch = 19, col = adjustcolor(input$pairs_color_scheme, alpha.f = input$pairs_alpha), cex = 0.7)
            if (input$show_ellipses) {
              # Vérifier que car::dataEllipse est disponible
              if (requireNamespace("car", quietly = TRUE)) {
                car::dataEllipse(x, y, add = TRUE, plot.points = FALSE, 
                                 levels = 0.95, col = "red", lwd = 1.5)
              }
            }
            if (input$show_smooth) {
              # Ajouter la ligne de régression avec gestion d'erreur
              tryCatch({
                abline(lm(y ~ x), col = "darkgreen", lwd = 1.5)
              }, error = function(e) {
                # Ne rien faire si la régression échoue
              })
            }
          },
          upper.panel = function(x, y, ...) {
            # Calculer la corrélation avec gestion d'erreur
            r <- tryCatch({
              cor(x, y, method = input$corr_method, use = "pairwise.complete.obs")
            }, error = function(e) {
              NA
            })
            
            # Vérifier que r n'est pas NA avant de l'utiliser
            if (!is.na(r)) {
              # Colorier selon la force de la corrélation
              text_col <- if (abs(r) > 0.7) {
                "darkred"
              } else if (abs(r) > 0.4) {
                "darkorange"
              } else {
                "darkblue"
              }
              text(mean(range(x)), mean(range(y)), sprintf("r = %.2f", r), 
                   cex = 1.25 + abs(r), col = text_col, font = 2)
            } else {
              # Afficher un message si la corrélation n'a pas pu être calculée
              text(mean(range(x)), mean(range(y)), "r = NA", 
                   cex = 1, col = "gray", font = 2)
            }
          }
        )
      }
      
      dev.off()
    }
  )
  
  # Afficher les données simulées
  output$simulated_data <- renderDT({
    req(sim_data())
    data_obj <- sim_data()
    
    # Limiter à un échantillon plus petit pour l'affichage
    response_data <- data_obj$response
    if (nrow(response_data) > 100) {
      display_data <- response_data[1:100, ]
    } else {
      display_data <- response_data
    }
    
    datatable(
      display_data,
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        dom = 'rtip'
      ),
      caption = ifelse(nrow(response_data) > 100, 
                       "Affichage des 100 premières observations (sur {nrow(response_data)})",
                       "Données de réponse simulées complètes")
    ) %>%
      formatRound(columns = 1:ncol(display_data), digits = 3)
  })
  
  # Afficher la structure des facteurs
  output$factor_structure <- renderPrint({
    req(sim_data())
    data_obj <- sim_data()
    
    cat("STRUCTURE DES FACTEURS ET ITEMS\n")
    cat("==============================\n\n")
    
    for (f in 1:length(data_obj$items)) {
      cat("Facteur", f, ":\n")
      cat("  Variables:", paste(names(data_obj$items[[f]]), collapse=", "), "\n\n")
    }
    
    cat("\nINFORMATIONS GÉNÉRALES\n")
    cat("=====================\n")
    cat("Nombre total de facteurs:", length(data_obj$items), "\n")
    cat("Nombre total de variables:", ncol(data_obj$response), "\n")
    cat("Nombre d'observations:", nrow(data_obj$response), "\n")
  })
  
  # Visualisation des données simulées
  output$visualization_plot <- renderPlot({
    req(sim_data())
    data_obj <- sim_data()
    
    response_data <- as.data.frame(data_obj$response)
    
    if (input$visualization_type == "hist") {
      # Histogrammes des variables
      # Limiter à 16 variables pour la lisibilité
      n_vars <- min(16, ncol(response_data))
      vars_to_plot <- names(response_data)[1:n_vars]
      
      # Définir la disposition de la grille
      par(mfrow = c(4, 4), mar = c(4, 4, 2, 1))
      
      for (var in vars_to_plot) {
        hist(response_data[, var], 
             main = var,
             xlab = "",
             col = "skyblue",
             border = "white",
             probability = TRUE)
        lines(density(response_data[, var]), col = "red", lwd = 2)
      }
      
    } else if (input$visualization_type == "heat") {
      # Heatmap des corrélations
      # Limiter à 30 variables pour la lisibilité
      n_vars <- min(30, ncol(response_data))
      cors <- cor(response_data[, 1:n_vars])
      
      corrplot(cors, method = "color", 
               type = "upper", 
               tl.col = "black", tl.srt = 45,
               tl.cex = 0.7,
               diag = FALSE,
               title = "Heatmap des corrélations entre variables",
               mar = c(0, 0, 2, 0))
      
    } else if (input$visualization_type == "pairs") {
      # Scatterplot matrix (limité à 6 variables)
      n_vars <- min(6, ncol(response_data))
      vars_to_plot <- names(response_data)[1:n_vars]
      
      pairs(response_data[, vars_to_plot], 
            pch = 19, 
            col = rgb(0, 0, 1, 0.3),
            main = "Matrice de nuages de points (6 premières variables)")
    }
  })

  output$structure_title <- renderUI({
    req(sim_data())
    HTML(paste0("<h4>Structure à ", attr(sim_data(), "n.factors"), " facteurs avec ", 
                nrow(attr(sim_data(), "loadings")), " variables</h4>"))
  })
  
  observeEvent(input$run_dimtests, {
    req(sim_data())
    
    # Vérifier que des tests sont sélectionnés
    if (length(input$dimtests) == 0) {
      showNotification("Veuillez sélectionner au moins un test.", type = "warning")
      return()
    }
    
    # Afficher un indicateur de chargement
    withProgress(message = 'Exécution des tests dimensionnels...', {
      # Récupérer les données
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Configurer les paramètres
      hull_method <- input$dimtests_hull_method
      hull_gof <- input$dimtests_hull_gof
      corkind <- input$dimtests_cor
      
      # Créer une structure pour stocker les résultats
      results <- list()
      
      # Stocker les valeurs propres avec gestion d'erreur
      tryCatch({
        results$eigenvalues <- svd(cor(response_data, method = corkind))$d
      }, error = function(e) {
        showNotification(paste("Erreur lors du calcul des valeurs propres:", e$message), 
                         type = "error", duration = 10)
        results$eigenvalues <- NULL
      })
      
      # Initialiser les résultats finaux
      final_results <- list()
      
      # Traiter le test NEVALSGT1 manuellement
      if ("NEVALSGT1" %in% input$dimtests && !is.null(results$eigenvalues)) {
        # Calculer directement
        n_factors_gt1 <- sum(results$eigenvalues > 1)
        final_results[["NEVALSGT1"]] <- n_factors_gt1
      }
      
      # Utiliser RAWPAR si sélectionné
      if ("PARALLEL" %in% input$dimtests) {
        tryCatch({
          # Récupérer les paramètres de l'interface utilisateur
          randtype <- input$rawpar_randtype     # "generated" ou "permuted"
          extraction <- input$rawpar_extraction # "PAF", "PCA" ou "image"
          ndatasets <- input$rawpar_ndatasets   # nombre de jeux de données aléatoires
          percentile <- input$rawpar_percentile # percentile à utiliser 
          corkindRAND <- input$rawpar_corkindRAND # type de corrélation pour données aléatoires
          
          # Exécuter RAWPAR
          rawpar_result <- EFA.dimensions::RAWPAR(
            data = response_data,
            randtype = randtype,
            extraction = extraction,
            Ndatasets = ndatasets,
            percentile = percentile,
            corkind = corkind,
            corkindRAND = corkindRAND
          )
          
          # Vérifier que rawpar_result est une liste avant d'accéder à ses éléments
          if (is.list(rawpar_result)) {
            # Stocker les résultats du test parallèle
            if (!is.null(rawpar_result$eigenvalues)) {
              results$parallel_data <- list(
                observed = rawpar_result$eigenvalues[, "Real Data"],
                Mean = rawpar_result$eigenvalues[, "Mean"],
                Percentile = rawpar_result$eigenvalues[, "Percentile"]
              )
            }
            
            # Stocker le nombre de facteurs trouvé
            if (!is.null(rawpar_result$NfactorsPA)) {
              final_results[["PARALLEL"]] <- rawpar_result$NfactorsPA
            } else {
              showNotification("Le test parallèle n'a pas pu déterminer le nombre de facteurs.", type = "warning", duration = 5)
            }
          } else {
            warning("Format de résultat du test parallèle inattendu")
          }
          
        }, error = function(e) {
          showNotification(paste("Erreur lors de l'exécution du test parallèle:", e$message), 
                           type = "error", duration = 10)
        })
      }
      
      # Pour les autres tests, capturer la sortie de DIMTESTS
      other_tests <- setdiff(input$dimtests, c("NEVALSGT1", "PARALLEL"))
      
      if (length(other_tests) > 0) {
        # Exécution de DIMTESTS avec gestion d'erreur
        tryCatch({
          # Capturer la sortie de la console
          output_lines <- capture.output({
            dim_tests_result <- EFA.dimensions::DIMTESTS(
              data = response_data,
              tests = other_tests,
              corkind = corkind,
              HULL_method = hull_method,
              HULL_gof = hull_gof,
              HULL_cor_method = corkind,
              CD_cor_method = corkind,
              display = 2  # Display=2 pour voir dans la console
            )
          })
          
          # Stocker la sortie console pour l'affichage détaillé
          results$console_output <- paste(output_lines, collapse = "\n")
          
          # Approche générique pour tous les tests
          for (test in other_tests) {
            # Pattern pour trouver les résultats dans le format "TEST              N"
            pattern <- paste0("^", test, "[[:space:]]+([0-9]+)$")
            test_line <- grep(pattern, output_lines, value = TRUE)
            
            if (length(test_line) > 0) {
              test_value <- suppressWarnings(as.numeric(gsub(paste0("^", test, "[[:space:]]+([0-9]+)$"), "\\1", test_line)))
              if (!is.na(test_value)) {
                final_results[[test]] <- test_value
              }
            }
          }
          
          # Extraire également les informations détaillées pour chaque test
          # CD (Coordonnées optimales)
          if ("CD" %in% other_tests) {
            cd_lines <- grep("CD method", output_lines, value = TRUE, fixed = TRUE)
            if (length(cd_lines) > 0) {
              results$CD_details <- cd_lines
            }
          }
          
          # EMPKC (Critère d'information empirique de Kaiser)
          if ("EMPKC" %in% other_tests) {
            empkc_lines <- grep("EMPKC", output_lines, value = TRUE, fixed = TRUE)
            if (length(empkc_lines) > 0) {
              results$EMPKC_details <- empkc_lines
            }
          }
          
          # HULL
          if ("HULL" %in% other_tests) {
            hull_start <- grep("HULL method", output_lines, fixed = TRUE)
            if (length(hull_start) > 0) {
              hull_end <- min(grep("^HULL[[:space:]]+[0-9]+$", output_lines), length(output_lines))
              if (length(hull_end) > 0 && hull_end > hull_start) {
                results$HULL_details <- output_lines[hull_start:hull_end]
              }
            }
          }
          
          # MAP (Minimum Average Partial)
          if ("MAP" %in% other_tests) {
            map_start <- grep("MAP method", output_lines, fixed = TRUE)
            if (length(map_start) > 0) {
              map_end <- min(grep("^MAP[[:space:]]+[0-9]+$", output_lines), length(output_lines))
              if (length(map_end) > 0 && map_end > map_start) {
                results$MAP_details <- output_lines[map_start:map_end]
              }
            }
          }
          
          # SESCREE (Scree test)
          if ("SESCREE" %in% other_tests) {
            sescree_lines <- grep("SESCREE", output_lines, value = TRUE, fixed = TRUE)
            if (length(sescree_lines) > 0) {
              results$SESCREE_details <- sescree_lines
            }
          }
        }, error = function(e) {
          showNotification(paste("Erreur lors de l'exécution des tests dimensionnels:", e$message), 
                           type = "error", duration = 10)
        })
      }
      
      # Pour le test SMT (Test de dimensionnalité séquentielle)
      if ("SMT" %in% input$dimtests) {
        tryCatch({
          # Utiliser directement SMT au lieu de passer par DIMTESTS
          smt_result <- EFA.dimensions::SMT(
            data = response_data,
            corkind = corkind
          )
          
          # Stocker les résultats détaillés pour SMT
          results$SMT_details <- capture.output(print(smt_result))
          
          # Stocker la matrice des p-values pour l'affichage détaillé
          results$SMT_pvalues <- smt_result$pvalues
          
          # Stocker le nombre de facteurs suggéré directement
          final_results[["SMT"]] <- smt_result$NfactorsSMT
          
        }, error = function(e) {
          showNotification(paste("Erreur lors de l'exécution du test SMT:", e$message), 
                           type = "error", duration = 10)
        })
      }
      
      # Vérifier si des résultats ont été trouvés
      if (length(final_results) == 0) {
        final_results <- list(NOTE = "Les tests ont été exécutés mais les résultats n'ont pas pu être extraits.")
      }
      
      # Stocker les résultats
      results$final <- final_results
      dimtests_results(results)
      
      # Notification
      if ("NOTE" %in% names(final_results)) {
        showNotification(final_results$NOTE, type = "warning")
      } else {
        showNotification("Tests dimensionnels complétés !", type = "message")
      }
    })
  })


  
  # Afficher le résumé des tests
  output$dimtests_summary <- renderUI({
    req(sim_data(), dimtests_results())
    results <- dimtests_results()
    
    # Si pas de résultats finaux
    if (is.null(results$final)) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"), 
        strong("Tests non disponibles"),
        p("Aucun résultat de test dimensionnel n'est disponible.")
      ))
    }
    
    # Récupérer le nombre de facteurs extrait par chaque test
    test_results <- results$final
    
    # Vérifier s'il y a uniquement une note d'erreur
    if (length(test_results) == 1 && !is.null(test_results$NOTE)) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"), 
        strong("Tests non disponibles"),
        p(test_results$NOTE)
      ))
    }
    
    # Récupérer le nombre réel de facteurs (avec gestion des erreurs)
    n_factors <- tryCatch({
      # Utiliser input$nfact comme solution de secours si l'attribut n'est pas disponible
      factor_count <- attr(sim_data(), "n.factors")
      if (is.null(factor_count)) {
        factor_count <- input$nfact  # Fallback sur le paramètre utilisé pour générer les données
      }
      factor_count
    }, error = function(e) {
      # En cas d'erreur, utiliser input$nfact ou NA
      if (!is.null(input$nfact)) {
        return(input$nfact)
      }
      return(NA)
    })
    
    # Créer un tableau de résumé formaté HTML
    output <- div(
      h4("Résumé des tests dimensionnels"),
      p(
        "Nombre réel de facteurs dans les données générées: ",
        if(is.na(n_factors)) {
          tags$em("Non disponible", class = "text-muted")
        } else {
          tags$strong(n_factors)
        }
      ),
      tags$table(
        class = "table table-striped table-bordered",
        tags$thead(
          tags$tr(
            tags$th("Test"),
            tags$th("Nombre de facteurs suggéré"),
            tags$th("Cohérence") 
          )
        ),
        tags$tbody(
          lapply(names(test_results), function(test_name) {
            # Ignorer la note si elle existe
            if (test_name == "NOTE") return(NULL)
            
            result_value <- if(is.null(test_results[[test_name]]) || is.na(test_results[[test_name]])) {
              "Non disponible"
            } else {
              test_results[[test_name]]
            }
            
            # Vérifier la cohérence avec le nombre réel seulement si n_factors est connu
            consistency <- if(is.numeric(result_value) && !is.na(result_value) && !is.na(n_factors)) {
              if(result_value == n_factors) {
                list(class = "text-success", icon = icon("check"), text = "Correct")
              } else {
                list(class = "text-danger", icon = icon("times"), text = "Incorrect")
              }
            } else {
              list(class = "text-muted", icon = icon("question"), text = "Indéterminé")
            }
            
            tags$tr(
              tags$td(test_name),
              tags$td(result_value),
              tags$td(
                class = consistency$class,
                consistency$icon,
                consistency$text
              )
            )
          })
        )
      ),
      hr(),
      p(
        class = "text-muted",
        "Note: Le critère \"Valeurs propres > 1\" est souvent considéré comme insuffisant à lui seul."
      )
    )
    
    return(output)
  }) 
  
  # Observer pour l'analyse de factorabilité indépendante
  observeEvent(input$run_factorability, {
    req(sim_data())
    
    # Afficher un indicateur de chargement
    withProgress(message = 'Analyse de factorabilité en cours...', {
      # Récupérer les données
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Configurer les paramètres
      corkind <- input$factor_cor
      
      # Capturer la sortie console de la fonction FACTORABILITY
      output_capture <- capture.output({
        fact_result <- EFA.dimensions::FACTORABILITY(
          data = response_data,
          corkind = corkind
        )
      })
      
      # Stocker les résultats
      factorability_results(list(
        result = fact_result,
        console_output = paste(output_capture, collapse = "\n")
      ))
      
      # Notification
      showNotification("Analyse de factorabilité complétée !", type = "message")
    })
  })
  
  
  # Output pour l'affichage indépendant des résultats de factorabilité
  output$factorability_main_results <- renderUI({
    req(factorability_results())
    result <- factorability_results()$result
    
    # Si pas de résultats
    if (is.null(result)) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"), 
        strong("Analyse non disponible"),
        p("Les résultats de l'analyse de factorabilité ne sont pas disponibles.")
      ))
    }
    
    # Vérifier soigneusement la présence de chaque élément nécessaire
    has_overall_msa <- !is.null(result$Overall_MSA) && is.numeric(result$Overall_MSA)
    has_bartlett <- !is.null(result$bartlett) && is.list(result$bartlett) && 
      !is.null(result$bartlett$chisq) && !is.null(result$bartlett$df) && !is.null(result$bartlett$p.value)
    has_det_corr <- !is.null(result$det_corr) && is.numeric(result$det_corr)
    # Ajouter cette ligne importante qui manquait
    has_msa_values <- !is.null(result$MSA_values) && is.vector(result$MSA_values) && length(result$MSA_values) > 0
    
    # Formatage des valeurs avec vérification préalable
    overall_msa <- if(has_overall_msa) round(result$Overall_MSA, 3) else "Non disponible"
    bartlett_chisq <- if(has_bartlett) round(result$bartlett$chisq, 2) else "Non disponible"
    bartlett_df <- if(has_bartlett) result$bartlett$df else "Non disponible"
    bartlett_p <- if(has_bartlett) result$bartlett$p.value else NA
    det_corr <- if(has_det_corr) format(result$det_corr, scientific = TRUE, digits = 4) else "Non disponible"
    
    # Créer un tableau de résumé formaté HTML
    output <- div(
      h4("Analyse de la factorabilité des données"),
      
      # KMO & Bartlett
      card(
        card_header("Tests généraux de factorabilité"),
        card_body(
          tags$table(
            class = "table table-bordered table-sm",
            tags$tbody(
              tags$tr(
                tags$th("Mesure KMO globale"),
                tags$td(
                  style = if(has_overall_msa) {
                    ifelse(result$Overall_MSA >= 0.8, "color: green; font-weight: bold;",
                           ifelse(result$Overall_MSA >= 0.6, "color: orange;", "color: red;"))
                  } else { "" },
                  overall_msa,
                  if(has_overall_msa) {
                    ifelse(result$Overall_MSA >= 0.8, " (Excellente)",
                           ifelse(result$Overall_MSA >= 0.7, " (Bonne)",
                                  ifelse(result$Overall_MSA >= 0.6, " (Acceptable)",
                                         ifelse(result$Overall_MSA >= 0.5, " (Médiocre)", " (Inacceptable)"))))
                  } else { " (Non évaluable)" }
                )
              ),
              tags$tr(
                tags$th("Test de sphéricité de Bartlett"),
                tags$td(
                  if(has_bartlett) {
                    paste("χ² =", bartlett_chisq, 
                          ", df =", bartlett_df, 
                          ", p", ifelse(bartlett_p < 0.001, " < 0.001", 
                                        paste0(" = ", round(bartlett_p, 3))))
                  } else {
                    "Résultats incomplets ou non disponibles"
                  }
                )
              ),
              tags$tr(
                tags$th("Déterminant de la matrice de corrélation"),
                tags$td(
                  style = if(has_det_corr) {
                    ifelse(result$det_corr < 0.00001, "color: red;",
                           ifelse(result$det_corr < 0.001, "color: orange;", "color: green;"))
                  } else { "" },
                  det_corr,
                  if(has_det_corr) {
                    ifelse(result$det_corr < 0.00001, " (Multicolinéarité potentielle)",
                           ifelse(result$det_corr < 0.001, " (Acceptable pour l'analyse factorielle)", 
                                  " (Bon)"))
                  } else { " (Non évaluable)" }
                )
              )
            )
          )
        )
      ),
      
      hr(),
      
      # Variables problématiques (si présentes et si MSA_values est disponible)
      if (has_msa_values && any(result$MSA_values < 0.5, na.rm = TRUE)) {
        card(
          card_header(
            HTML("<span style='color: red;'><i class='fa fa-exclamation-triangle'></i> Variables problématiques</span>")
          ),
          card_body(
            p("Les variables suivantes ont des indices KMO inférieurs à 0.5, ce qui peut indiquer qu'elles ne sont pas adaptées à l'analyse factorielle :"),
            tags$ul(
              lapply(names(which(result$MSA_values < 0.5)), function(var_name) {
                tags$li(
                  HTML(paste0("<strong>", var_name, "</strong>: KMO = ", 
                              round(result$MSA_values[var_name], 3)))
                )
              })
            ),
            p("Il peut être intéressant de considérer l'exclusion de ces variables.")
          )
        )
      } else if (has_msa_values) {
        div(
          class = "alert alert-success",
          icon("check-circle"),
          strong("Toutes les variables sont adaptées à l'analyse factorielle"),
          p("Tous les indices KMO spécifiques aux variables sont supérieurs à 0.5.")
        )
      } else {
        div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          strong("Impossible d'évaluer les variables individuelles"),
          p("Les indices KMO spécifiques aux variables ne sont pas disponibles.")
        )
      },
      
      hr(),
      
      # Interprétation générale (uniquement si les valeurs requises sont disponibles)
      if (has_overall_msa && has_bartlett) {
        div(
          class = ifelse(result$Overall_MSA >= 0.6 && result$bartlett$p.value < 0.05, 
                         "alert alert-success", "alert alert-warning"),
          icon(ifelse(result$Overall_MSA >= 0.6 && result$bartlett$p.value < 0.05, 
                      "thumbs-up", "exclamation-triangle")),
          strong("Conclusion générale :"),
          p(HTML(ifelse(result$Overall_MSA >= 0.6 && result$bartlett$p.value < 0.05,
                        "Les données sont <strong>adaptées</strong> pour une analyse factorielle.",
                        "Les données présentent des <strong>problèmes potentiels</strong> pour l'analyse factorielle.")))
        )
      } else {
        div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          strong("Conclusion générale :"),
          p("Impossible de déterminer la factorabilité des données en raison de valeurs manquantes dans l'analyse.")
        )
      }
    )
    
    return(output)
  })
  

  # Affichage brut des résultats de factorabilité
  output$factorability_raw_results <- renderPrint({
    req(sim_data())
    
    # Exécuter FACTORABILITY directement pour obtenir sa sortie console complète
    data_obj <- sim_data()
    response_data <- data_obj$response  
    
    cat("RÉSULTATS DÉTAILLÉS DE L'ANALYSE DE FACTORABILITÉ\n")
    cat("==============================================\n")
    result <- EFA.dimensions::FACTORABILITY(
      data = response_data,
      corkind = input$factor_cor,
      verbose = FALSE
    )
    
    # Vérifier correctement la présence des éléments avant de les utiliser
    has_kmo <- !is.null(result$KMO)
    has_pvalue <- !is.null(result$pvalue)
    
    # 1. Test de Bartlett
    cat("\nTEST DE BARTLETT\n")
    cat("---------------\n")
    cat("Chi-carré: ", if(!is.null(result$chisq)) result$chisq else "Non disponible", "\n")
    cat("Degrés de liberté: ", if(!is.null(result$df)) result$df else "Non disponible", "\n")
    cat("Valeur p: ", if(!is.null(result$pvalue)) result$pvalue else "Non disponible", "\n\n")
    
    # 2. Calculer et afficher le déterminant
    cat("DÉTERMINANT DE LA MATRICE DE CORRÉLATION\n")
    cat("--------------------------------------\n")
    corr_matrix <- cor(response_data, method = input$factor_cor)
    det_corr <- det(corr_matrix)
    cat("Valeur: ", format(det_corr, scientific = TRUE, digits = 4), "\n")
    cat("Interprétation: ", ifelse(det_corr > 0.00001, 
                                   "Factorabilité indiquée (> 0.00001)", 
                                   "Multicolinéarité potentielle (< 0.00001)"), "\n\n")
    
    
    # 3. KMO global
    cat("MESURE KMO GLOBALE\n")
    cat("----------------\n")
    if(has_kmo) {
      cat("Valeur: ", result$KMO, "\n")
      cat("Interprétation: ", 
          ifelse(result$KMO >= 0.9, "Merveilleuse (≥ 0.9)",
                 ifelse(result$KMO >= 0.8, "Méritoire (0.8-0.89)",
                        ifelse(result$KMO >= 0.7, "Moyenne (0.7-0.79)",
                               ifelse(result$KMO >= 0.6, "Médiocre (0.6-0.69)",
                                      ifelse(result$KMO >= 0.5, "Misérable (0.5-0.59)",
                                             "Inacceptable (< 0.5)"))))), "\n\n")
    } else {
      cat("Valeur: Non disponible\n\n")
    }

    # 5. KMO par variable
    cat("MESURES KMO PAR VARIABLE\n")
    cat("----------------------\n")
    if(!is.null(result$KMOvars)) {
      print(result$KMOvars)
    } else {
      cat("Non disponibles\n\n")
    }
    
    # Conclusion générale - Maintenant nous utilisons det_corr au lieu de result$det_corr
    cat("\nCONCLUSION GÉNÉRALE\n")
    cat("-----------------\n")
    if(has_kmo && has_pvalue) {
      if(result$KMO >= 0.6 && result$pvalue < 0.05 && det_corr > 0.00001) {
        cat("Les données sont adéquates pour l'analyse factorielle.\n")
      } else {
        cat("Les données présentent des caractéristiques qui pourraient compromettre l'analyse factorielle.\n")
        if(result$KMO < 0.6) cat("- KMO insuffisant (< 0.6)\n")
        if(result$pvalue >= 0.05) cat("- Test de Bartlett non significatif\n")
        if(det_corr <= 0.00001) cat("- Déterminant trop faible (multicolinéarité possible)\n")
      }
    } else {
      cat("Impossible de déterminer complètement la factorabilité des données en raison de résultats incomplets.\n")
    }
  })
  
  output$analysis_result_title <- renderUI({
    if(input$analysis_type == "efa") {
      "Résultats de l'analyse factorielle exploratoire"
    } else {
      "Résultats de l'analyse en composantes principales"
    }
  })
  
  # Observer pour l'analyse de cohérence interne
  observeEvent(input$run_internal_consistency, {
    req(sim_data())
    
    # Afficher un indicateur de chargement
    withProgress(message = 'Analyse de cohérence interne en cours...', {
      # Récupérer les données
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Configurer les paramètres
      extraction_method <- input$ic_extraction
      
      # Exécuter l'analyse avec gestion d'erreur
      tryCatch({
        # Capturer la sortie console
        output_capture <- capture.output({
          ic_result <- EFA.dimensions::INTERNAL.CONSISTENCY(
            data = response_data,
            extraction = extraction_method,
            verbose = FALSE  
          )
        })
        
        # Stocker les résultats pour affichage formaté
        internal_consistency_results(ic_result)
        
        # Notification
        showNotification("Analyse de cohérence interne complétée !", type = "message")
        
      }, error = function(e) {
        # En cas d'erreur
        showNotification(paste("Erreur lors de l'analyse de cohérence interne:", e$message), 
                         type = "error", 
                         duration = 10)
      })
    })
  })

  # Fonction pour formater les résultats de cohérence interne
  output$internal_consistency_results <- renderPrint({
    req(sim_data())
    
    # Vérifier si l'utilisateur a cliqué sur le bouton d'exécution
    if (is.null(internal_consistency_results())) {
      cat("Cliquez sur 'Exécuter l'analyse' pour calculer les indices de cohérence interne.\n")
      return(invisible())
    }
    
    # Récupérer les résultats
    ic_result <- internal_consistency_results()
    
    # Afficher les résultats dans un format structuré
    cat("RÉSULTATS DE L'ANALYSE DE COHÉRENCE INTERNE\n")
    cat("==========================================\n\n")
    
    # Vérifier que int.consist_scale existe
    if (!is.null(ic_result$int.consist_scale)) {
      # Extraire les valeurs correctement - en utilisant la syntaxe pour les matrices
      omega_value <- ic_result$int.consist_scale[1, "omega"]
      alpha_value <- ic_result$int.consist_scale[1, "alpha"]
      alpha_z_value <- ic_result$int.consist_scale[1, "alpha.z"]
      r_mean_value <- ic_result$int.consist_scale[1, "r_mean"]
      r_median_value <- ic_result$int.consist_scale[1, "r_median"]
      rmsr_value <- ic_result$int.consist_scale[1, "rmsr"]
      
      # Alpha de Cronbach
      cat("ALPHA DE CRONBACH\n")
      cat("----------------\n")
      cat("Alpha:", round(alpha_value, 3), "\n")
      cat("Alpha standardisé:", round(alpha_z_value, 3), "\n")
      
      # Interpretation standard de l'alpha
      cat("Interprétation: ", 
          ifelse(alpha_value >= 0.9, "Excellent (≥ 0.9)",
                 ifelse(alpha_value >= 0.8, "Bon (0.8-0.89)",
                        ifelse(alpha_value >= 0.7, "Acceptable (0.7-0.79)",
                               ifelse(alpha_value >= 0.6, "Questionnable (0.6-0.69)",
                                      ifelse(alpha_value >= 0.5, "Faible (0.5-0.59)",
                                             "Inacceptable (< 0.5)"))))), "\n\n")
      
      # Oméga
      cat("COEFFICIENT OMÉGA\n")
      cat("-----------------\n")
      cat("Oméga total:", round(omega_value, 3), "\n\n")
      
      # Autres statistiques
      cat("AUTRES INDICES\n")
      cat("-------------\n")
      cat("Moyenne des corrélations inter-items:", round(r_mean_value, 3), "\n")
      cat("Médiane des corrélations inter-items:", round(r_median_value, 3), "\n")
      cat("RMSR (Root Mean Square Residual):", round(rmsr_value, 3), "\n\n")
      
      # Afficher également les statistiques détaillées par item
      if (!is.null(ic_result$int.consist_dropped)) {
        cat("STATISTIQUES DÉTAILLÉES PAR ITEM\n")
        cat("-------------------------------\n")
        cat("Indices de cohérence interne lorsque l'item est supprimé:\n\n")
        
        # En-tête du tableau
        cat(sprintf("%-15s  %8s  %8s  %8s  %8s  %8s  %8s\n", 
                    "Item", "Omega", "Alpha", "Alpha-z", "r_mean", "r_median", "RMSR"))
        cat(paste(rep("-", 80), collapse = ""), "\n")
        
        # Afficher les valeurs pour chaque item
        for (i in 1:nrow(ic_result$int.consist_dropped)) {
          item_name <- rownames(ic_result$int.consist_dropped)[i]
          omega_drop <- ic_result$int.consist_dropped[i, "omega"]
          alpha_drop <- ic_result$int.consist_dropped[i, "alpha"]
          alpha_z_drop <- ic_result$int.consist_dropped[i, "alpha.z"]
          r_mean_drop <- ic_result$int.consist_dropped[i, "r_mean"]
          r_median_drop <- ic_result$int.consist_dropped[i, "r_median"]
          rmsr_drop <- ic_result$int.consist_dropped[i, "rmsr"]
          
          cat(sprintf("%-15s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n", 
                      item_name, omega_drop, alpha_drop, alpha_z_drop, 
                      r_mean_drop, r_median_drop, rmsr_drop))
        }
        cat("\n")
      }
      
      # Conclusion
      cat("CONCLUSION\n")
      cat("----------\n")
      
      if (alpha_value >= 0.7 || omega_value >= 0.7) {
        cat("L'échelle présente une cohérence interne satisfaisante.\n")
      } else if (alpha_value >= 0.6 || omega_value >= 0.6) {
        cat("L'échelle présente une cohérence interne acceptable mais pourrait être améliorée.\n")
      } else {
        cat("L'échelle présente une cohérence interne insuffisante pour la plupart des applications.\n")
      }
      
      # Recommandations
      if (!is.null(ic_result$int.consist_dropped)) {
        # Identifier les items qui amélioreraient l'alpha
        better_items <- names(which(alpha_drop > alpha_value))
        
        if (length(better_items) > 0) {
          cat("\nRECOMMANDATIONS\n")
          cat("--------------\n")
          cat("La suppression des items suivants pourrait améliorer la cohérence interne :\n")
          
          for (item in better_items) {
            old_alpha <- alpha_value
            new_alpha <- alpha_drop[item]
            improvement <- (new_alpha - old_alpha)
            cat(sprintf("- %s: Augmentation de l'alpha de %.3f à %.3f (+%.3f)\n", 
                        item, old_alpha, new_alpha, improvement))
          }
        }
      }
    } else {
      cat("Aucun résultat de cohérence interne n'est disponible.\n")
      cat("Vérifiez que les données sont appropriées pour cette analyse.\n")
    }
  })
  
  # Afficher le plot des résultats avec plotly
  output$dimtests_plot <- renderPlotly({
    req(dimtests_results())
    results <- dimtests_results()
    
    # Récupérer les valeurs propres
    # Utiliser les valeurs propres appropriées selon la méthode d'extraction
    if (!is.null(results$parallel_data) && !is.null(results$parallel_data$observed)) {
      # Utiliser les eigenvalues spécifiques à la méthode d'extraction (PAF, PCA, etc.)
      eigenvalues <- results$parallel_data$observed
    } else {
      # Fallback sur les eigenvalues standards de la matrice de corrélation
      eigenvalues <- results$eigenvalues
    }
    
    # Créer un dataframe pour le plot
    n_factors <- length(eigenvalues)
    factor_numbers <- 1:n_factors
    
    # Initialiser le plotly avec un dataframe explicite
    p <- plot_ly(data.frame(x = factor_numbers, y = eigenvalues), 
                 x = ~x, y = ~y, height = 600) %>%
      # Ajouter la trace pour les eigenvalues observées
      add_trace(
        type = "scatter", 
        mode = "lines+markers",
        name = "Eigenvalues observées",
        line = list(color = "blue", width = 3),
        marker = list(color = "red", size = 8),
        hoverinfo = "text",
        text = paste("Facteur:", factor_numbers, "<br>Eigenvalue:", round(eigenvalues, 3))
      )
    
    # Ajouter la ligne pour le critère de Kaiser
    kaiser_df <- data.frame(x = c(1, n_factors), y = c(1, 1))
    p <- p %>% add_trace(
      data = kaiser_df,
      x = ~x, y = ~y,
      type = "scatter", 
      mode = "lines",
      name = "Critère de Kaiser (λ = 1)",
      line = list(color = "darkgreen", width = 2, dash = "dash"),
      hoverinfo = "text",
      text = "Critère de Kaiser: eigenvalue = 1"
    )
    
    # Ajouter les données du test parallèle si disponibles
    if (!is.null(results$parallel_data)) {
      # Vérifier que les composants nécessaires sont disponibles
      if (!is.null(results$parallel_data$Mean)) {
        # Extraire les valeurs propres aléatoires (moyennes)
        rand_eigenvalues <- results$parallel_data$Mean
        
        # Limiter à la longueur des valeurs propres réelles
        rand_eigenvalues <- rand_eigenvalues[1:min(length(rand_eigenvalues), n_factors)]
        rand_factors <- 1:length(rand_eigenvalues)
        
        # Créer un dataframe pour les moyennes
        mean_df <- data.frame(x = rand_factors, y = rand_eigenvalues)
        
        # Ajouter la trace pour le test parallèle - moyennes
        p <- p %>% add_trace(
          data = mean_df,
          x = ~x, y = ~y,
          type = "scatter", 
          mode = "lines",
          name = "Test parallèle (moyennes)",
          line = list(color = "orange", width = 2, dash = "dot"),
          hoverinfo = "text",
          text = paste("Facteur:", rand_factors, 
                       "<br>Eigenvalue aléatoire (moyenne):", round(rand_eigenvalues, 3))
        )
      }
      
      # Ajouter également la ligne du percentile si disponible
      if (!is.null(results$parallel_data$Percentile)) {
        # Extraire les valeurs propres au percentile spécifié
        percentile_eigenvalues <- results$parallel_data$Percentile
        
        # Limiter à la longueur des valeurs propres réelles
        percentile_eigenvalues <- percentile_eigenvalues[1:min(length(percentile_eigenvalues), n_factors)]
        percentile_factors <- 1:length(percentile_eigenvalues)
        
        # Créer un dataframe pour les percentiles
        percentile_df <- data.frame(x = percentile_factors, y = percentile_eigenvalues)
        
        # Obtenir le percentile à partir de l'input
        percentile_value <- input$rawpar_percentile
        
        # Ajouter la trace pour le test parallèle - percentile
        p <- p %>% add_trace(
          data = percentile_df,
          x = ~x, y = ~y,
          type = "scatter", 
          mode = "lines",
          name = paste0("Test parallèle (", percentile_value, "e percentile)"),
          line = list(color = "red", width = 2, dash = "dashdot"),
          hoverinfo = "text",
          text = paste("Facteur:", percentile_factors, 
                       "<br>Eigenvalue aléatoire (percentile):", round(percentile_eigenvalues, 3))
        )
      }
    }
    
    # Ajouter le nombre réel de facteurs s'il est disponible
    n_factors_original <- NULL
    tryCatch({
      if (!is.null(sim_data()) && !is.null(attr(sim_data(), "n.factors"))) {
        n_factors_original <- attr(sim_data(), "n.factors")
      }
    }, error = function(e) {
      message("Erreur lors de l'accès à n.factors: ", e$message)
    })
    
    if (!is.null(n_factors_original)) {
      # Dataframe pour la ligne verticale
      real_factors_df <- data.frame(
        x = c(n_factors_original, n_factors_original), 
        y = c(0, max(eigenvalues))
      )
      
      p <- p %>% add_trace(
        data = real_factors_df,
        x = ~x, y = ~y,
        type = "scatter", 
        mode = "lines",
        name = paste("Nombre réel de facteurs:", n_factors_original),
        line = list(color = "purple", width = 2, dash = "dot"),
        hoverinfo = "text",
        text = paste("Nombre réel de facteurs:", n_factors_original)
      )
    }
    
    # Configuration du layout
    p <- p %>% layout(
      title = list(
        text = "Scree Plot avec résultats des tests",
        font = list(size = 22)
      ),
      xaxis = list(
        title = "Numéro de facteur",
        titlefont = list(size = 16),
        tickfont = list(size = 14),
        tickmode = "array",
        tickvals = 1:n_factors
      ),
      yaxis = list(
        title = "Valeur propre",
        titlefont = list(size = 16),
        tickfont = list(size = 14)
      ),
      legend = list(
        x = 0.5,                         # Déplacé plus à gauche (était 0.7)
        y = 0.9,                         # Position près du haut
        font = list(size = 16),          # Taille du texte de la légende
        bgcolor = "rgba(255, 255, 255, 0.8)", # Fond semi-transparent
        bordercolor = "rgba(0, 0, 0, 0.3)",  # Bordure légère
        borderwidth = 1
      ),
      hovermode = "closest",
      plot_bgcolor = "rgba(255, 255, 255, 1)",    # Fond blanc
      paper_bgcolor = "rgba(255, 255, 255, 1)",   # Fond blanc
      margin = list(t = 100, b = 80, l = 80, r = 40)  # Marges
    )
    
    # Configuration des outils
    p <- p %>% config(
      displayModeBar = TRUE,
      modeBarButtonsToRemove = c("lasso2d", "select2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = "screeplot",
        width = 1200,
        height = 800
      )
    )
    
    return(p)
  })
  
  # Afficher le résumé de factorabilité
  output$factorability_summary <- renderUI({
    req(factorability_results())
    result <- factorability_results()
    
    # Si pas de résultats
    if (is.null(result)) {
      return(div(
        class = "alert alert-warning",
        icon("exclamation-triangle"), 
        strong("Analyse non disponible"),
        p("Les résultats de l'analyse de factorabilité ne sont pas disponibles.")
      ))
    }
    
    # Créer un tableau de résumé formaté HTML
    output <- div(
      h4("Analyse de la factorabilité des données"),
      
      # KMO & Bartlett
      card(
        card_header("Tests généraux de factorabilité"),
        card_body(
          tags$table(
            class = "table table-bordered table-sm",
            tags$tbody(
              tags$tr(
                tags$th("Mesure KMO globale"),
                tags$td(
                  style = ifelse(result$Overall_MSA >= 0.8, "color: green; font-weight: bold;",
                                 ifelse(result$Overall_MSA >= 0.6, "color: orange;", "color: red;")),
                  round(result$Overall_MSA, 3),
                  ifelse(result$Overall_MSA >= 0.8, " (Excellente)",
                         ifelse(result$Overall_MSA >= 0.7, " (Bonne)",
                                ifelse(result$Overall_MSA >= 0.6, " (Acceptable)",
                                       ifelse(result$Overall_MSA >= 0.5, " (Médiocre)", " (Inacceptable)"))))
                )
              ),
              tags$tr(
                tags$th("Test de sphéricité de Bartlett"),
                tags$td(
                  paste("χ² =", round(result$bartlett$chisq, 2), 
                        ", df =", result$bartlett$df, 
                        ", p", ifelse(result$bartlett$p.value < 0.001, " < 0.001", 
                                      paste0(" = ", round(result$bartlett$p.value, 3))))
                )
              ),
              tags$tr(
                tags$th("Déterminant de la matrice de corrélation"),
                tags$td(
                  style = ifelse(result$det_corr < 0.00001, "color: red;",
                                 ifelse(result$det_corr < 0.001, "color: orange;", "color: green;")),
                  format(result$det_corr, scientific = TRUE, digits = 4),
                  ifelse(result$det_corr < 0.00001, " (Multicolinéarité potentielle)",
                         ifelse(result$det_corr < 0.001, " (Acceptable pour l'analyse factorielle)", 
                                " (Bon)"))
                )
              )
            )
          )
        )
      ),
      
      hr(),
      
      # Variables problématiques (si présentes)
      if (any(result$MSA_values < 0.5)) {
        card(
          card_header(
            HTML("<span style='color: red;'><i class='fa fa-exclamation-triangle'></i> Variables problématiques</span>")
          ),
          card_body(
            p("Les variables suivantes ont des indices KMO inférieurs à 0.5, ce qui peut indiquer qu'elles ne sont pas adaptées à l'analyse factorielle :"),
            tags$ul(
              lapply(names(which(result$MSA_values < 0.5)), function(var_name) {
                tags$li(
                  HTML(paste0("<strong>", var_name, "</strong>: KMO = ", 
                              round(result$MSA_values[var_name], 3)))
                )
              })
            ),
            p("Il peut être intéressant de considérer l'exclusion de ces variables.")
          )
        )
      } else {
        div(
          class = "alert alert-success",
          icon("check-circle"),
          strong("Toutes les variables sont adaptées à l'analyse factorielle"),
          p("Tous les indices KMO spécifiques aux variables sont supérieurs à 0.5.")
        )
      },
      
      hr(),
      
      # Interprétation générale
      div(
        class = ifelse(result$Overall_MSA >= 0.6 && result$bartlett$p.value < 0.05, 
                       "alert alert-success", "alert alert-warning"),
        icon(ifelse(result$Overall_MSA >= 0.6 && result$bartlett$p.value < 0.05, 
                    "thumbs-up", "exclamation-triangle")),
        strong("Conclusion générale :"),
        p(HTML(ifelse(result$Overall_MSA >= 0.6 && result$bartlett$p.value < 0.05,
                      "Les données sont <strong>adaptées</strong> pour une analyse factorielle.",
                      "Les données présentent des <strong>problèmes potentiels</strong> pour l'analyse factorielle.")))
      )
    )
    
    return(output)
  })
  
  output$factor_structure_detailed <- renderUI({
    req(sim_data())
    data_obj <- sim_data()
    
    # Récupérer les informations nécessaires
    loadings_matrix <- data_obj$loadings
    factor_names <- colnames(loadings_matrix)
    var_names <- rownames(loadings_matrix)
    
    # Déterminer les variables associées à chaque facteur
    # en utilisant la charge la plus élevée
    factor_vars_list <- list()
    for (i in 1:ncol(loadings_matrix)) {
      # Trouver les variables qui chargent principalement sur ce facteur
      max_loadings <- apply(loadings_matrix, 1, function(row) which.max(abs(row)))
      factor_vars_list[[i]] <- which(max_loadings == i)
    }
    
    html_content <- tagList()
    
    for (i in 1:length(factor_names)) {
      factor_name <- factor_names[i]
      variables <- factor_vars_list[[i]]
      
      # Ne créer une carte que s'il y a des variables pour ce facteur
      if (length(variables) > 0) {
        factor_card <- card(
          card_header(
            tags$span(
              tags$span(class = "badge bg-primary", style = "font-size: 1rem; margin-right: 10px;", paste0("F", i)),
              tags$strong(paste0("Facteur ", i))
            )
          ),
          card_body(
            tags$ul(
              style = "list-style-type: none; padding-left: 0;",
              lapply(variables, function(var_idx) {
                var_name <- var_names[var_idx]
                loading_value <- round(loadings_matrix[var_idx, i], 3)
                loading_class <- ifelse(abs(loading_value) >= 0.6, "text-success", 
                                        ifelse(abs(loading_value) >= 0.4, "text-primary", "text-muted"))
                
                tags$li(
                  style = "margin-bottom: 5px; border-bottom: 1px solid #eee; padding-bottom: 5px;",
                  tags$div(
                    tags$span(var_name, style = "font-weight: bold;"),
                    tags$span(class = loading_class, style = "float: right;", 
                              paste0("λ = ", loading_value))
                  )
                )
              })
            )
          )
        )
        
        html_content <- tagAppendChild(html_content, factor_card)
        html_content <- tagAppendChild(html_content, tags$br())
      }
    }
    
    return(html_content)
  })

  # Fonction pour l'analyse factorielle exploratoire (EFA)
  output$efa_results <- renderPrint({
    req(sim_data())
    
    # Vérifier si l'utilisateur a cliqué sur le bouton d'exécution
    if (is.null(efa_result())) {
      cat("Cliquez sur 'Exécuter l'analyse' pour effectuer l'analyse factorielle exploratoire.\n")
      return(invisible())
    }
    
    # Récupérer les résultats
    results <- efa_result()
    
    # Déterminer si c'est une rotation orthogonale ou oblique
    is_none <- input$rotation_type == "none"
    is_orthogonal <- input$rotation_type == "orthogonal"
    
    # Afficher les résultats dans un format structuré
    cat("RÉSULTATS DE L'ANALYSE FACTORIELLE EXPLORATOIRE\n")
    cat("=============================================\n\n")
    
    # Méthode d'extraction et de rotation
    cat("MÉTHODE D'EXTRACTION ET DE ROTATION\n")
    cat("-----------------------------------\n")
    cat("Extraction: ", input$extraction_method, "\n")
    rotation <- if(is_none){"none"
      } else {input$rotation_method}
    cat("Rotation: ", rotation, "\n\n")
    
    # Nombre de facteurs
    cat("NOMBRE DE FACTEURS\n")
    cat("----------------\n")
    
    # Déterminer le nombre de facteurs à partir des résultats disponibles
    if (input$factor_method == "auto") {
      cat("Nombre de facteurs (EMPCK): ", results$nfactors, "\n\n")
    } else {
      cat("Nombre de facteurs spécifié: ", input$num_factors, "\n\n")
    }
    
    # Afficher les valeurs propres initiales et la variance expliquée (varexplNOROT1)
    if (!is.null(results$varexplNOROT1) && is.matrix(results$varexplNOROT1)) {
      cat("VALEURS PROPRES INITIALES ET VARIANCE EXPLIQUÉE\n")
      cat(paste(rep("-", 60), collapse = ""), "\n")
      
      varexpl <- results$varexplNOROT1
      
      # En-tête du tableau
      cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                  "Facteur", "Valeur propre", "% de variance", "% cumulatif"))
      cat(paste(rep("-", 60), collapse = ""), "\n")
      
      # Contenu du tableau
      for (i in 1:nrow(varexpl)) {
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                    paste0("Facteur ", i),
                    varexpl[i, 1],  # Eigenvalues
                    varexpl[i, 2] * 100,  # Proportion de variance
                    varexpl[i, 3] * 100   # Proportion cumulative
        ))
      }
      cat("\n")
    } 
    
    # Saturations avant rotation (si disponibles)
    if (!is.null(results$loadingsNOROT) && is.matrix(results$loadingsNOROT)) {
      cat("MATRICE DES SATURATIONS AVANT ROTATION\n")
      sep_length <- 15 + 10*ncol(results$loadingsNOROT)
      cat(paste(rep("-", sep_length), collapse = ""), "\n")

      # En-tête avec les noms des facteurs
      cat(sprintf("%-15s", "Variable"))
      for (i in 1:ncol(results$loadingsNOROT)) {
        cat(sprintf("  %8s", paste0("F", i)))
      }
      cat("\n")
      
      # Ligne de séparation
      cat(paste(rep("-", sep_length), collapse = ""), "\n")
      
      # Afficher les saturations pour chaque variable
      for (i in 1:nrow(results$loadingsNOROT)) {
        var_name <- rownames(results$loadingsNOROT)[i]
        if (is.null(var_name)) var_name <- paste0("Var", i)
        
        cat(sprintf("%-15s", var_name))
        for (j in 1:ncol(results$loadingsNOROT)) {
          loading <- results$loadingsNOROT[i, j]
          cat(sprintf("  %8.3f", loading))
        }
        cat("\n")
      }
      cat("\n")
    }
    
    # Afficher la variance expliquée après extraction (sans rotation) (varexplNOROT2)
    if (!is.null(results$varexplNOROT2) && is.matrix(results$varexplNOROT2)) {
      cat("VARIANCE EXPLIQUÉE APRÈS EXTRACTION (SANS ROTATION)\n")
      cat(paste(rep("-", 60), collapse = ""), "\n")
      
      varexpl <- results$varexplNOROT2
      
      # En-tête du tableau
      cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                  "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
      cat(paste(rep("-", 60), collapse = ""), "\n")
      
      # Contenu du tableau
      for (i in 1:nrow(varexpl)) {
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                    paste0("Facteur ", i),
                    varexpl[i, 1],  # Sommes des carrés
                    varexpl[i, 2] * 100,  # Proportion de variance
                    varexpl[i, 3] * 100   # Proportion cumulative
        ))
      }
      cat("\n")
    }
    
    # Afficher la variance expliquée après rotation (varexplROT)
    if (!is.null(results$varexplROT)) {
      cat("VARIANCE EXPLIQUÉE APRÈS ROTATION\n")
      cat("---------------------------------\n")
      
      # Pour les rotations obliques, nous devons utiliser la matrice de structure
      if (!is.null(results$structure) && is.matrix(results$structure)) {
        # Calculer la variance expliquée à partir de la matrice de structure
        structure_matrix <- results$structure
        
        # Calculer les sommes des carrés par facteur (somme des carrés des corrélations)
        ss_loadings <- colSums(structure_matrix^2)
        
        # Nombre total de variables pour le calcul des pourcentages
        total_variance <- nrow(structure_matrix)  # Nombre total de variables
        
        # En-tête du tableau
        cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                    "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
        cat(paste(rep("-", 60), collapse = ""), "\n")
        
        # Contenu du tableau - pourcentages par rapport au nombre total de variables
        cumulative_var <- 0
        for (i in 1:length(ss_loadings)) {
          var_percent <- (ss_loadings[i] / total_variance) * 100
          cumulative_var <- cumulative_var + var_percent
          
          cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                      paste0("Facteur ", i),
                      ss_loadings[i],
                      var_percent,
                      cumulative_var))
        }
        
        # Ajouter une ligne pour le total
        total_ss <- sum(ss_loadings)
        total_percent <- (total_ss / total_variance) * 100
        cat(paste(rep("-", 60), collapse = ""), "\n")
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14s\n", 
                    "Total",
                    total_ss,
                    total_percent,
                    ""))
        
        # Ajouter une note sur l'utilisation de la matrice de structure
        cat("\nNote: Pour la rotation oblique, les variances sont calculées à partir de la matrice de structure\n")
        cat("(corrélations entre variables et facteurs). En raison de la non-orthogonalité des facteurs,\n")
        cat("les pourcentages peuvent ne pas s'additionner exactement.\n")
      }
      # Sinon, utiliser les valeurs disponibles
      else if (is.matrix(results$varexplROT)) {
        # Format matriciel (rotation orthogonale)
        varexpl <- results$varexplROT
        
        # En-tête du tableau
        cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                    "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
        cat(paste(rep("-", 60), collapse = ""), "\n")
        
        # Contenu du tableau - pourcentages déjà calculés dans la matrice
        for (i in 1:nrow(varexpl)) {
          cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                      paste0("Facteur ", i),
                      varexpl[i, 1],  # Sommes des carrés
                      varexpl[i, 2] * 100,  # Proportion de variance
                      varexpl[i, 3] * 100   # Proportion cumulative
          ))
        }
      } 
      else if (is.vector(results$varexplROT)) {
        # Format vectoriel - dans ce cas, on a seulement les SS loadings
        varexpl <- results$varexplROT
        total_variance <- nrow(results$pattern) # Nombre total de variables
        
        # En-tête du tableau
        cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                    "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
        cat(paste(rep("-", 60), collapse = ""), "\n")
        
        # Contenu du tableau - calculer les pourcentages à partir des SS loadings
        cumulative_var <- 0
        for (i in 1:length(varexpl)) {
          var_percent <- (varexpl[i] / total_variance) * 100
          cumulative_var <- cumulative_var + var_percent
          
          cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                      paste0("Facteur ", i),
                      varexpl[i],
                      var_percent,
                      cumulative_var))
        }
        
        # Ajouter une ligne pour le total
        total_ss <- sum(varexpl)
        total_percent <- (total_ss / total_variance) * 100
        cat(paste(rep("-", 60), collapse = ""), "\n")
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14s\n", 
                    "Total",
                    total_ss,
                    total_percent,
                    ""))
        
        # Ajouter une note si la matrice de structure n'est pas disponible
        cat("\nNote: Les pourcentages sont calculés en divisant les sommes des carrés par le nombre total de variables.\n")
        cat("Pour une rotation oblique, une estimation plus précise nécessiterait la matrice de structure.\n")
      }
      
      cat("\n")
    }
    
    # Ratio des deux premières valeurs propres (utile pour évaluer l'unidimensionnalité)
    if (!is.null(results$evals12ratio)) {
      cat("RATIO DES DEUX PREMIÈRES VALEURS PROPRES\n")
      cat("----------------------------------------\n")
      cat(sprintf("Ratio: %.4f\n", results$evals12ratio))
      
      # Interprétation du ratio
      if (results$evals12ratio > 4) {
        cat("(Ratio > 4: Forte indication d'unidimensionnalité)\n")
      } else if (results$evals12ratio > 3) {
        cat("(Ratio > 3: Indication modérée d'unidimensionnalité)\n")
      } else if (results$evals12ratio > 2) {
        cat("(Ratio > 2: Faible indication d'unidimensionnalité)\n")
      } else {
        cat("(Ratio ≤ 2: Structure probablement multidimensionnelle)\n")
      }
      cat("\n")
    }
    
    # Déterminer quelle matrice de saturations utiliser
    loadings_matrix <- NULL
    if (input$rotation_type == "none") {
      loadings_matrix <- results$loadingsNOROT
      cat("MATRICE DES SATURATIONS (SANS ROTATION)\n")
    } else if (input$rotation_type == "orthogonal") {
      loadings_matrix <- results$loadingsROT
      cat("MATRICE DES SATURATIONS (ROTATION ORTHOGONALE)\n")
    } else if (input$rotation_type == "oblique") {
      loadings_matrix <- results$pattern
      cat("MATRICE DES SATURATIONS (PATTERN)\n")
    }
    
    # Si aucune des matrices ci-dessus n'est disponible, afficher un message et sortir
    if (is.null(loadings_matrix)) {
      cat("MATRICE DES SATURATIONS NON DISPONIBLE\n")
      return()  # Sortir de la fonction pour éviter les erreurs
    }
    
    # Vérifier explicitement que loadings_matrix est une matrice valide
    if (!is.matrix(loadings_matrix) || ncol(loadings_matrix) <= 0) {
      cat("MATRICE DES SATURATIONS MAL FORMATÉE\n")
      cat("Type:", class(loadings_matrix), "\n")
      if (is.matrix(loadings_matrix)) {
        cat("Dimensions:", nrow(loadings_matrix), "×", ncol(loadings_matrix), "\n")
      }
      return()  # Sortir de la fonction pour éviter les erreurs
    }
    
    # Afficher la matrice des saturations si disponible
    if (!is.null(loadings_matrix)) {
      sep_length <- 15 + 10*ncol(loadings_matrix) 
      cat(paste(rep("-", sep_length), collapse = ""), "\n")
      
      # En-tête avec les noms des facteurs
      cat(sprintf("%-15s", "Variable"))
      for (i in 1:ncol(loadings_matrix)) {
        cat(sprintf("  %8s", paste0("F", i)))
      }
      if (!is.null(results$h2)) {
        cat(sprintf("  %8s", "h²"))
      }
      cat("\n")
      
      # Ligne de séparation
      sep_length <- 15 + 10*ncol(loadings_matrix) 
      cat(paste(rep("-", sep_length), collapse = ""), "\n")
      
      # Récupérer le seuil des saturations à afficher depuis l'interface
      threshold <- input$loading_threshold
      
      # Afficher les saturations pour chaque variable
      for (i in 1:nrow(loadings_matrix)) {
        var_name <- rownames(loadings_matrix)[i]
        if (is.null(var_name)) var_name <- paste0("Var", i)
        
        cat(sprintf("%-15s", var_name))
        for (j in 1:ncol(loadings_matrix)) {
          # Mettre en évidence les saturations importantes (> threshold en valeur absolue)
          loading <- loadings_matrix[i, j]
          if (abs(loading) >= threshold) {
            cat(sprintf("  %8.3f", loading))
          } else {
            cat(sprintf("  %8s", ""))  # Laisser vide si < threshold
          }
        }
        if (!is.null(results$h2)) {
          cat(sprintf("  %8.3f", results$h2[i]))
        }
        cat("\n")
      }
      cat("\n")
      
      # Ajouter une note sur le seuil appliqué
      cat(sprintf("Note: Seules les saturations factorielles dont la valeur absolue est ≥ %.3f sont affichées.\n\n", threshold))
    }
    
    # Afficher la variance expliquée après rotation (varexplROT)
    if (!is.null(results$varexplROT)) {
      cat("VARIANCE EXPLIQUÉE APRÈS ROTATION\n")
      cat("---------------------------------\n")
      
      if (is.matrix(results$varexplROT)) {
        # Format matriciel
        varexpl <- results$varexplROT
        total_variance <- nrow(results$loadingsNOROT)
        
        # En-tête du tableau
        cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                    "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
        cat(paste(rep("-", 60), collapse = ""), "\n")
        
        # Contenu du tableau - pourcentages par rapport au nombre total de variables
        cumulative_var <- 0
        for (i in 1:nrow(varexpl)) {
          ss_loading <- varexpl[i, 1]  # Sommes des carrés
          var_percent <- (ss_loading / total_variance) * 100  # Division par le nombre de variables
          cumulative_var <- cumulative_var + var_percent
          
          cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                      paste0("Facteur ", i),
                      ss_loading,
                      var_percent,
                      cumulative_var))
        }
        
        # Ajouter une ligne pour le total
        total_ss <- sum(varexpl[, 1])
        total_percent <- (total_ss / total_variance) * 100
        cat(paste(rep("-", 60), collapse = ""), "\n")
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14s\n", 
                    "Total",
                    total_ss,
                    total_percent,
                    ""))
        
      } else if (is.vector(results$varexplROT)) {
        # Format vectoriel
        varexpl <- results$varexplROT
        total_variance <- ncol(data)  # Nombre total de variables = variance totale
        
        # En-tête du tableau
        cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                    "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
        cat(paste(rep("-", 60), collapse = ""), "\n")
        
        # Contenu du tableau - pourcentages par rapport au nombre total de variables
        cumulative_var <- 0
        for (i in 1:length(varexpl)) {
          ss_loading <- varexpl[i]
          var_percent <- (ss_loading / total_variance) * 100  # Division par le nombre de variables
          cumulative_var <- cumulative_var + var_percent
          
          cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                      paste0("Facteur ", i),
                      ss_loading,
                      var_percent,
                      cumulative_var))
        }
        
        # Ajouter une ligne pour le total
        total_ss <- sum(varexpl)
        total_percent <- (total_ss / total_variance) * 100
        cat(paste(rep("-", 60), collapse = ""), "\n")
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14s\n", 
                    "Total",
                    total_ss,
                    total_percent,
                    ""))
      }
      
      cat("\n")
    }
    
    # Afficher la matrice des corrélations entre facteurs (phi) pour rotation oblique
    if (!is.null(results$phi) && is.matrix(results$phi)) {
      cat("\nMATRICE DE CORRÉLATION ENTRE FACTEURS (PHI)\n")
      cat("-------------------------------------------\n")
      
      # En-tête du tableau
      cat(sprintf("%-10s", ""))
      for (i in 1:ncol(results$phi)) {
        cat(sprintf(" %8s", paste0("F", i)))
      }
      cat("\n")
      
      # Contenu du tableau
      for (i in 1:nrow(results$phi)) {
        cat(sprintf("%-10s", paste0("F", i)))
        for (j in 1:ncol(results$phi)) {
          val <- results$phi[i, j]
          cat(sprintf(" %8.3f", val))
        }
        cat("\n")
      }
      cat("\n")
    }

    # Afficher la matrice de corrélation reproduite
    if (!is.null(results$cormat_reprod) && is.matrix(results$cormat_reprod)) {
      cat("MATRICE DE CORRÉLATION REPRODUITE\n")
      cat("---------------------------------\n")
      
      # Sélectionner un sous-ensemble de variables si la matrice est grande
      cormat_reprod <- results$cormat_reprod
      max_display <- 10  # Maximum de variables à afficher pour la lisibilité
      
      if (nrow(cormat_reprod) > max_display) {
        cat(sprintf("Note: Affichage des %d premières variables (sur %d au total)\n\n", 
                    max_display, nrow(cormat_reprod)))
        
        cormat_sub <- cormat_reprod[1:max_display, 1:max_display]
      } else {
        cormat_sub <- cormat_reprod
      }
      
      # En-tête du tableau
      cat(sprintf("%-10s", ""))
      for (i in 1:ncol(cormat_sub)) {
        var_name <- colnames(cormat_sub)[i]
        if (is.null(var_name)) var_name <- paste0("Var", i)
        # Prendre les 8 premiers caractères du nom de variable pour l'affichage
        short_name <- substr(var_name, 1, 8)
        cat(sprintf(" %8s", short_name))
      }
      cat("\n")
      
      # Contenu du tableau
      for (i in 1:nrow(cormat_sub)) {
        var_name <- rownames(cormat_sub)[i]
        if (is.null(var_name)) var_name <- paste0("Var", i)
        # Prendre les 10 premiers caractères pour l'affichage
        short_name <- substr(var_name, 1, 10)
        cat(sprintf("%-10s", short_name))
        
        for (j in 1:ncol(cormat_sub)) {
          cat(sprintf(" %8.3f", cormat_sub[i, j]))
        }
        cat("\n")
      }
      
      # Note sur l'interprétation
      cat("\nNote: La matrice de corrélation reproduite montre les corrélations reconstruites à partir du modèle factoriel.\n")
      cat("Les différences entre cette matrice et la matrice de corrélation originale sont les résidus.\n\n")
    } 
    
    # Indices d'ajustement du modèle
    if (!is.null(results$fit_coefs)) {
      cat("CRITÈRES D'AJUSTEMENT DÉTAILLÉS\n")
      cat("------------------------------\n")
      
      # Créer un tableau organisé des indices d'ajustement
      if (!is.null(results$fit_coefs$RMSR)) cat(sprintf("%-10s: %7.4f\n", "RMSR", results$fit_coefs$RMSR))
      if (!is.null(results$fit_coefs$GFI)) cat(sprintf("%-10s: %7.4f\n", "GFI", results$fit_coefs$GFI))
      if (!is.null(results$fit_coefs$CAF)) cat(sprintf("%-10s: %7.4f\n", "CAF", results$fit_coefs$CAF))
      if (!is.null(results$fit_coefs$TLI)) cat(sprintf("%-10s: %7.4f\n", "TLI", results$fit_coefs$TLI))
      if (!is.null(results$fit_coefs$CFI)) cat(sprintf("%-10s: %7.4f\n", "CFI", results$fit_coefs$CFI))
      if (!is.null(results$fit_coefs$MFI)) cat(sprintf("%-10s: %7.4f\n", "MFI", results$fit_coefs$MFI))
      
      # Critères d'information
      cat("\nCritères d'information:\n")
      if (!is.null(results$fit_coefs$AIC)) cat(sprintf("%-10s: %7.2f\n", "AIC", results$fit_coefs$AIC))
      if (!is.null(results$fit_coefs$BIC)) cat(sprintf("%-10s: %7.2f\n", "BIC", results$fit_coefs$BIC))
      if (!is.null(results$fit_coefs$CAIC)) cat(sprintf("%-10s: %7.2f\n", "CAIC", results$fit_coefs$CAIC))
      if (!is.null(results$fit_coefs$SABIC)) cat(sprintf("%-10s: %7.2f\n", "SABIC", results$fit_coefs$SABIC))
      
      cat("\n")
    }
    
    # Communalités et unicités
    if (!is.null(results$communalities) || !is.null(results$uniquenesses)) {
      cat("COMMUNALITÉS ET UNICITÉS\n")
      cat(paste(rep("-", 45), collapse = ""), "\n")
      
      # En-tête du tableau
      cat(sprintf("%-15s  %12s  %12s\n", "Variable", "Communalité", "Unicité"))
      cat(paste(rep("-", 45), collapse = ""), "\n")
      
      # Afficher communalités et unicités pour chaque variable
      if (!is.null(results$communalities) && is.matrix(results$communalities)) {
        n_vars <- nrow(results$communalities)
        
        for (i in 1:n_vars) {
          var_name <- rownames(results$communalities)[i]
          if (is.null(var_name)) var_name <- paste0("Var", i)
          
          # Récupérer les valeurs d'extraction
          comm_value <- results$communalities[i, "Extraction"]
          unique_value <- results$uniquenesses[i, "Extraction"]
          
          # Formatage avec interprétation visuelle de la communalité
          if (comm_value >= 0.5) {
            quality <- "Élevée"
          } else if (comm_value >= 0.3) {
            quality <- "Acceptable"
          } else {
            quality <- "Faible"
          }
          cat(sprintf("%-15s  %10.4f  %10.4f  (%s)\n", var_name, comm_value, unique_value, quality))
        }
      }
      cat("\n")
    }
    
    # Afficher la variance expliquée après rotation (varexplROT)
    if (!is.null(results$varexplROT) && is.matrix(results$varexplROT)) {
      cat("VARIANCE EXPLIQUÉE APRÈS ROTATION\n")
      cat("---------------------------------\n")
      
      varexpl <- results$varexplROT
      
      # En-tête du tableau
      cat(sprintf("%-10s  %14s  %14s  %14s\n", 
                  "Facteur", "Somme carrés", "% de variance", "% cumulatif"))
      cat(paste(rep("-", 60), collapse = ""), "\n")
      
      # Contenu du tableau
      for (i in 1:nrow(varexpl)) {
        cat(sprintf("%-10s  %14.4f  %14.4f%%  %14.4f%%\n", 
                    paste0("Facteur ", i),
                    varexpl[i, 1],  # Sommes des carrés
                    varexpl[i, 2] * 100,  # Proportion de variance
                    varexpl[i, 3] * 100   # Proportion cumulative
        ))
      }
      cat("\n")
    }
    
    # Conclusion
    cat("CONCLUSION ET INTERPRÉTATION\n")
    cat("----------------------------\n")
    
    # Recommandation générale sur le nombre de facteurs
    if (input$factor_method == "auto") {
      cat(sprintf("L'analyse factorielle exploratoire a identifié %d facteur(s) en utilisant la méthode automatique.\n", results$nfactors))
    } else {
      cat(sprintf("L'analyse factorielle a été réalisée avec %d facteur(s) spécifié(s) par l'utilisateur.\n", input$num_factors))
    }
    
    # Synthèse de la variance expliquée
    # Pour la variance totale expliquée
    if (!is.null(results$varexplROT)) {
      if (is.matrix(results$varexplROT)) {
        total_var_explained <- results$varexplROT[nrow(results$varexplROT), 3] * 100
      } else if (is.vector(results$varexplROT)) {
        total_var_explained <- sum(results$varexplROT) * 100
      }
    } else if (!is.null(results$varexplNOROT2) && is.matrix(results$varexplNOROT2)) {
      total_var_explained <- results$varexplNOROT2[nrow(results$varexplNOROT2), 3] * 100
    } else if (!is.null(results$eigenvalues)) {
      eigenvalues <- results$eigenvalues[1:nfactors]
      total_var_explained <- sum(eigenvalues) / length(results$eigenvalues) * 100
    } else {
      cat("Information sur la variance expliquée non disponible.\n")
    }

    # Commentaire sur la qualité d'ajustement du modèle
    if (!is.null(results$rmsea)) {
      cat("\nQualité d'ajustement du modèle: ")
      if (results$rmsea < 0.05) {
        cat("Excellente (RMSEA < 0.05)\n")
      } else if (results$rmsea < 0.08) {
        cat("Bonne (RMSEA < 0.08)\n")
      } else if (results$rmsea < 0.10) {
        cat("Médiocre (RMSEA < 0.10)\n")
      } else {
        cat("Insuffisante (RMSEA ≥ 0.10)\n")
      }
    }
    
    # Recommandations pour l'interprétation
    cat("\nRECOMMANDATIONS POUR L'INTERPRÉTATION:\n")
    cat("1. Examinez les saturations factorielles supérieures à 0.3 en valeur absolue (affichées ci-dessus).\n")
    cat("2. Recherchez les variables qui saturent fortement sur un seul facteur pour l'interprétation thématique.\n")
    cat("3. Les variables à saturation faible (< 0.3) sur tous les facteurs peuvent être considérées pour élimination.\n")
    cat("4. Examinez les communalités pour identifier les variables mal représentées par la solution factorielle.\n")
    
    if (input$rotation_type == "oblique") {
      cat("5. Tenez compte des corrélations entre facteurs dans votre interprétation théorique.\n")
    }
    
    if (!is.null(results$h2)) {
      # Identifier les variables avec communalités faibles
      low_h2_vars <- which(results$h2 < 0.3)
      if (length(low_h2_vars) > 0) {
        cat("\nVariables avec communalités faibles (< 0.3):\n")
        for (i in low_h2_vars) {
          var_name <- rownames(loadings_matrix)[i]
          if (is.null(var_name)) var_name <- paste0("Var", i)
          cat(sprintf("- %s (h² = %.3f): Considérer pour élimination\n", var_name, results$h2[i]))
        }
      }
    }
    
    # Identifier les variables à saturations complexes
    if (!is.null(loadings_matrix)) {
      complex_vars <- c()
      for (i in 1:nrow(loadings_matrix)) {
        sig_loadings <- sum(abs(loadings_matrix[i,]) >= 0.3)
        if (sig_loadings > 1) {
          complex_vars <- c(complex_vars, i)
        }
      }
      
      if (length(complex_vars) > 0) {
        cat("\nVariables à saturations complexes (saturant sur plusieurs facteurs):\n")
        for (i in complex_vars) {
          var_name <- rownames(loadings_matrix)[i]
          if (is.null(var_name)) var_name <- paste0("Var", i)
          cat(sprintf("- %s: Saturations significatives sur %d facteurs\n", 
                      var_name, sum(abs(loadings_matrix[i,]) >= 0.3)))
        }
        cat("Ces variables pourraient être ambiguës dans l'interprétation des facteurs.\n")
      }
    }
    
    # Conseils finaux
    cat("\nRAPPEL: \n")
    cat("- L'interprétation des facteurs est subjective et doit être guidée par la théorie.\n")
    cat("- La fiabilité de la solution factorielle dépend de la taille de l'échantillon et de la qualité des données.\n")
    cat("- Une analyse confirmatoire peut être nécessaire pour valider la structure factorielle identifiée.\n")
  })
  
  # Téléchargement des données simulées
  output$download_data <- downloadHandler(
    filename = function() {
      paste("donnees_simulees_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
    },
    content = function(file) {
      req(sim_data())
      write.csv(sim_data()$response, file, row.names = FALSE)
    }
  )
  
  # UI dynamique pour sélectionner les facteurs pour le biplot
  output$factor_selector_ui <- renderUI({
    req(efa_result())
    results <- efa_result()
    
    # Déterminer le nombre de facteurs
    n_factors <- 0
    if (!is.null(results$nfactors)) {
      n_factors <- results$nfactors
    } else if (!is.null(results$loadings) && is.matrix(results$loadings)) {
      n_factors <- ncol(results$loadings)
    } else if (!is.null(results$loadingsROT) && is.matrix(results$loadingsROT)) {
      n_factors <- ncol(results$loadingsROT)
    } else if (!is.null(results$pattern) && is.matrix(results$pattern)) {
      n_factors <- ncol(results$pattern)
    }
    
    # Si moins de 2 facteurs, afficher un message
    if (n_factors < 2) {
      return(div(
        class = "alert alert-warning",
        "Le biplot nécessite au moins 2 facteurs. La solution actuelle ne contient qu'un seul facteur."
      ))
    }
    
    # Si exactement 2 facteurs, pas besoin de sélection
    if (n_factors == 2) {
      return(NULL)
    }
    
    # Sinon, créer des contrôles pour sélectionner les facteurs
    fluidRow(
      column(6,
             selectInput("biplot_x_factor", "Facteur pour l'axe X:",
                         choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
                         selected = 1)
      ),
      column(6,
             selectInput("biplot_y_factor", "Facteur pour l'axe Y:",
                         choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
                         selected = 2)
      )
    )
  })
  
  # Biplot interactif avec plotly
  output$efa_biplot <- renderPlotly({
    req(efa_result())
    results <- efa_result()
    
    # Déterminer quelle matrice de loadings utiliser
    loadings_matrix <- NULL
    is_none <- input$rotation_type == "none"
    is_orthogonal <- input$rotation_type == "orthogonal"

    if (is_none) {loadings_matrix <- results$loadingsNOROT 
    } else if (is_orthogonal){loadings_matrix <- results$loadingsROT
    } else {loadings_matrix <- results$pattern}

    # Si aucune matrice de loadings n'est disponible
    if (is.null(loadings_matrix)) {
      return(plotly_empty() %>% 
               layout(title = "Données insuffisantes pour générer un biplot"))
    }
    
    # Déterminer le nombre de facteurs
    n_factors <- ncol(loadings_matrix)
    
    # Si un seul facteur, ne pas générer de biplot
    if (n_factors < 2) {
      return(plotly_empty() %>% 
               layout(title = "Le biplot nécessite au moins 2 facteurs"))
    }
    
    # Déterminer quels facteurs utiliser pour le biplot
    factor_x <- 1
    factor_y <- 2
    
    if (n_factors > 2) {
      factor_x <- as.numeric(input$biplot_x_factor)
      factor_y <- as.numeric(input$biplot_y_factor)
    }
    
    # Récupérer les pourcentages de variance expliquée
    variance_explained <- NULL
    if (!is.null(results$varexplROT) && is.matrix(results$varexplROT)) {
      variance_explained <- results$varexplROT[, 2] * 100
    } else if (!is.null(results$eigenvalues)) {
      variance_explained <- results$eigenvalues / sum(results$eigenvalues) * 100
    }
    
    # Créer un dataframe pour le biplot
    biplot_data <- data.frame(
      Variable = rownames(loadings_matrix),
      X = round(loadings_matrix[, factor_x], 4),
      Y = round(loadings_matrix[, factor_y], 4)
    )
    
    # Labels des axes
    x_label <- paste0("Facteur ", factor_x)
    y_label <- paste0("Facteur ", factor_y)
    
    if (!is.null(variance_explained)) {
      x_label <- paste0(x_label, " (", round(variance_explained[factor_x], 1), "%)")
      y_label <- paste0(y_label, " (", round(variance_explained[factor_y], 1), "%)")
    }
    
    # Calculer les limites du graphique pour avoir un aspect carré
    max_abs_value <- max(abs(biplot_data$X), abs(biplot_data$Y)) * 1.1
    
    # Créer le biplot interactif
    p <- plot_ly() %>%
      add_trace(
        data = biplot_data,
        x = ~X, 
        y = ~Y,
        type = "scatter",
        mode = "markers+text",
        text = ~Variable,
        textposition = "top center",
        marker = list(size = 10, color = "blue"),
        textfont = list(size = 12),
        hoverinfo = "text",
        hovertext = ~paste("Variable:", Variable, 
                           "<br>", x_label, ":", round(X, 3),
                           "<br>", y_label, ":", round(Y, 3))
      ) %>%
      add_segments(
        x = 0, xend = 0,
        y = -max_abs_value, yend = max_abs_value,
        line = list(color = "gray", width = 1, dash = "dash"),
        showlegend = FALSE
      ) %>%
      add_segments(
        x = -max_abs_value, xend = max_abs_value,
        y = 0, yend = 0,
        line = list(color = "gray", width = 1, dash = "dash"),
        showlegend = FALSE
      ) %>%
      layout(
        title = list(
          text = paste("Biplot des Facteurs", factor_x, "et", factor_y),
          font = list(size = 16)
        ),
        xaxis = list(
          title = x_label,
          zeroline = FALSE,
          range = c(-max_abs_value, max_abs_value)
        ),
        yaxis = list(
          title = y_label,
          zeroline = FALSE,
          range = c(-max_abs_value, max_abs_value),
          scaleanchor = "x",  # Pour maintenir un aspect ratio égal
          scaleratio = 1
        ),
        showlegend = FALSE,
        annotations = list(
          list(
            x = 0,
            y = max_abs_value,
            xref = "x",
            yref = "y",
            text = "Variables et leurs charges factorielles",
            showarrow = FALSE,
            font = list(size = 14)
          )
        )
      )
    
    return(p)
  })
  
  output$complexity_results <- renderPrint({
    req(efa_result())
    results <- efa_result()
    
    # Déterminer quelle matrice de loadings utiliser
    loadings_matrix <- NULL
    
    if (!is.null(results$pattern)) {
      loadings_matrix <- results$pattern
    } else if (!is.null(results$loadingsROT)) {
      loadings_matrix <- results$loadingsROT
    } else if (!is.null(results$loadings)) {
      loadings_matrix <- results$loadings
    }
    
    # Si une matrice de loadings est disponible, calculer la complexité
    if (!is.null(loadings_matrix)) {
      cat("ANALYSE DE LA COMPLEXITÉ DE LA SOLUTION FACTORIELLE\n")
      cat("=================================================\n\n")
      
      # Vérifier que loadings_matrix est bien numérique
      if (!is.numeric(loadings_matrix)) {
        cat("Erreur: La matrice des charges factorielles n'est pas numérique.\n")
        return()
      }
      
      # Vérifier s'il y a des valeurs NA dans la matrice
      if (any(is.na(loadings_matrix))) {
        cat("Attention: La matrice contient des valeurs manquantes (NA).\n")
        loadings_matrix[is.na(loadings_matrix)] <- 0  # Remplacer les NA par 0
      }
      
      tryCatch({
        # Vérifier qu'il y a au moins 2 facteurs pour calculer la complexité
        if (ncol(loadings_matrix) < 2) {
          cat("Impossible de calculer la complexité: au moins 2 facteurs sont nécessaires.\n")
          return()
        }
        
        # Calculer la complexité manuellement
        n_vars <- nrow(loadings_matrix)
        complexity_values <- numeric(n_vars)
        
        for (i in 1:n_vars) {
          loadings_squared <- loadings_matrix[i,]^2
          complexity_values[i] <- sum(loadings_squared)^2 / sum(loadings_squared^2)
        }
        
        names(complexity_values) <- rownames(loadings_matrix)
        
        # Afficher les valeurs de complexité pour chaque variable
        cat("VALEURS DE COMPLEXITÉ PAR VARIABLE\n")
        cat("--------------------------------\n")
        # En-tête du tableau avec ajout de la colonne de pourcentage
        cat(sprintf("%-20s  %10s  %10s\n", "Variable", "Complexité", "Pourcentage"))
        cat(paste(rep("-", 44), collapse = ""), "\n")
        
        # Afficher les valeurs pour chaque variable
        n_factors <- ncol(loadings_matrix)  # Nombre de facteurs
        
        for (i in 1:length(complexity_values)) {
          var_name <- names(complexity_values)[i]
          if (is.null(var_name) || var_name == "") var_name <- paste0("Var", i)
          
          # Calculer le pourcentage de complexité pour chaque variable
          # La formule est: (complexité - 1) / (nombre de facteurs - 1) * 100
          pct_complexity <- (complexity_values[i] - 1) / (n_factors - 1) * 100
          
          cat(sprintf("%-20s  %10.3f  %10.1f%%\n", var_name, complexity_values[i], pct_complexity))
        }
        cat("\n")
        
        # Calculer et afficher les statistiques de résumé
        mean_c <- mean(complexity_values)
        sd_c <- sd(complexity_values)
        min_c <- min(complexity_values)
        max_c <- max(complexity_values)
        
        # Afficher un résumé
        cat("RÉSUMÉ DE LA COMPLEXITÉ\n")
        cat("----------------------\n")
        cat("Complexité moyenne: ", round(mean_c, 3), "\n")
        cat("Écart-type de la complexité: ", round(sd_c, 3), "\n")
        cat("Complexité minimale: ", round(min_c, 3), "\n")
        cat("Complexité maximale: ", round(max_c, 3), "\n\n")
        
        # Afficher le pourcentage de complexité
        percent_complexity <- (mean_c - 1) / (ncol(loadings_matrix) - 1) * 100
        cat("Pourcentage de complexité: ", round(percent_complexity, 1), "%\n\n")
        
        # Interprétation
        cat("INTERPRÉTATION\n")
        cat("-------------\n")
        if (mean_c < 1.2) {
          cat("Très bonne solution: La plupart des variables chargent sur un seul facteur.\n")
        } else if (mean_c < 1.5) {
          cat("Bonne solution: La structure factorielle est claire.\n")
        } else if (mean_c < 2) {
          cat("Solution acceptable: Quelques variables ont des charges croisées.\n")
        } else {
          cat("Solution complexe: Nombreuses variables avec des charges croisées.\n")
        }
        
        # Explication de la complexité
        cat("\nEXPLICATION DE LA COMPLEXITÉ\n")
        cat("-------------------------\n")
        cat("- Une complexité de 1 indique qu'une variable charge sur un seul facteur\n")
        cat("- Une complexité de", ncol(loadings_matrix), "indiquerait des charges égales sur tous les facteurs\n")
        cat("- Plus la complexité est élevée, plus la variable est ambiguë dans la structure factorielle\n")
        cat("\n")
        cat("Le pourcentage de complexité par variable indique à quel point la solution factorielle est \n ")
        cat("complexe. Une solution où toutes les variables sont fortement associées à un seul facteur \n")
        cat("donnera un pourcentage de complexité faible, tandis qu'une solution où plusieurs variables\n")
        cat("sont chargées sur plusieurs facteurs aura une complexité plus élevée.\n\n")
        cat("Interprétation des pourcentages :\n")
        cat("- 0-20% : Excellente simplicité factorielle\n")
        cat("- 20-40% : Bonne simplicité factorielle\n")
        cat("- 40-60% : Solution factorielle moyennement complexe\n")
        cat("- 60-80% : Solution factorielle complexe\n")
        cat("- 80-100% : Solution factorielle très complexe (ambiguïté élevée)\n")
        
      }, error = function(e) {
        cat("Erreur lors du calcul de la complexité:", e$message, "\n")
        cat("\nDétails techniques pour le débogage:\n")
        cat("Dimensions de la matrice:", paste(dim(loadings_matrix), collapse=" x "), "\n")
        cat("Classe de la matrice:", paste(class(loadings_matrix), collapse=", "), "\n")
        cat("Type des éléments:", typeof(loadings_matrix), "\n")
        
        # Afficher un échantillon des données pour vérifier
        cat("\nÉchantillon de la matrice des charges (5 premières lignes si disponibles):\n")
        if (nrow(loadings_matrix) > 0) {
          sample_rows <- min(5, nrow(loadings_matrix))
          print(loadings_matrix[1:sample_rows, , drop=FALSE])
        }
      })
    } else {
      cat("Aucune matrice de loadings disponible pour calculer la complexité.")
    }
  })
  
  # Observer pour le calcul des scores factoriels
  observeEvent(input$run_factor_scores, {
    req(sim_data(), efa_result())
    
    # Afficher un indicateur de chargement
    withProgress(message = 'Calcul des scores factoriels...', {
      # Récupérer les données
      data_obj <- sim_data()
      response_data <- data_obj$response
      
      # Récupérer les résultats de l'analyse factorielle
      results <- efa_result()
      
      # Déterminer quelle matrice de loadings utiliser
      loadings_matrix <- NULL
      loadings_type <- "structure"  # Par défaut
      phi_matrix <- NULL
 
      if(input$rotation_type == "none"){loadings_matrix <- results$loadingsNOROT}   
      if(input$rotation_type == "orthogonal"){loadings_matrix <- results$loadingsROT}
      if(input$rotation_type == "oblique"){
        loadings_matrix <- results$pattern
        phi_matrix <- results$phi
      }

      # Méthode de calcul des scores
      method <- input$factor_scores_method
      
      # Vérification des conditions nécessaires
      if (is.null(loadings_matrix)) {
        showNotification("Erreur: Aucune matrice de charges factorielles disponible.", 
                         type = "error", duration = 10)
        return()
      }
      
      # Exécuter le calcul des scores avec gestion d'erreur
      tryCatch({
        # Utiliser EFA_SCORES de EFA.dimensions
        scores_result <- EFA.dimensions::EFA_SCORES(
          loadings = loadings_matrix,
          loadings_type = loadings_type,
          data = response_data,
          corkind = input$efa_cor,
          phi = phi_matrix,
          method = method,
          verbose = FALSE
        )
        
        # Stocker les résultats
        factor_scores_result(scores_result)
        
        # Notification
        showNotification("Calcul des scores factoriels complété !", type = "message")
        
      }, error = function(e) {
        # En cas d'erreur
        showNotification(paste("Erreur lors du calcul des scores factoriels:", e$message), 
                         type = "error", duration = 10)
      })
    })
  })
  
  # Output pour les statistiques des scores factoriels
  output$factor_scores_stats <- renderPrint({
    req(factor_scores_result())
    results <- factor_scores_result()
    
    cat("RÉSULTATS DU CALCUL DES SCORES FACTORIELS\n")
    cat("========================================\n\n")
    
    cat("MÉTHODE UTILISÉE\n")
    cat("--------------\n")
    cat("Méthode:", input$factor_scores_method, "\n\n")
    
    # Coefficients des scores factoriels
    if (!is.null(results$FSCoef)) {
      cat("COEFFICIENTS DES SCORES FACTORIELS\n")
      cat("--------------------------------\n")
      # En-tête du tableau
      cat(sprintf("%-15s", "Variable"))
      for (i in 1:ncol(results$FSCoef)) {
        cat(sprintf("  %8s", paste0("F", i)))
      }
      cat("\n")
      
      # Séparateur
      cat(paste(rep("-", 15 + 10*ncol(results$FSCoef)), collapse = ""), "\n")
      
      # Valeurs
      for (i in 1:nrow(results$FSCoef)) {
        var_name <- rownames(results$FSCoef)[i]
        if (is.null(var_name)) var_name <- paste0("Var", i)
        
        cat(sprintf("%-15s", var_name))
        for (j in 1:ncol(results$FSCoef)) {
          cat(sprintf("  %8.3f", results$FSCoef[i, j]))
        }
        cat("\n")
      }
      cat("\n")
    }
    
    # Indicateurs de qualité
    cat("INDICATEURS DE QUALITÉ DES SCORES\n")
    cat("------------------------------\n")
    
    # Multiple R
    if (!is.null(results$MULTR)) {
      cat("\nCorrélations multiples (MULTR) - Qualité de prédiction:\n")
      for (i in 1:length(results$MULTR)) {
        cat(sprintf("Facteur %d: %.3f\n", i, results$MULTR[i]))
      }
    }
    
    # R-carré
    if (!is.null(results$RSQR)) {
      cat("\nCoefficients de détermination (RSQR) - Variance expliquée:\n")
      for (i in 1:length(results$RSQR)) {
        cat(sprintf("Facteur %d: %.3f\n", i, results$RSQR[i]))
      }
    }
    
    # MINCOR (indice d'indétermination de Guttman)
    if (!is.null(results$MINCOR)) {
      cat("\nIndice d'indétermination de Guttman (MINCOR):\n")
      for (i in 1:length(results$MINCOR)) {
        cat(sprintf("Facteur %d: %.3f\n", i, results$MINCOR[i]))
      }
    }
    
    # Validité
    if (!is.null(results$VALIDITY)) {
      cat("\nValidité (corrélation entre scores estimés et facteurs):\n")
      for (i in 1:length(results$VALIDITY)) {
        cat(sprintf("Facteur %d: %.3f\n", i, results$VALIDITY[i]))
      }
    }
    
    # Univocité
    if (!is.null(results$UNIVOCALITY)) {
      cat("\nUnivocité (pureté des scores factoriels):\n")
      for (i in 1:length(results$UNIVOCALITY)) {
        cat(sprintf("Facteur %d: %.3f\n", i, results$UNIVOCALITY[i]))
      }
    }
    
    # Afficher les corrélations entre scores factoriels
    if (!is.null(results$FactorScore_Correls)) {
      cat("\nCORRÉLATIONS ENTRE LES SCORES FACTORIELS\n")
      cat("-------------------------------------\n")
      
      # En-tête du tableau
      cat(sprintf("%-10s", ""))
      for (i in 1:ncol(results$FactorScore_Correls)) {
        cat(sprintf("  %8s", paste0("F", i)))
      }
      cat("\n")
      
      # Séparateur
      cat(paste(rep("-", 10 + 10*ncol(results$FactorScore_Correls)), collapse = ""), "\n")
      
      # Valeurs
      for (i in 1:nrow(results$FactorScore_Correls)) {
        cat(sprintf("%-10s", paste0("F", i)))
        for (j in 1:ncol(results$FactorScore_Correls)) {
          val <- results$FactorScore_Correls[i, j]
          if (i == j) {
            # Diagonale (toujours 1)
            cat(sprintf("  %8s", "1.000"))
          } else {
            cat(sprintf("  %8.3f", val))
          }
        }
        cat("\n")
      }
      cat("\n")
    }
  })

# Affichage des scores factoriels individuels dans une table
  output$factor_scores_table <- renderDT({
    req(factor_scores_result())
    
    results <- factor_scores_result()
    
    # Créer un dataframe pour la table
    scores_df <- as.data.frame(results$FactorScores)
    
    # Renommer les colonnes pour plus de clarté
    colnames(scores_df) <- paste0("Facteur ", 1:ncol(scores_df))
    
    # Ajouter une colonne d'ID
    scores_df <- cbind(ID = 1:nrow(scores_df), scores_df)
    
    # Afficher avec datatable
    datatable(
      scores_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = 'tip'
      ),
      caption = "Scores factoriels individuels"
    ) %>%
      formatRound(columns = 2:(ncol(scores_df)), digits = 3)
  })
  
  # La correction complète serait donc:
  output$dimtests_results <- renderPrint({
    req(dimtests_results())
    results <- dimtests_results()
    
    cat("RÉSULTATS DÉTAILLÉS DES TESTS DIMENSIONNELS\n")
    cat("==========================================\n\n")
    
    # Afficher les résultats de chaque test
    if (!is.null(results$final)) {
      cat("RÉSUMÉ DES TESTS\n")
      cat("--------------\n")
      
      # Récupérer le nombre réel de facteurs
      n_factors_real <- attr(sim_data(), "n.factors")
      
      # Parcourir tous les tests
      for (test_name in names(results$final)) {
        if (test_name != "NOTE") {
          test_value <- results$final[[test_name]]
          
          # Vérifier que test_value n'est ni NULL ni NA
          if (!is.null(test_value) && !is.na(test_value)) {
            cat(sprintf("%-15s: %d facteur(s)", test_name, test_value))
            
            # Comparer avec le nombre réel de facteurs
            if (!is.null(n_factors_real) && !is.na(n_factors_real)) {
              if (test_value == n_factors_real) {
                cat(" ✓ (correspond au nombre réel)\n")
              } else {
                cat(" ✗ (diffère du nombre réel)\n")
              }
            } else {
              cat("\n")
            }
          } else {
            cat(sprintf("%-15s: Non disponible\n", test_name))
          }
        }
      }
      cat("\n")
    }
    
    # Afficher les valeurs propres
    if (!is.null(results$eigenvalues)) {
      cat("\n==============================================\n")
      cat("| VALEURS PROPRES (EIGENVALUES)              |\n")
      cat("==============================================\n\n")
      
      n_eigenvalues <- length(results$eigenvalues)
      cumulative <- cumsum(results$eigenvalues) / sum(results$eigenvalues) * 100
      
      cat(sprintf("%-10s  %12s  %12s  %12s\n", 
                  "Facteur", "Valeur propre", "% de variance", "% cumulatif"))
      cat(paste(rep("-", 50), collapse = ""), "\n")
      
      for (i in 1:n_eigenvalues) {
        cat(sprintf("%-10d  %12.4f  %12.2f%%  %12.2f%%\n", 
                    i, results$eigenvalues[i], 
                    results$eigenvalues[i] / sum(results$eigenvalues) * 100,
                    cumulative[i]))
      }
      cat("\n")
    }
    
    # Afficher les résultats détaillés du test parallèle si disponibles
    if (!is.null(results$parallel_data)) {
      cat("\n==============================================\n")
      cat("| TEST PARALLÈLE (PARALLEL ANALYSIS)         |\n")
      cat("==============================================\n\n")
      
      # Récupérer les données du test parallèle
      observed <- NULL
      if (!is.null(results$eigenvalues)) {
        observed <- results$eigenvalues
      } else if (!is.null(results$parallel_data$observed)) {
        observed <- results$parallel_data$observed
      }
      
      # Utiliser les noms corrects pour les champs Mean et Percentile
      mean_evals <- results$parallel_data$Mean
      percentile_evals <- results$parallel_data$Percentile
      
      if (!is.null(observed) && !is.null(mean_evals) && !is.null(percentile_evals)) {
        # Limiter aux 10 premiers facteurs pour la lisibilité
        max_display <- min(10, length(observed))
        
        cat(sprintf("%-10s  %15s  %15s  %15s\n", 
                    "Facteur", "Valeur observée", "Moyenne aléatoire", 
                    paste0(input$rawpar_percentile, "e percentile")))
        cat(paste(rep("-", 60), collapse = ""), "\n")
        
        # Afficher les valeurs pour chaque facteur jusqu'à max_display
        for (i in 1:max_display) {
          cat(sprintf("%-10d  %15.4f  %15.4f  %15.4f\n", 
                      i, 
                      observed[i], 
                      mean_evals[i], 
                      percentile_evals[i]))
        }
        cat("\n")
        
        # Afficher le résultat de la comparaison
        n_factors_parallel <- sum(observed > percentile_evals)
        cat(sprintf("Nombre de facteurs suggéré par le test parallèle: %d\n", n_factors_parallel))
        cat("(Basé sur le nombre de valeurs propres observées > valeurs propres aléatoires au percentile spécifié)\n\n")
        
        # Explications additionnelles
        cat("INTERPRÉTATION:\n")
        cat("- Les valeurs propres observées (eigenvalues) représentent la variance expliquée par chaque facteur\n")
        cat("- Les valeurs propres aléatoires sont générées à partir de données aléatoires ayant la même taille\n")
        cat("- Le critère de rétention: conserver les facteurs dont la valeur propre observée\n")
        cat("  dépasse la valeur correspondante au percentile spécifié des eigenvalues aléatoires\n")
      } else {
        cat("Données insuffisantes pour afficher les résultats détaillés du test parallèle.\n")
      }
    }
    
    # Test CD (Coordonnées optimales)
    if (!is.null(results$CD_details) || "CD" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST CD (COORDONNÉES OPTIMALES)            |\n")
      cat("==============================================\n\n")
      
      if (!is.null(results$CD_details)) {
        cat(paste(results$CD_details, collapse = "\n"), "\n\n")
      } else if ("CD" %in% names(results$final)) {
        cat("Description: Ce test utilise les coordonnées optimales pour déterminer le nombre de dimensions.\n")
        cat("Il évalue la reproductibilité des corrélations observées à partir d'un nombre donné de dimensions.\n")
        cat(sprintf("Nombre de facteurs suggéré: %d\n\n", results$final[["CD"]]))
      }
    }
    
    # Test EMPKC (Critère d'information empirique de Kaiser)
    if (!is.null(results$EMPKC_details) || "EMPKC" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST EMPKC (CRITÈRE DE KAISER EMPIRIQUE)   |\n")
      cat("==============================================\n\n")
      
      if (!is.null(results$EMPKC_details)) {
        cat(paste(results$EMPKC_details, collapse = "\n"), "\n\n")
      } else if ("EMPKC" %in% names(results$final)) {
        cat("Description: Extension du critère de Kaiser qui considère non seulement les eigenvalues > 1\n")
        cat("mais aussi la proportion de variance expliquée par chaque facteur.\n")
        cat(sprintf("Nombre de facteurs suggéré: %d\n\n", results$final[["EMPKC"]]))
      }
    }
    
    # Test HULL
    if (!is.null(results$HULL_details) || "HULL" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST HULL (MÉTHODE HULL)                   |\n")
      cat("==============================================\n\n")
      
      if (!is.null(results$HULL_details)) {
        cat(paste(results$HULL_details, collapse = "\n"), "\n\n")
      } else if ("HULL" %in% names(results$final)) {
        cat("Description: Ce test recherche le point optimal sur l'enveloppe convexe (hull) des points\n")
        cat("représentant les degrés de liberté vs l'ajustement pour différents nombres de facteurs.\n")
        cat(sprintf("Méthode d'extraction utilisée: %s\n", input$dimtests_hull_method))
        cat(sprintf("Indice d'ajustement utilisé: %s\n", input$dimtests_hull_gof))
        cat(sprintf("Nombre de facteurs suggéré: %d\n\n", results$final[["HULL"]]))
      }
    }
    
    # Test MAP (Minimum Average Partial)
    if (!is.null(results$MAP_details) || "MAP" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST MAP (MINIMUM AVERAGE PARTIAL)         |\n")
      cat("==============================================\n\n")
      
      if (!is.null(results$MAP_details)) {
        cat(paste(results$MAP_details, collapse = "\n"), "\n\n")
      } else if ("MAP" %in% names(results$final)) {
        cat("Description: Ce test examine les corrélations partielles après élimination successive des facteurs.\n")
        cat("Le nombre optimal de facteurs est déterminé par le minimum de la moyenne des corrélations partielles au carré.\n")
        cat(sprintf("Nombre de facteurs suggéré: %d\n\n", results$final[["MAP"]]))
      }
    }
    
    # Test NEVALSGT1 (Eigenvalues > 1)
    if ("NEVALSGT1" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST NEVALSGT1 (EIGENVALUES > 1)           |\n")
      cat("==============================================\n\n")
      
      cat("Description: Critère de Kaiser classique qui retient les facteurs dont la valeur propre est supérieure à 1.\n")
      
      if (!is.null(results$eigenvalues)) {
        cat("Valeurs propres observées:\n")
        for (i in 1:min(10, length(results$eigenvalues))) {
          cat(sprintf("Facteur %d: %.4f %s\n", 
                      i, 
                      results$eigenvalues[i],
                      ifelse(results$eigenvalues[i] > 1, "(Retenu)", "(Non retenu)")))
        }
      }
      
      cat(sprintf("\nNombre de facteurs suggéré: %d\n\n", results$final[["NEVALSGT1"]]))
    }
    
    # Test SESCREE (Scree test)
    if (!is.null(results$SESCREE_details)) {
      cat("\n==============================================\n")
      cat("| TEST SESCREE (SCREE TEST)                |\n")
      cat("==============================================\n\n")
      
      cat(paste(results$SESCREE_details, collapse = "\n"), "\n\n")
    } else if ("SESCREE" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST SESCREE (SCREE TEST)                |\n")
      cat("==============================================\n\n")
      
      cat("Description: Le test du coude (scree test) examine le graphique des valeurs propres\n")
      cat("et identifie le point où la courbe s'aplatit (le 'coude').\n")
      cat(sprintf("Nombre de facteurs suggéré: %d\n\n", results$final[["SESCREE"]]))
    }
    
    # Test SMT (Test de dimensionnalité séquentielle)
    if (!is.null(results$SMT_details) || "SMT" %in% names(results$final)) {
      cat("\n==============================================\n")
      cat("| TEST SMT (DIMENSIONNALITÉ SÉQUENTIELLE)  |\n")
      cat("==============================================\n\n")
      
      if (!is.null(results$SMT_details)) {
        cat(paste(results$SMT_details, collapse = "\n"), "\n\n")
        
        # Afficher la matrice des p-values si disponible
        if (!is.null(results$SMT_pvalues)) {
          cat("DÉTAILS DES P-VALUES PAR MODÈLE:\n")
          cat("-------------------------------\n")
          print(results$SMT_pvalues)
          cat("\n")
          cat("Interprétation: Chaque ligne teste l'hypothèse que le modèle à N facteurs est suffisant.\n")
          cat("                 Une p-value > 0.05 indique que le modèle n'est pas rejeté.\n")
        }
      } else if ("SMT" %in% names(results$final)) {
        cat("Description: Ce test compare le modèle à k facteurs avec le modèle à k+1 facteurs\n")
        cat("en utilisant une approche séquentielle pour déterminer le nombre optimal de facteurs.\n")
        cat(sprintf("Nombre de facteurs suggéré: %d\n\n", results$final[["SMT"]]))
      }
    }
    
    # Synthèse et conclusion
    cat("\n==============================================\n")
    cat("| SYNTHÈSE DES RÉSULTATS                    |\n")
    cat("==============================================\n\n")
    
    # Compiler les résultats de tous les tests
    all_tests <- names(results$final)[names(results$final) != "NOTE"]
    n_tests <- length(all_tests)
    
    if (n_tests > 0) {
      test_values <- unlist(results$final[all_tests])
      
      # Mode statistique (valeur la plus fréquente)
      freq_table <- table(test_values)
      mode_factors <- as.numeric(names(freq_table)[which.max(freq_table)])
      mode_count <- max(freq_table)
      mode_percent <- round(mode_count / n_tests * 100)
      
      cat(sprintf("Synthèse des %d tests effectués:\n", n_tests))
      cat(sprintf("- Nombre de facteurs suggéré le plus souvent: %d (%d tests sur %d, soit %d%%)\n", 
                  mode_factors, mode_count, n_tests, mode_percent))
      
      # Comparer avec le nombre réel si disponible
      n_factors_real <- attr(sim_data(), "n.factors")
      if (!is.null(n_factors_real) && !is.na(n_factors_real)) {
        match_count <- sum(test_values == n_factors_real)
        match_percent <- round(match_count / n_tests * 100)
        
        cat(sprintf("- Tests suggérant le nombre réel de facteurs (%d): %d sur %d (%d%%)\n", 
                    n_factors_real, match_count, n_tests, match_percent))
        
        # Recommandation
        if (match_count / n_tests >= 0.5) {
          cat("\nRECOMMANDATION: Les tests confirment majoritairement la structure à", 
              n_factors_real, "facteurs des données.\n")
        } else if (mode_count / n_tests >= 0.5) {
          cat("\nRECOMMANDATION: Bien que les données aient été générées avec", n_factors_real, 
              "facteurs, la majorité des tests suggère une structure à", mode_factors, "facteurs.\n")
          cat("Cette divergence pourrait s'expliquer par des facteurs faibles ou par des spécificités\n")
          cat("dans la génération des données qui n'ont pas produit une structure factorielle claire.\n")
        } else {
          cat("\nRECOMMANDATION: Les tests ne convergent pas vers une solution unique,\n")
          cat("ce qui suggère une structure factorielle complexe ou ambiguë.\n")
          cat("Examinez attentivement le scree plot et considérez plusieurs solutions factorielles.\n")
        }
      } else {
        # Recommandation sans connaissance du nombre réel
        if (mode_count / n_tests >= 0.5) {
          cat("\nRECOMMANDATION: La majorité des tests suggère une structure à", 
              mode_factors, "facteurs.\n")
        } else {
          cat("\nRECOMMANDATION: Les tests ne convergent pas vers une solution unique.\n")
          cat("Considérez plusieurs solutions factorielles et utilisez des critères théoriques\n")
          cat("pour déterminer le nombre optimal de facteurs.\n")
        }
      }
    } else {
      cat("Aucun résultat de test dimensionnel n'est disponible pour fournir une interprétation.\n")
    }
    
    # Conseils généraux pour les tests
    cat("\nRAPPEL SUR L'INTERPRÉTATION DES TESTS:\n")
    cat("- Aucun test n'est parfait et chacun a ses propres limites.\n")
    cat("- La convergence de plusieurs tests renforce la confiance dans le résultat.\n")
    cat("- L'interprétation théorique reste essentielle pour déterminer le nombre final de facteurs.\n")
    cat("- Le critère des eigenvalues > 1 (Kaiser) tend à surestimer le nombre de facteurs.\n")
    cat("- Le test parallèle et la méthode Hull sont souvent considérés comme plus robustes.\n")
})

# Fonction de téléchargement des scores factoriels
output$download_factor_scores <- downloadHandler(
  filename = function() {
    paste("scores_factoriels_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv", sep = "")
  },
  content = function(file) {
    req(factor_scores_result())
    
    # Créer un dataframe de scores
    scores_df <- as.data.frame(factor_scores_result()$FactorScores)
    
    # Renommer les colonnes pour plus de clarté
    colnames(scores_df) <- paste0("Facteur_", 1:ncol(scores_df))
    
    # Ajouter une colonne d'ID
    scores_df <- cbind(ID = 1:nrow(scores_df), scores_df)
    
    # Écrire le CSV
    write.csv(scores_df, file, row.names = FALSE)
  }
)

# Fonction pour créer un graphique pour une seule métrique
create_single_metric_plot <- function(metric, comparison_df) {
  # Vérifier si la métrique est meilleure quand elle est plus basse
  is_lower_better <- function(metric) {
    return(metric %in% c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR"))
  }
  
  # Extraire les données pour cette métrique
  x_vals <- comparison_df$Factors
  y_vals <- comparison_df[[metric]]
  
  # Filtrer les valeurs NA
  valid_idx <- which(!is.na(y_vals))
  if (length(valid_idx) == 0) {
    return(plotly_empty() %>% 
             layout(title = paste("Aucune donnée disponible pour", metric)))
  }
  
  x_vals <- x_vals[valid_idx]
  y_vals <- y_vals[valid_idx]
  
  # Déterminer la meilleure valeur
  if (is_lower_better(metric)) {
    best_idx <- which.min(y_vals)
    best_text <- "Minimum (optimal)"
  } else {
    best_idx <- which.max(y_vals)
    best_text <- "Maximum (optimal)"
  }
  
  # S'assurer que best_idx est valide
  if (length(best_idx) == 0 || best_idx > length(x_vals)) {
    best_idx <- 1  # Par défaut
  }
  
  # Création du graphique
  p <- plot_ly() %>%
    add_trace(
      x = x_vals, 
      y = y_vals, 
      type = "scatter", 
      mode = "lines+markers",
      line = list(color = "blue", width = 2),
      marker = list(color = "blue", size = 10),
      name = metric,
      hoverinfo = "text",
      text = paste(
        "Facteurs:", x_vals, 
        "<br>", metric, ":", round(y_vals, 4)
      )
    ) %>%
    add_trace(
      x = x_vals[best_idx], 
      y = y_vals[best_idx], 
      type = "scatter", 
      mode = "markers",
      marker = list(color = "red", size = 14, symbol = "star"),
      name = best_text,
      hoverinfo = "text",
      text = paste(
        best_text, "<br>",
        "Facteurs:", x_vals[best_idx], "<br>",
        metric, ":", round(y_vals[best_idx], 4)
      )
    ) %>%
    layout(
      title = list(
        text = paste("Indice d'ajustement:", metric),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Nombre de facteurs",
        tickmode = "array",
        tickvals = x_vals
      ),
      yaxis = list(
        title = ifelse(is_lower_better(metric), 
                       "Plus petit = meilleur", 
                       "Plus grand = meilleur")
      ),
      showlegend = TRUE,
      legend = list(orientation = "h", y = 1.1)
    )
  
  return(p)
}

# Graphiques comparatifs - solution avec facet_wrap plutôt que subplot
output$compare_plots <- renderPlotly({
  req(compare_results())
  req(input$compare_metrics)
  
  # Vérifier qu'il y a au moins un indice sélectionné
  if (length(input$compare_metrics) == 0) {
    return(plotly_empty() %>% 
             layout(title = "Veuillez sélectionner au moins un indice à afficher"))
  }
  
  results <- compare_results()
  comparison_df <- results$comparison_df
  
  # Filtrer pour n'inclure que les métriques sélectionnées par l'utilisateur
  key_metrics <- input$compare_metrics
  
  # Vérifier que les métriques sélectionnées existent dans les données
  available_metrics <- key_metrics[key_metrics %in% names(comparison_df)]
  
  if (length(available_metrics) == 0) {
    return(plotly_empty() %>% 
             layout(title = "Aucun des indices sélectionnés n'est disponible dans les résultats"))
  }
  
  # Si on a un seul métrique, on utilise l'approche dédiée
  if (length(available_metrics) == 1) {
    return(create_single_metric_plot(available_metrics[1], comparison_df))
  }
  
  # ---------- NOUVELLE SOLUTION AVEC GGPLOT2 ET FACET_WRAP ----------
  
  # Préparer un dataframe long pour tous les indices sélectionnés
  plot_data <- data.frame()
  
  for (metric in available_metrics) {
    if (!all(is.na(comparison_df[[metric]]))) {
      metric_data <- data.frame(
        Factors = comparison_df$Factors,
        Value = comparison_df[[metric]],
        Metric = metric,
        stringsAsFactors = FALSE
      )
      
      # Déterminer le point optimal pour ce métrique
      is_lower_better <- metric %in% c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR")
      valid_values <- metric_data$Value[!is.na(metric_data$Value)]
      
      if (length(valid_values) > 0) {
        if (is_lower_better) {
          best_idx <- which.min(valid_values)
        } else {
          best_idx <- which.max(valid_values)
        }
        
        # Marquer les points optimaux
        metric_data$Optimal <- FALSE
        if (length(best_idx) > 0 && best_idx <= nrow(metric_data)) {
          rows_not_na <- which(!is.na(metric_data$Value))
          if (length(rows_not_na) >= best_idx) {
            metric_data$Optimal[rows_not_na[best_idx]] <- TRUE
          }
        }
      }
      
      # Ajouter une information sur le sens de l'optimisation
      metric_data$Better <- ifelse(is_lower_better, "Plus petit = meilleur", "Plus grand = meilleur")
      
      plot_data <- rbind(plot_data, metric_data)
    }
  }
  
  # Si aucune donnée valide
  if (nrow(plot_data) == 0) {
    return(plotly_empty() %>% 
             layout(title = "Aucune donnée disponible pour les indices sélectionnés"))
  }
  
  # Créer le graphique avec ggplot2 et facet_wrap
  require(ggplot2)
  
  g <- ggplot(plot_data, aes(x = Factors, y = Value, group = 1)) +
    geom_line(size = 0.8, color = "blue") + 
    geom_point(size = 2.5, color = "blue") +
    geom_point(data = subset(plot_data, Optimal == TRUE), 
               aes(x = Factors, y = Value),
               color = "red", size = 4, shape = 8) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    labs(
      title = "",
      x = "Nombre de facteurs",
      y = "Valeur"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "lightblue", color = "black"),
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 11),
      panel.spacing = unit(1, "lines")
    )
  
  # Convertir en plotly pour l'interactivité
  p <- ggplotly(g, tooltip = c("x", "y"))
  
  # Améliorer le texte au survol
  for (i in 1:length(p$x$data)) {
    if (!is.null(p$x$data[[i]]$marker) && !is.null(p$x$data[[i]]$marker$color)) {
      if (p$x$data[[i]]$marker$color == "rgba(255,0,0,1)") {
        # C'est un point optimal
        current_metric <- p$x$data[[i]]$name
        current_metric <- gsub("\\(", "", current_metric)
        current_metric <- gsub("\\)", "", current_metric)
        current_metric <- gsub("~", "", current_metric)
        
        for (metric in unique(plot_data$Metric)) {
          if (grepl(metric, current_metric)) {
            p$x$data[[i]]$text <- paste(
              "Point optimal<br>",
              "Facteurs:", p$x$data[[i]]$x,
              "<br>", metric, ":", round(p$x$data[[i]]$y, 4)
            )
            p$x$data[[i]]$hoverinfo <- "text"
            break
          }
        }
      } else if (p$x$data[[i]]$marker$color == "rgba(0,0,255,1)") {
        # C'est un point standard
        current_metric <- p$x$data[[i]]$name
        current_metric <- gsub("\\(", "", current_metric)
        current_metric <- gsub("\\)", "", current_metric)
        current_metric <- gsub("~", "", current_metric)
        
        for (metric in unique(plot_data$Metric)) {
          if (grepl(metric, current_metric)) {
            p$x$data[[i]]$text <- paste(
              "Facteurs:", p$x$data[[i]]$x,
              "<br>", metric, ":", round(p$x$data[[i]]$y, 4)
            )
            p$x$data[[i]]$hoverinfo <- "text"
            break
          }
        }
      }
    }
  }
  
  return(p)
})

# Fonction pour créer un graphique pour une seule métrique
create_single_metric_plot <- function(metric, comparison_df) {
  # Vérifier si la métrique est meilleure quand elle est plus basse
  is_lower_better <- function(metric) {
    return(metric %in% c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR"))
  }
  
  # Extraire les données pour cette métrique
  x_vals <- comparison_df$Factors
  y_vals <- comparison_df[[metric]]
  
  # Filtrer les valeurs NA
  valid_idx <- which(!is.na(y_vals))
  if (length(valid_idx) == 0) {
    return(plotly_empty() %>% 
             layout(title = paste("Aucune donnée disponible pour", metric)))
  }
  
  x_vals <- x_vals[valid_idx]
  y_vals <- y_vals[valid_idx]
  
  # Déterminer la meilleure valeur
  if (is_lower_better(metric)) {
    best_idx <- which.min(y_vals)
    best_text <- "Minimum (optimal)"
  } else {
    best_idx <- which.max(y_vals)
    best_text <- "Maximum (optimal)"
  }
  
  # S'assurer que best_idx est valide
  if (length(best_idx) == 0 || best_idx > length(x_vals)) {
    best_idx <- 1  # Par défaut
  }
  
  # Création du graphique
  p <- plot_ly() %>%
    add_trace(
      x = x_vals, 
      y = y_vals, 
      type = "scatter", 
      mode = "lines+markers",
      line = list(color = "blue", width = 2),
      marker = list(color = "blue", size = 10),
      name = metric,
      hoverinfo = "text",
      text = paste(
        "Facteurs:", x_vals, 
        "<br>", metric, ":", round(y_vals, 4)
      )
    ) %>%
    add_trace(
      x = x_vals[best_idx], 
      y = y_vals[best_idx], 
      type = "scatter", 
      mode = "markers",
      marker = list(color = "red", size = 14, symbol = "star"),
      name = best_text,
      hoverinfo = "text",
      text = paste(
        best_text, "<br>",
        "Facteurs:", x_vals[best_idx], "<br>",
        metric, ":", round(y_vals[best_idx], 4)
      )
    ) %>%
    layout(
      title = list(
        text = paste("Indice d'ajustement:", metric),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Nombre de facteurs",
        tickmode = "array",
        tickvals = x_vals
      ),
      yaxis = list(
        title = ifelse(is_lower_better(metric), 
                       "Plus petit = meilleur", 
                       "Plus grand = meilleur")
      ),
      showlegend = TRUE,
      legend = list(orientation = "h", y = 1.1)
    )
  
  return(p)
}

# Fonction pour créer un graphique pour une seule métrique
create_single_metric_plot <- function(metric, comparison_df) {
  # Vérifier si la métrique est meilleure quand elle est plus basse
  is_lower_better <- function(metric) {
    return(metric %in% c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR"))
  }
  
  # Extraire les données pour cette métrique
  x_vals <- comparison_df$Factors
  y_vals <- comparison_df[[metric]]
  
  # Filtrer les valeurs NA
  valid_idx <- which(!is.na(y_vals))
  if (length(valid_idx) == 0) {
    return(plotly_empty() %>% 
             layout(title = paste("Aucune donnée disponible pour", metric)))
  }
  
  x_vals <- x_vals[valid_idx]
  y_vals <- y_vals[valid_idx]
  
  # Déterminer la meilleure valeur
  if (is_lower_better(metric)) {
    best_idx <- which.min(y_vals)
    best_text <- "Minimum (optimal)"
  } else {
    best_idx <- which.max(y_vals)
    best_text <- "Maximum (optimal)"
  }
  
  # S'assurer que best_idx est valide
  if (length(best_idx) == 0 || best_idx > length(x_vals)) {
    best_idx <- 1  # Par défaut
  }
  
  # Création du graphique
  p <- plot_ly() %>%
    add_trace(
      x = x_vals, 
      y = y_vals, 
      type = "scatter", 
      mode = "lines+markers",
      line = list(color = "blue", width = 2),
      marker = list(color = "blue", size = 10),
      name = metric,
      hoverinfo = "text",
      text = paste(
        "Facteurs:", x_vals, 
        "<br>", metric, ":", round(y_vals, 4)
      )
    ) %>%
    add_trace(
      x = x_vals[best_idx], 
      y = y_vals[best_idx], 
      type = "scatter", 
      mode = "markers",
      marker = list(color = "red", size = 14, symbol = "star"),
      name = best_text,
      hoverinfo = "text",
      text = paste(
        best_text, "<br>",
        "Facteurs:", x_vals[best_idx], "<br>",
        metric, ":", round(y_vals[best_idx], 4)
      )
    ) %>%
    layout(
      title = list(
        text = paste("Indice d'ajustement:", metric),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Nombre de facteurs",
        tickmode = "array",
        tickvals = x_vals
      ),
      yaxis = list(
        title = ifelse(is_lower_better(metric), 
                       "Plus petit = meilleur", 
                       "Plus grand = meilleur")
      ),
      showlegend = TRUE,
      legend = list(orientation = "h", y = 1.1)
    )
  
  return(p)
}

# Graphiques comparatifs avec sélection d'indices et titres visibles
output$compare_plots_2 <- renderPlotly({
  req(compare_results())
  
  results <- compare_results()
  comparison_df <- results$comparison_df
  
  # Vérifier qu'on a bien des indices sélectionnés
  if (length(input$compare_metrics) == 0) {
    return(plotly_empty() %>% 
             layout(title = "Veuillez sélectionner au moins un indice à afficher"))
  }
  
  # Filtrer pour n'inclure que les métriques sélectionnées par l'utilisateur
  key_metrics <- input$compare_metrics
  
  # Vérifier que les métriques sélectionnées existent dans les données
  available_metrics <- key_metrics[key_metrics %in% names(comparison_df)]
  
  if (length(available_metrics) == 0) {
    return(plotly_empty() %>% 
             layout(title = "Aucun des indices sélectionnés n'est disponible dans les résultats"))
  }
  
  # Fonction pour déterminer si une métrique est meilleure quand elle est plus basse
  is_lower_better <- function(metric) {
    lower_better_metrics <- c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR")
    return(metric %in% lower_better_metrics)
  }
  
  # Créer des graphiques individuels pour chaque métrique
  plot_list <- list()
  
  for (i in seq_along(available_metrics)) {
    metric <- available_metrics[i]
    
    # Vérifier si la métrique a des valeurs non-NA
    if (!all(is.na(comparison_df[[metric]]))) {
      # Extraire les données pour cette métrique
      x_vals <- comparison_df$Factors
      y_vals <- comparison_df[[metric]]
      
      # Filtrer les valeurs NA
      valid_idx <- which(!is.na(y_vals))
      if (length(valid_idx) > 0) {
        x_vals <- x_vals[valid_idx]
        y_vals <- y_vals[valid_idx]
        
        # Déterminer la meilleure valeur
        if (is_lower_better(metric)) {
          best_idx <- which.min(y_vals)
          best_text <- "Valeur minimale (meilleure)"
        } else {
          best_idx <- which.max(y_vals)
          best_text <- "Valeur maximale (meilleure)"
        }
        
        # Créer un graphique pour cette métrique avec titre court
        p <- plot_ly() %>%
          add_trace(
            x = x_vals, 
            y = y_vals, 
            type = "scatter", 
            mode = "lines+markers",
            line = list(color = "blue"),
            marker = list(color = "blue", size = 10),
            name = metric,
            hoverinfo = "text",
            text = paste(
              "Facteurs:", x_vals, 
              "<br>", metric, ":", round(y_vals, 4),
              "<br>", ifelse(is_lower_better(metric), "Plus petit = meilleur", "Plus grand = meilleur")
            )
          ) %>%
          add_trace(
            x = x_vals[best_idx], 
            y = y_vals[best_idx], 
            type = "scatter", 
            mode = "markers",
            marker = list(color = "red", size = 15, symbol = "star"),
            name = best_text,
            hoverinfo = "text",
            text = paste(
              best_text, "<br>",
              "Facteurs:", x_vals[best_idx], "<br>",
              metric, ":", round(y_vals[best_idx], 4)
            )
          ) %>%
          layout(
            # Titre court et visible
            title = list(
              text = metric,
              font = list(size = 16),
              y = 0.9
            ),
            margin = list(t = 50),  # Marge supérieure pour le titre
            xaxis = list(
              title = "Nombre de facteurs", 
              tickmode = "array", 
              tickvals = x_vals
            ),
            yaxis = list(
              title = ifelse(is_lower_better(metric), "Plus petit = meilleur", "Plus grand = meilleur")
            ),
            showlegend = FALSE
          )
        
        plot_list[[i]] <- p
      }
    }
  }
  
  # Si aucun graphique n'a pu être créé, afficher un message
  if (length(plot_list) == 0) {
    return(plotly_empty() %>% 
             layout(title = "Aucune donnée disponible pour les indices sélectionnés"))
  }
  
  # Configurer le layout de base
  n_plots <- length(plot_list)
  subplot_rows <- ceiling(n_plots / 2)  # 2 colonnes par défaut
  
  # Combiner les graphiques individuels en sous-graphiques
  subplot_plot <- subplot(
    plot_list, 
    nrows = subplot_rows, 
    shareX = TRUE,
    margin = 0.12  # Espace plus grand entre les graphiques
  ) %>%
    layout(
      title = list(
        text = "Indices d'ajustement par nombre de facteurs",
        font = list(size = 18)
      ),
      showlegend = FALSE,
      margin = list(t = 100, b = 80, l = 80, r = 40)
    )
  
  # Retourner le graphique combiné
  subplot_plot
})

# Fonction pour créer un graphique pour une seule métrique
create_single_metric_plot <- function(metric, comparison_df) {
  # Vérifier si la métrique est meilleure quand elle est plus basse
  is_lower_better <- function(metric) {
    return(metric %in% c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR"))
  }
  
  # Extraire les données pour cette métrique
  x_vals <- comparison_df$Factors
  y_vals <- comparison_df[[metric]]
  
  # Filtrer les valeurs NA
  valid_idx <- which(!is.na(y_vals))
  if (length(valid_idx) == 0) {
    return(plotly_empty() %>% 
             layout(title = paste("Aucune donnée disponible pour", metric)))
  }
  
  x_vals <- x_vals[valid_idx]
  y_vals <- y_vals[valid_idx]
  
  # Déterminer la meilleure valeur
  if (is_lower_better(metric)) {
    best_idx <- which.min(y_vals)
    best_text <- "Minimum (optimal)"
  } else {
    best_idx <- which.max(y_vals)
    best_text <- "Maximum (optimal)"
  }
  
  # S'assurer que best_idx est valide
  if (length(best_idx) == 0 || best_idx > length(x_vals)) {
    best_idx <- 1  # Par défaut
  }
  
  # Création du graphique
  p <- plot_ly() %>%
    add_trace(
      x = x_vals, 
      y = y_vals, 
      type = "scatter", 
      mode = "lines+markers",
      line = list(color = "blue", width = 2),
      marker = list(color = "blue", size = 10),
      name = metric,
      hoverinfo = "text",
      text = paste(
        "Facteurs:", x_vals, 
        "<br>", metric, ":", round(y_vals, 4)
      )
    ) %>%
    add_trace(
      x = x_vals[best_idx], 
      y = y_vals[best_idx], 
      type = "scatter", 
      mode = "markers",
      marker = list(color = "red", size = 14, symbol = "star"),
      name = best_text,
      hoverinfo = "text",
      text = paste(
        best_text, "<br>",
        "Facteurs:", x_vals[best_idx], "<br>",
        metric, ":", round(y_vals[best_idx], 4)
      )
    ) %>%
    layout(
      title = list(
        text = paste("Indice d'ajustement:", metric),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Nombre de facteurs",
        tickmode = "array",
        tickvals = x_vals
      ),
      yaxis = list(
        title = ifelse(is_lower_better(metric), 
                       "Plus petit = meilleur", 
                       "Plus grand = meilleur")
      ),
      showlegend = TRUE,
      legend = list(orientation = "h", y = 1.1)
    )
  
  return(p)
}

# Fonction pour créer un graphique pour une seule métrique
create_single_metric_plot <- function(metric, comparison_df) {
  # Vérifier si la métrique est meilleure quand elle est plus basse
  is_lower_better <- function(metric) {
    return(metric %in% c("AIC", "BIC", "CAIC", "SABIC", "RMSEA", "RMSR"))
  }
  
  # Extraire les données pour cette métrique
  x_vals <- comparison_df$Factors
  y_vals <- comparison_df[[metric]]
  
  # Filtrer les valeurs NA
  valid_idx <- which(!is.na(y_vals))
  if (length(valid_idx) == 0) {
    return(plotly_empty() %>% 
             layout(title = paste("Aucune donnée disponible pour", metric)))
  }
  
  x_vals <- x_vals[valid_idx]
  y_vals <- y_vals[valid_idx]
  
  # Déterminer la meilleure valeur
  if (is_lower_better(metric)) {
    best_idx <- which.min(y_vals)
    best_text <- "Minimum (optimal)"
  } else {
    best_idx <- which.max(y_vals)
    best_text <- "Maximum (optimal)"
  }
  
  # S'assurer que best_idx est valide
  if (length(best_idx) == 0 || best_idx > length(x_vals)) {
    best_idx <- 1  # Par défaut
  }
  
  # Création du graphique
  p <- plot_ly() %>%
    add_trace(
      x = x_vals, 
      y = y_vals, 
      type = "scatter", 
      mode = "lines+markers",
      line = list(color = "blue", width = 2),
      marker = list(color = "blue", size = 10),
      name = metric,
      hoverinfo = "text",
      text = paste(
        "Facteurs:", x_vals, 
        "<br>", metric, ":", round(y_vals, 4)
      )
    ) %>%
    add_trace(
      x = x_vals[best_idx], 
      y = y_vals[best_idx], 
      type = "scatter", 
      mode = "markers",
      marker = list(color = "red", size = 14, symbol = "star"),
      name = best_text,
      hoverinfo = "text",
      text = paste(
        best_text, "<br>",
        "Facteurs:", x_vals[best_idx], "<br>",
        metric, ":", round(y_vals[best_idx], 4)
      )
    ) %>%
    layout(
      title = list(
        text = paste("Indice d'ajustement:", metric),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Nombre de facteurs",
        tickmode = "array",
        tickvals = x_vals
      ),
      yaxis = list(
        title = ifelse(is_lower_better(metric), 
                       "Plus petit = meilleur", 
                       "Plus grand = meilleur")
      ),
      showlegend = TRUE,
      legend = list(orientation = "h", y = 1.1)
    )
  
  return(p)
}

# Pour les statistiques descriptives
output$summary_stats <- renderPrint({
  req(sim_data())
  
  data_obj <- sim_data()
  response_data <- data_obj$response
  
  # Déterminer quelles variables analyser
  vars_to_analyze <- NULL
  if (input$stats_vars == "all") {
    vars_to_analyze <- colnames(response_data)
  } else if (input$stats_vars == "first") {
    n_vars <- min(input$num_stats_vars, ncol(response_data))
    vars_to_analyze <- colnames(response_data)[1:n_vars]
  } else if (input$stats_vars == "random") {
    n_vars <- min(input$num_stats_vars, ncol(response_data))
    vars_to_analyze <- sample(colnames(response_data), n_vars)
  }
  
  cat("STATISTIQUES DESCRIPTIVES DES DONNÉES\n")
  cat("===================================\n\n")
  
  # Afficher des statistiques de base
  cat("MESURES DE TENDANCE CENTRALE ET DE DISPERSION\n")
  cat("------------------------------------------\n")
  
  n_vars <- length(vars_to_analyze)
  subset_data <- response_data[, vars_to_analyze, drop = FALSE]
  
  # Calculer et afficher les statistiques pour chaque variable
  for (i in 1:n_vars) {
    var_name <- vars_to_analyze[i]
    var_data <- subset_data[, i]
    
    cat(sprintf("Variable %s:\n", var_name))
    cat(sprintf("  Moyenne: %.4f\n", mean(var_data)))
    cat(sprintf("  Médiane: %.4f\n", median(var_data)))
    cat(sprintf("  Écart-type: %.4f\n", sd(var_data)))
    cat(sprintf("  Min: %.4f\n", min(var_data)))
    cat(sprintf("  Max: %.4f\n", max(var_data)))
    
    # Statistiques avancées conditionnelles
    if (input$include_advanced_stats) {
      cat(sprintf("  Asymétrie (skewness): %.4f\n", moments::skewness(var_data)))
      cat(sprintf("  Aplatissement (kurtosis): %.4f\n", moments::kurtosis(var_data)))
      cat(sprintf("  Intervalle de confiance (95%%): [%.4f, %.4f]\n", 
                  mean(var_data) - 1.96 * sd(var_data)/sqrt(length(var_data)),
                  mean(var_data) + 1.96 * sd(var_data)/sqrt(length(var_data))))
    }
    cat("\n")
  }
})

# Génération automatique des tests de normalité (sans bouton)
observe({
  req(sim_data())
  req(input$normality_vars) # S'assurer que le contrôle existe et a une valeur
  
  # N'exécuter cette fonction que lorsque les données ou les sélections changent
  data_obj <- sim_data()
  response_data <- data_obj$response
  
  # Sélectionner les variables à tester
  selected_vars <- c()
  if (input$normality_vars == "all") {
    selected_vars <- colnames(response_data)
  } else if (input$normality_vars == "first") {
    n_vars <- min(input$num_norm_vars, ncol(response_data))
    selected_vars <- colnames(response_data)[1:n_vars]
  } else if (input$normality_vars == "random") {
    n_vars <- min(input$num_norm_vars, ncol(response_data))
    selected_vars <- sample(colnames(response_data), n_vars)
  }
  
  # Vérifier qu'on a des variables sélectionnées
  if (length(selected_vars) == 0) return()
  
  # Préparer un tableau pour stocker les résultats
  results <- data.frame(
    Variable = selected_vars,
    Shapiro_p = NA,
    Normalité = NA,
    stringsAsFactors = FALSE
  )
  
  # Exécuter le test de Shapiro-Wilk pour chaque variable
  for (i in 1:length(selected_vars)) {
    var_name <- selected_vars[i]
    var_data <- response_data[, var_name]
    
    # Test de Shapiro-Wilk
    test_result <- tryCatch(
      shapiro.test(var_data),
      error = function(e) {
        return(list(p.value = NA))
      }
    )
    
    # Stocker les résultats
    results$Shapiro_p[i] <- test_result$p.value
    results$Normalité[i] <- ifelse(test_result$p.value > 0.05, 
                                   "Normale", "Non-normale")
  }
  
  # Stocker les résultats pour l'affichage
  normality_results_data(list(
    results = results,
    selected_vars = selected_vars,
    data = response_data[, selected_vars, drop = FALSE]
  ))
})

observeEvent(input$deselect_all_distrib_vars, {
  updateCheckboxGroupInput(session, "selected_distrib_vars", selected = character(0))
})

# UI pour la sélection des variables de distribution
output$distrib_var_selector <- renderUI({
  req(sim_data())
  data_obj <- sim_data()
  var_names <- colnames(data_obj$response)
  
  # Limiter à 20 variables maximum pour éviter une interface trop chargée
  max_choices <- min(20, length(var_names))
  var_choices <- var_names[1:max_choices]
  
  # Sélectionner par défaut les 8 premières variables ou moins
  default_selected <- var_choices[1:min(8, length(var_choices))]
  
  tagList(
    checkboxGroupInput(
      "selected_distrib_vars",
      "Sélectionner les variables:",
      choices = var_choices,
      selected = default_selected,
      inline = FALSE
    ),
    actionLink("select_all_distrib_vars", "Sélectionner tout"),
    actionLink("deselect_all_distrib_vars", "Désélectionner tout")
  )
})

# Observateurs pour les boutons Sélectionner/Désélectionner tout
observeEvent(input$select_all_distrib_vars, {
  req(sim_data())
  data_obj <- sim_data()
  var_names <- colnames(data_obj$response)
  max_choices <- min(20, length(var_names))
  var_choices <- var_names[1:max_choices]
  updateCheckboxGroupInput(session, "selected_distrib_vars", selected = var_choices)
})

observeEvent(input$deselect_all_distrib_vars, {
  updateCheckboxGroupInput(session, "selected_distrib_vars", selected = character(0))
})

# Visualisation des distributions avec les variables sélectionnées manuellement
output$distribution_plot <- renderPlot({
  req(sim_data(), input$selected_distrib_vars)
  
  data_obj <- sim_data()
  response_data <- data_obj$response
  
  # Utiliser les variables sélectionnées directement
  vars_to_display <- input$selected_distrib_vars
  
  # Vérifier que vars_to_display n'est pas vide
  if (length(vars_to_display) == 0) {
    plot(0, 0, type = "n", xlab = "", ylab = "", main = "Aucune variable sélectionnée")
    return()
  }
  
  # Définir un nombre de colonnes raisonnable
  n_cols <- min(4, length(vars_to_display))
  n_rows <- ceiling(length(vars_to_display) / n_cols)
  
  # Créer les graphiques en fonction du type choisi
  if (input$distribution_plot_type == "hist") {
    # Histogrammes
    par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
    
    for (var in vars_to_display) {
      # Sécurité pour éviter les erreurs d'indice
      if (var %in% colnames(response_data)) {
        hist(response_data[, var], 
             main = var,
             xlab = "Valeur",
             col = "skyblue",
             border = "white",
             probability = TRUE)
        
        # Ajouter la densité seulement si possible
        tryCatch({
          lines(density(response_data[, var]), col = "red", lwd = 2)
        }, error = function(e) {
          # Ne rien faire si la densité ne peut pas être calculée
        })
      }
    }
  } else if (input$distribution_plot_type == "boxplot") {
    if (length(vars_to_display) == 1) {
      # Cas spécial pour une seule variable
      par(mar = c(4, 4, 3, 1))
      boxplot(response_data[, vars_to_display], 
              main = vars_to_display,
              ylab = "Valeur",
              col = "lightblue")
    } else {
      # Cas normal pour plusieurs variables
      par(mar = c(8, 4, 3, 1))
      # S'assurer que toutes les variables sélectionnées existent
      valid_vars <- vars_to_display[vars_to_display %in% colnames(response_data)]
      if (length(valid_vars) > 0) {
        boxplot(response_data[, valid_vars, drop = FALSE], 
                col = "lightblue", 
                main = "Boîtes à moustaches des variables sélectionnées",
                las = 2,
                cex.axis = 0.8)
      } else {
        plot(0, 0, type = "n", xlab = "", ylab = "", 
             main = "Aucune variable valide sélectionnée")
      }
    }
  } else if (input$distribution_plot_type == "density") {
    # Courbes de densité
    par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
    
    for (var in vars_to_display) {
      # Sécurité pour éviter les erreurs d'indice
      if (var %in% colnames(response_data)) {
        # Essayer de calculer la densité avec gestion d'erreur
        tryCatch({
          d <- density(response_data[, var])
          plot(d, 
               main = var,
               xlab = "Valeur",
               col = "blue",
               lwd = 2,
               fill = adjustcolor("skyblue", alpha.f = 0.3))
        }, error = function(e) {
          # Afficher un message d'erreur si la densité ne peut pas être calculée
          plot(0, 0, type = "n", 
               main = paste(var, "- Densité non calculable"),
               xlab = "", ylab = "")
          text(0, 0, "Impossible de calculer la densité\npour cette variable", col = "red")
        })
      }
    }
  } else if (input$distribution_plot_type == "qqplot") {
    # QQ-plots
    par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
    
    for (var in vars_to_display) {
      # Sécurité pour éviter les erreurs d'indice
      if (var %in% colnames(response_data)) {
        var_data <- response_data[, var]
        
        # Essayer de créer le QQ-plot avec gestion d'erreur
        tryCatch({
          qqnorm(var_data, main = paste("QQ Plot:", var),
                 pch = 16, col = adjustcolor("blue", alpha.f = 0.7))
          qqline(var_data, col = "red", lwd = 2)
          
          # Calculer et afficher le résultat du test de Shapiro-Wilk
          test_result <- shapiro.test(var_data)
          shapiro_p <- test_result$p.value
          
          shapiro_text <- paste("p =", round(shapiro_p, 3))
          verdict <- ifelse(shapiro_p > 0.05, "Normale", "Non-normale")
          text_color <- ifelse(shapiro_p > 0.05, "green4", "red")
          
          # Ajouter les résultats du test dans le coin du graphique
          legend("topleft", legend = c(shapiro_text, verdict),
                 text.col = c("black", text_color),
                 bty = "n", cex = 0.8)
        }, error = function(e) {
          # Afficher un message d'erreur si le QQ-plot ne peut pas être créé
          plot(0, 0, type = "n", 
               main = paste(var, "- QQ-plot non calculable"),
               xlab = "", ylab = "")
          text(0, 0, "Impossible de créer le QQ-plot\npour cette variable", col = "red")
        })
      }
    }
  }
})

# Mettre à jour le handler de téléchargement
output$download_distrib_plot <- downloadHandler(
  filename = function() {
    paste0("distribution_", input$distribution_plot_type, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  },
  content = function(file) {
    # Créer un fichier PNG plus grand pour une meilleure qualité
    png(file, width = 1200, height = 800, res = 120)
    
    # Utiliser le même code que pour l'affichage
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Utiliser les variables sélectionnées directement
    vars_to_display <- input$selected_distrib_vars
    
    # Vérifier que vars_to_display n'est pas vide
    if (length(vars_to_display) == 0) {
      plot(0, 0, type = "n", xlab = "", ylab = "", main = "Aucune variable sélectionnée")
      dev.off()
      return()
    }
    
    # Définir un nombre de colonnes raisonnable
    n_cols <- min(4, length(vars_to_display))
    n_rows <- ceiling(length(vars_to_display) / n_cols)
    
    # Créer les graphiques en fonction du type choisi - même logique que dans renderPlot
    if (input$distribution_plot_type == "hist") {
      par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
      for (var in vars_to_display) {
        if (var %in% colnames(response_data)) {
          hist(response_data[, var], 
               main = var,
               xlab = "Valeur",
               col = "skyblue",
               border = "white",
               probability = TRUE)
          tryCatch({
            lines(density(response_data[, var]), col = "red", lwd = 2)
          }, error = function(e) {})
        }
      }
    } else if (input$distribution_plot_type == "boxplot") {
      if (length(vars_to_display) == 1) {
        par(mar = c(4, 4, 3, 1))
        boxplot(response_data[[vars_to_display]], 
                main = vars_to_display,
                ylab = "Valeur",
                col = "lightblue")
      } else {
        par(mar = c(8, 4, 3, 1))
        valid_vars <- vars_to_display[vars_to_display %in% colnames(response_data)]
        if (length(valid_vars) > 0) {
          boxplot(response_data[, valid_vars, drop = FALSE], 
                  col = "lightblue", 
                  main = "Boîtes à moustaches des variables sélectionnées",
                  las = 2,
                  cex.axis = 0.8)
        } else {
          plot(0, 0, type = "n", xlab = "", ylab = "", 
               main = "Aucune variable valide sélectionnée")
        }
      }
    } else if (input$distribution_plot_type == "density") {
      par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
      for (var in vars_to_display) {
        # Sécurité pour éviter les erreurs d'indice
        if (var %in% colnames(response_data)) {
          # Essayer de calculer la densité avec gestion d'erreur
          tryCatch({
            d <- density(response_data[, var])
            plot(d, 
                 main = var,
                 xlab = "Valeur",
                 col = "blue",
                 lwd = 2,
                 fill = adjustcolor("skyblue", alpha.f = 0.3))
          }, error = function(e) {
            # Afficher un message d'erreur si la densité ne peut pas être calculée
            plot(0, 0, type = "n", 
                 main = paste(var, "- Densité non calculable"),
                 xlab = "", ylab = "")
            text(0, 0, "Impossible de calculer la densité\npour cette variable", col = "red")
          })
        }
      }
    } else if (input$distribution_plot_type == "qqplot") {
      par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
      for (var in vars_to_display) {
        # Sécurité pour éviter les erreurs d'indice
        if (var %in% colnames(response_data)) {
          var_data <- response_data[, var]
          
          # Essayer de créer le QQ-plot avec gestion d'erreur
          tryCatch({
            qqnorm(var_data, main = paste("QQ Plot:", var),
                   pch = 16, col = adjustcolor("blue", alpha.f = 0.7))
            qqline(var_data, col = "red", lwd = 2)
            
            # Calculer et afficher le résultat du test de Shapiro-Wilk
            test_result <- shapiro.test(var_data)
            shapiro_p <- test_result$p.value
            
            shapiro_text <- paste("p =", round(shapiro_p, 3))
            verdict <- ifelse(shapiro_p > 0.05, "Normale", "Non-normale")
            text_color <- ifelse(shapiro_p > 0.05, "green4", "red")
            
            # Ajouter les résultats du test dans le coin du graphique
            legend("topleft", legend = c(shapiro_text, verdict),
                   text.col = c("black", text_color),
                   bty = "n", cex = 0.8)
          }, error = function(e) {
            # Afficher un message d'erreur si le QQ-plot ne peut pas être créé
            plot(0, 0, type = "n", 
                 main = paste(var, "- QQ-plot non calculable"),
                 xlab = "", ylab = "")
            text(0, 0, "Impossible de créer le QQ-plot\npour cette variable", col = "red")
          })
        }
      }
    }
    
    dev.off()
  }
)

# UI pour sélectionner les variables à corréler
output$corr_var_selector <- renderUI({
  req(sim_data())
  data_obj <- sim_data()
  var_names <- colnames(data_obj$response)
  
  # Limiter à 25 variables maximum pour éviter une interface trop chargée
  max_choices <- min(25, length(var_names))
  var_choices <- var_names[1:max_choices]
  
  # Sélectionner par défaut les 10 premières variables ou moins
  default_selected <- var_choices[1:min(10, length(var_choices))]
  
  # Déterminer le nombre de colonnes en fonction du nombre de variables
  n_columns <- if(length(var_choices) > 15) 3 else if(length(var_choices) > 8) 2 else 1
  
  tagList(
    checkboxGroupInput(
      "selected_corr_vars",
      NULL, # Pas de label car déjà indiqué au-dessus
      choices = var_choices,
      selected = default_selected,
      inline = FALSE
    ) %>% tagAppendAttributes(class = paste0("multicol-", n_columns)),  # Ajouter une classe pour le style multi-colonnes
    
    # Ajouter du CSS pour le style multi-colonnes
    tags$style(HTML(paste0("
      .multicol-2 .shiny-options-group {
        -webkit-column-count: 2;
        -moz-column-count: 2;
        column-count: 2;
        -webkit-column-gap: 10px;
        -moz-column-gap: 10px;
        column-gap: 10px;
      }
      .multicol-3 .shiny-options-group {
        -webkit-column-count: 3;
        -moz-column-count: 3;
        column-count: 3;
        -webkit-column-gap: 10px;
        -moz-column-gap: 10px;
        column-gap: 10px;
      }
    "))),
    
    div(style = "display: flex; justify-content: space-between; margin-top: 10px;",
        actionLink("select_all_corr_vars", "Tout sélectionner"),
        actionLink("deselect_all_corr_vars", "Tout désélectionner")
    )
  )
})

# Observateurs pour les boutons Sélectionner/Désélectionner tout pour corrélations
observeEvent(input$select_all_corr_vars, {
  req(sim_data())
  data_obj <- sim_data()
  var_names <- colnames(data_obj$response)
  max_choices <- min(25, length(var_names))
  var_choices <- var_names[1:max_choices]
  updateCheckboxGroupInput(session, "selected_corr_vars", selected = var_choices)
})

observeEvent(input$deselect_all_corr_vars, {
  updateCheckboxGroupInput(session, "selected_corr_vars", selected = character(0))
})

output$item_corr_plot <- renderPlot({
  req(sim_data(), input$selected_corr_vars)
  
  data_obj <- sim_data()
  response_data <- data_obj$response
  
  # Vérifier qu'au moins 2 variables sont sélectionnées
  if (length(input$selected_corr_vars) < 2) {
    plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
    text(0, 0, "Sélectionnez au moins 2 variables pour afficher la matrice de corrélation", cex = 1.5)
    return()
  }
  
  # Sélectionner les variables à inclure
  var_subset <- response_data[, input$selected_corr_vars, drop = FALSE]
  
  # Calculer la matrice de corrélation
  cor_matrix <- cor(var_subset, method = input$corr_method)
  
  # Déterminer le titre
  n_vars_total <- ncol(response_data)
  n_vars_shown <- ncol(var_subset)
  title <- if (n_vars_shown < n_vars_total) {
    paste0("Matrice de corrélation (", n_vars_shown, " variables sur ", n_vars_total, ")")
  } else {
    "Matrice de corrélation complète"
  }
  
  # Définir la méthode de clustering à utiliser - CORRECTION DE LA PARTIE CLÉ
  clust_method <- "complete"  # Valeur par défaut
  if (input$corr_plot_order == "hclust" && !is.null(input$corr_clust_method)) {
    clust_method <- input$corr_clust_method
  }
  
  # Important: permettre le rafraîchissement quand la méthode de clustering change
  if (input$corr_plot_order == "hclust") {
    req(input$corr_clust_method)
  }
  
  # Visualiser avec corrplot et les options améliorées
  corrplot(
    cor_matrix, 
    method = input$corr_plot_method,
    type = input$corr_plot_type,
    order = input$corr_plot_order,
    addCoef.col = if(input$show_corr_values) "black" else NULL,
    number.cex = 0.7,
    tl.col = "black",
    tl.srt = 45,
    tl.cex = 0.8,
    diag = FALSE,
    title = title,
    mar = c(0, 0, 0.5, 0),
    cl.ratio = 0.2,
    col = COL2("RdBu", 200),
    hclust.method = clust_method  # Utilisation de la variable définie plus haut
  )
})

output$download_corr_plot <- downloadHandler(
  filename = function() {
    paste0("correlation_matrix_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  },
  content = function(file) {
    # Créer un fichier PNG plus grand pour une meilleure qualité
    png(file, width = 1200, height = 1000, res = 120)
    
    # Reproduire le même graphique
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Vérifier qu'au moins 2 variables sont sélectionnées
    if (length(input$selected_corr_vars) < 2) {
      plot(0, 0, type = "n", ann = FALSE, axes = FALSE)
      text(0, 0, "Sélectionnez au moins 2 variables pour afficher la matrice de corrélation", cex = 1.5)
    } else {
      # Sélectionner les variables à inclure
      var_subset <- response_data[, input$selected_corr_vars, drop = FALSE]
      
      # Calculer la matrice de corrélation
      cor_matrix <- cor(var_subset, method = input$corr_method)
      
      # Déterminer le titre
      n_vars_total <- ncol(response_data)
      n_vars_shown <- ncol(var_subset)
      title <- if (n_vars_shown < n_vars_total) {
        paste0("Matrice de corrélation (", n_vars_shown, " variables sur ", n_vars_total, ")")
      } else {
        "Matrice de corrélation complète"
      }
      
      # Définir la méthode de clustering à utiliser
      clust_method <- "complete"  # Valeur par défaut
      if (input$corr_plot_order == "hclust" && !is.null(input$corr_clust_method)) {
        clust_method <- input$corr_clust_method
      }
      
      # Visualiser avec corrplot et les options améliorées
      corrplot(
        cor_matrix, 
        method = input$corr_plot_method,
        type = input$corr_plot_type,
        order = input$corr_plot_order,
        addCoef.col = if(input$show_corr_values) "black" else NULL,
        number.cex = 0.7,
        tl.col = "black",
        tl.srt = 45,
        tl.cex = 0.8,
        diag = input$corr_plot_diag,
        title = title,
        mar = c(0, 0, 0.5, 0),
        cl.ratio = 0.2,
        col = COL2("RdBu", 200),
        hclust.method = clust_method
      )
    }
    
    dev.off()
  }
)




# Téléchargement de la matrice de corrélation au format CSV
output$download_corr_matrix <- downloadHandler(
  filename = function() {
    paste0("correlation_matrix_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  },
  content = function(file) {
    # Récupérer les données
    data_obj <- sim_data()
    response_data <- data_obj$response
    
    # Vérifier qu'au moins 2 variables sont sélectionnées
    if (length(input$selected_corr_vars) < 2) {
      # Créer un fichier vide avec un message d'erreur
      write.csv(data.frame(Erreur = "Sélectionnez au moins 2 variables pour la matrice de corrélation"), file)
    } else {
      # Sélectionner les variables à inclure
      var_subset <- response_data[, input$selected_corr_vars, drop = FALSE]
      
      # Calculer la matrice de corrélation
      cor_matrix <- cor(var_subset, method = input$corr_method)
      
      # Sauvegarder au format CSV
      write.csv(cor_matrix, file)
    }
  }
)

# Définition de l'UI conditionnelle pour les méthodes de clustering
output$clust_method_ui <- renderUI({
  # Ne montrer que si l'utilisateur a choisi le mode "hclust"
  if (input$corr_plot_order == "hclust") {
    selectInput("corr_clust_method", "Méthode de clustering:",
                choices = c("Complete" = "complete", 
                            "Ward.D" = "ward.D", 
                            "Ward.D2" = "ward.D2",
                            "Single" = "single", 
                            "Average (UPGMA)" = "average",
                            "McQuitty (WPGMA)" = "mcquitty", 
                            "Median (WPGMC)" = "median",
                            "Centroid (UPGMC)" = "centroid"),
                selected = "complete")
  }
})


# Réactive pour stocker le contenu du rapport
rapport_data <- reactiveVal(NULL)

# Observer pour la génération du rapport final
observeEvent(input$generate_report, {
  req(sim_data())
  
  # Afficher un indicateur de chargement
  withProgress(message = 'Génération du rapport en cours...', {
    # Récupérer les données de simulation et analyse
    data_obj <- sim_data()
    efa_res <- efa_result()
    dimtests_res <- dimtests_results()
    factor_scores_res <- factor_scores_result()
    reliability_res <- internal_consistency_results()
    
    # Déterminer le nombre optimal de facteurs à partir des tests de dimensionnalité
    optimal_nfactors <- input$rapport_nfactors  # Valeur par défaut
    
    if (!is.null(dimtests_res) && !is.null(dimtests_res$final)) {
      # Compter l'occurrence de chaque nombre de facteurs suggéré
      test_results <- unlist(dimtests_res$final[names(dimtests_res$final) != "NOTE"])
      
      if (length(test_results) > 0) {
        # Utiliser le mode (valeur la plus fréquente)
        freq_table <- table(test_results)
        optimal_nfactors <- as.numeric(names(freq_table)[which.max(freq_table)])
      }
    }
    
    # Construire le contenu du rapport
    report_content <- list(
      title = "Rapport d'analyse factorielle exploratoire",
      timestamp = Sys.time(),
      nfactors = optimal_nfactors,  # Utiliser le nombre optimal de facteurs
      user_selected_nfactors = input$rapport_nfactors,  # Garder une trace de la sélection de l'utilisateur
      data_summary = list(
        n_variables = ncol(data_obj$response),
        n_observations = nrow(data_obj$response),
        factorability = !is.null(factorability_results())
      ),
      structure_type = input$rapport_structure,
      dimtests_available = !is.null(dimtests_res),
      efa_available = !is.null(efa_res),
      scores_available = !is.null(factor_scores_res),
      reliability_available = !is.null(reliability_res) && input$rapport_include_reliability,
      include_communalities = input$rapport_include_communalities,
      include_methodology = input$rapport_include_methodology
    )
    
    # Stocker pour affichage
    rapport_data(report_content)
    
    # Notification
    showNotification("Rapport généré avec succès !", type = "message")
  })
})

# Rendu du contenu du rapport
output$rapport_contenu <- renderUI({
  req(rapport_data())
  report <- rapport_data()
  
  # Si aucun rapport n'est disponible
  if (is.null(report)) {
    return(div(
      class = "alert alert-info",
      icon("info-circle"), 
      strong("Aucun rapport disponible"),
      p("Configurez et générez un rapport en utilisant le panneau de gauche.")
    ))
  }
  
  # Header du rapport
  output <- div(
    h2(report$title),
    p(class = "text-muted", "Généré le ", format(report$timestamp, "%d/%m/%Y à %H:%M:%S")),
    hr(),
    
    # Résumé exécutif
    card(
      card_header(
        div(class = "d-flex align-items-center",
            span(icon("chart-pie"), style = "margin-right: 10px;"),
            h4("Résumé exécutif", class = "m-0"))
      ),
      card_body(
        p("Cette analyse factorielle exploratoire a été réalisée sur un ensemble de données simulées comportant", 
          strong(report$data_summary$n_variables), "variables et", 
          strong(report$data_summary$n_observations), "observations."),
        p("Le nombre optimal de facteurs retenu est de", strong(report$nfactors), 
          "sur la base des tests de dimensionnalité et de considérations théoriques."),
        p("La solution factorielle présente une structure", 
          if(report$structure_type == "simple") {
            "relativement simple, avec peu de charges croisées significatives entre les facteurs."
          } else {
            "complexe, avec plusieurs variables présentant des charges significatives sur plus d'un facteur."
          })
      )
    ),
    
    # Méthodologie (si demandée)
    if(report$include_methodology) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("tasks"), style = "margin-right: 10px;"),
              h4("Méthodologie", class = "m-0"))
        ),
        card_body(
          h5("Extraction des facteurs"),
          p("La méthode d'extraction utilisée est la factorisation en axes principaux (Principal Axis Factoring). Cette méthode est adaptée pour identifier les facteurs latents sous-jacents aux variables observées."),
          
          h5("Rotation"),
          p("Une rotation", 
            if(report$efa_available && !is.null(efa_result()$rotation)) {
              paste0(efa_result()$rotation, " ")
            } else {
              "Promax "
            },
            "a été appliquée pour optimiser la structure factorielle et faciliter l'interprétation des facteurs."),
          
          h5("Critères de rétention des facteurs"),
          p("Le nombre de facteurs a été déterminé en combinant plusieurs approches :"),
          tags$ul(
            tags$li("Tests statistiques : analyse parallèle, méthode Hull, et critère d'information empirique de Kaiser"),
            tags$li("Scree plot : examen visuel du graphique des valeurs propres"),
            tags$li("Interprétabilité : évaluation de la cohérence théorique de chaque facteur")
          )
        )
      )
    } else NULL,
    
    # Résultats des tests dimensionnels
    if(report$dimtests_available) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("ruler"), style = "margin-right: 10px;"),
              h4("Tests de dimensionnalité", class = "m-0"))
        ),
        card_body(
          p("Les tests de dimensionnalité ont été réalisés pour déterminer le nombre optimal de facteurs à extraire."),
          
          # Ajouter le screeplot
          h5("Scree Plot"),
          div(
            style = "margin: 20px 0;",
            screeplot_img()
          ),
          
          # Tableau récapitulatif des tests
          h5("Résultats des tests statistiques"),
          tags$table(
            class = "table table-striped table-bordered",
            tags$thead(
              tags$tr(
                tags$th("Test"),
                tags$th("Nombre de facteurs suggéré")
              )
            ),
            tags$tbody(
              lapply(names(dimtests_results()$final), function(test_name) {
                if (test_name != "NOTE") {
                  tags$tr(
                    tags$td(test_name),
                    tags$td(dimtests_results()$final[[test_name]])
                  )
                }
              })
            )
          ),
          
          p("Sur la base de ces résultats, une solution à", strong(report$nfactors), "facteurs a été retenue pour l'analyse."),
          
          # Explication du choix si différent des résultats des tests
          if (report$user_selected_nfactors != report$nfactors) {
            div(
              class = "alert alert-info",
              icon("info-circle"),
              strong("Note : "),
              "Bien que les tests suggèrent majoritairement", strong(report$nfactors), "facteurs, l'utilisateur a spécifié", 
              strong(report$user_selected_nfactors), "facteurs dans la configuration du rapport."
            )
          } else NULL
        )
      )
    } else NULL,
    
    # Résultats de l'analyse factorielle
    if(report$efa_available) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("project-diagram"), style = "margin-right: 10px;"),
              h4("Solution factorielle", class = "m-0"))
        ),
        card_body(
          h5("Matrice des saturations"),
          p("Cette matrice présente les charges factorielles, qui indiquent la force et la direction des relations entre les variables observées et les facteurs latents."),
          
          # Matrice simplifiée avec les charges principales
          p("Un seuil de 0.3 a été utilisé pour identifier les charges significatives."),
          
          # Section sur les communalités si demandée
          if(report$include_communalities) {
            tagList(
              h5("Communalités"),
              p("Les communalités indiquent la proportion de variance de chaque variable expliquée par l'ensemble des facteurs retenus."),
              p("Une communalité élevée (> 0.5) indique que la variable est bien représentée dans la solution factorielle."),
              # Ici on pourrait afficher un tableau des communalités si disponibles
              br()
            )
          } else NULL,
          
          h5("Interprétation des facteurs"),
          p("Sur la base des charges factorielles, les facteurs identifiés peuvent être interprétés comme suit :"),
          
          # Interprétation des facteurs (à personnaliser)
          tags$ul(
            tags$li(strong("Facteur 1:"), "Ce facteur regroupe des variables liées aux compétences cognitives et à la performance intellectuelle."),
            tags$li(strong("Facteur 2:"), "Ce facteur semble capturer des aspects de la personnalité liés à l'extraversion et aux interactions sociales."),
            if(report$nfactors > 2) tags$li(strong("Facteur 3:"), "Ce facteur reflète des dimensions émotionnelles comme l'anxiété et la stabilité émotionnelle.") else NULL,
            if(report$nfactors > 3) tags$li(strong("Facteur 4:"), "Ce facteur rassemble des variables liées à l'ouverture à l'expérience et à la créativité.") else NULL,
            if(report$nfactors > 4) tags$li(strong("Facteur 5:"), "Ce facteur semble représenter des aspects de la conscience et de la rigueur.") else NULL
          )
        )
      )
    } else NULL,
    
    # Résultats de l'analyse de fiabilité
    if(report$reliability_available) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("shield-alt"), style = "margin-right: 10px;"),
              h4("Analyse de fiabilité", class = "m-0"))
        ),
        card_body(
          if(!is.null(internal_consistency_results())) {
            ic_result <- internal_consistency_results()
            
            if(!is.null(ic_result$int.consist_scale)) {
              omega_value <- ic_result$int.consist_scale[1, "omega"]
              alpha_value <- ic_result$int.consist_scale[1, "alpha"]
              
              tagList(
                p("L'analyse de la cohérence interne des items a été réalisée afin d'évaluer la fiabilité des scores factoriels."),
                
                h5("Coefficients de fiabilité globaux"),
                tags$table(
                  class = "table table-striped table-sm",
                  tags$thead(
                    tags$tr(
                      tags$th("Coefficient"),
                      tags$th("Valeur"),
                      tags$th("Interprétation")
                    )
                  ),
                  tags$tbody(
                    tags$tr(
                      tags$td("Alpha de Cronbach"),
                      tags$td(round(alpha_value, 3)),
                      tags$td(
                        if(alpha_value >= 0.9) {
                          "Excellente (≥ 0.9)"
                        } else if(alpha_value >= 0.8) {
                          "Bonne (0.8-0.89)"
                        } else if(alpha_value >= 0.7) {
                          "Acceptable (0.7-0.79)"
                        } else if(alpha_value >= 0.6) {
                          "Questionnable (0.6-0.69)"
                        } else if(alpha_value >= 0.5) {
                          "Faible (0.5-0.59)"
                        } else {
                          "Inacceptable (< 0.5)"
                        }
                      )
                    ),
                    tags$tr(
                      tags$td("Oméga de McDonald"),
                      tags$td(round(omega_value, 3)),
                      tags$td(
                        if(omega_value >= 0.9) {
                          "Excellente (≥ 0.9)"
                        } else if(omega_value >= 0.8) {
                          "Bonne (0.8-0.89)"
                        } else if(omega_value >= 0.7) {
                          "Acceptable (0.7-0.79)"
                        } else if(omega_value >= 0.6) {
                          "Questionnable (0.6-0.69)"
                        } else {
                          "Insuffisante (< 0.6)"
                        }
                      )
                    )
                  )
                ),
                
                # Conclusion sur la fiabilité
                if(alpha_value >= 0.7 || omega_value >= 0.7) {
                  div(
                    class = "alert alert-success",
                    icon("check-circle"),
                    strong("Conclusion : "),
                    "Les scores factoriels présentent une fiabilité satisfaisante."
                  )
                } else if(alpha_value >= 0.6 || omega_value >= 0.6) {
                  div(
                    class = "alert alert-warning",
                    icon("exclamation-triangle"),
                    strong("Conclusion : "),
                    "Les scores factoriels présentent une fiabilité acceptable mais pourraient être améliorés."
                  )
                } else {
                  div(
                    class = "alert alert-danger",
                    icon("times-circle"),
                    strong("Conclusion : "),
                    "Les scores factoriels présentent une fiabilité insuffisante."
                  )
                }
              )
            } else {
              p("Résultats de fiabilité indisponibles ou incomplets.")
            }
          } else {
            p("L'analyse de fiabilité n'a pas été effectuée dans cette session.")
          }
        )
      )
    } else NULL,
    
    # Conclusion générale
    card(
      card_header(
        div(class = "d-flex align-items-center",
            span(icon("flag-checkered"), style = "margin-right: 10px;"),
            h4("Conclusion générale", class = "m-0"))
      ),
      card_body(
        p("L'analyse factorielle exploratoire a identifié une structure à", strong(report$nfactors), "facteurs."),
        
        p("Cette solution est soutenue par plusieurs tests de dimensionnalité et présente des propriétés psychométriques",
          if(!is.null(efa_result()) && !is.null(efa_result()$fit_coefs) && !is.null(efa_result()$fit_coefs$CFI) && efa_result()$fit_coefs$CFI > 0.9) {
            "satisfaisantes"
          } else {
            "qui nécessitent une interprétation prudente"
          },
          "."
        ),
        
        if(report$structure_type == "simple") {
          p("La structure factorielle est simple et bien définie, avec peu de charges croisées significatives entre les facteurs, ce qui facilite l'interprétation théorique des dimensions identifiées.")
        } else {
          p("La structure factorielle présente une certaine complexité, avec plusieurs variables présentant des charges significatives sur plus d'un facteur. Cette structure complexe peut refléter la nature multidimensionnelle des construits mesurés ou indiquer des problèmes de spécification du modèle.")
        },
        
        p("Les résultats de cette analyse peuvent être utilisés pour :", 
          tags$ul(
            tags$li("Développer ou raffiner des instruments de mesure"),
            tags$li("Explorer la structure sous-jacente des données"),
            tags$li("Guider des analyses confirmatoires ultérieures")
          )
        ),
        
        # Recommandations finales - Version simplifiée pour éviter l'erreur
        h5("Recommandations"),
        tags$ul(
          tags$li("Valider cette structure factorielle sur un échantillon indépendant"),
          tags$li("Comparer cette solution à des modèles théoriques existants dans la littérature"),
          tags$li("Envisager une analyse factorielle confirmatoire pour tester plus rigoureusement cette structure")
        )
      )
    )
  )
  
  return(output)
})

# Téléchargement du rapport au format HTML
# Modification de la fonction de téléchargement du rapport au format HTML
output$download_report <- downloadHandler(
  filename = function() {
    paste0("rapport_factoriel_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
  },
  content = function(file) {
    # Isoler les valeurs réactives
    report <- isolate(rapport_data())
    
    # Créer un contenu HTML formaté sans utiliser isolate(output$rapport_contenu())
    html_content <- tags$html(
      tags$head(
        tags$title("Rapport d'analyse factorielle"),
        tags$style(HTML("
          body { font-family: Arial, sans-serif; margin: 30px; line-height: 1.6; }
          h1, h2, h3, h4 { color: #2c3e50; }
          .card { border: 1px solid #ddd; border-radius: 5px; margin-bottom: 20px; }
          .card-header { background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #ddd; }
          .card-body { padding: 15px; }
          .table { width: 100%; border-collapse: collapse; }
          .table th, .table td { padding: 8px; border: 1px solid #ddd; text-align: left; }
          .table th { background-color: #f2f2f2; }
          .alert { padding: 15px; border-radius: 4px; margin-bottom: 15px; }
          .alert-success { background-color: #d4edda; color: #155724; }
          .alert-warning { background-color: #fff3cd; color: #856404; }
          .alert-danger { background-color: #f8d7da; color: #721c24; }
          .text-success { color: #28a745; }
          .text-danger { color: #dc3545; }
        "))
      ),
      tags$body(
        tags$h1("Rapport d'analyse factorielle exploratoire"),
        tags$p(class = "text-muted", "Généré le ", format(Sys.time(), "%d/%m/%Y à %H:%M:%S")),
        tags$hr(),
        
        # Fonction isolate pour accéder aux valeurs réactives sans déclencher de réaction
        # Recréer manuellement le contenu du rapport
        isolate({
          if (is.null(report)) {
            tags$div(
              class = "alert alert-info",
              tags$strong("Aucun rapport disponible"),
              tags$p("Configurez et générez un rapport en utilisant le panneau de gauche.")
            )
          } else {
            # Récréer le contenu du rapport (similaire à la fonction renderUI)
            # Header du rapport
            div(
              h2(report$title),
              p(class = "text-muted", "Généré le ", format(report$timestamp, "%d/%m/%Y à %H:%M:%S")),
              hr(),
              
              # Résumé exécutif
              div(
                class = "card",
                div(
                  class = "card-header",
                  h4("Résumé exécutif")
                ),
                div(
                  class = "card-body",
                  p("Cette analyse factorielle exploratoire a été réalisée sur un ensemble de données simulées comportant", 
                    strong(report$data_summary$n_variables), "variables et", 
                    strong(report$data_summary$n_observations), "observations."),
                  p("Le nombre optimal de facteurs retenu est de", strong(report$nfactors), 
                    "sur la base des tests de dimensionnalité et de considérations théoriques."),
                  p("La solution factorielle présente une structure", 
                    if(report$structure_type == "simple") {
                      "relativement simple, avec peu de charges croisées significatives entre les facteurs."
                    } else {
                      "complexe, avec plusieurs variables présentant des charges significatives sur plus d'un facteur."
                    })
                )
              ),
              
              # Méthodologie (si demandée)
              if(report$include_methodology) {
                div(
                  class = "card",
                  div(
                    class = "card-header",
                    h4("Méthodologie")
                  ),
                  div(
                    class = "card-body",
                    h5("Extraction des facteurs"),
                    p("La méthode d'extraction utilisée est la factorisation en axes principaux (Principal Axis Factoring). Cette méthode est adaptée pour identifier les facteurs latents sous-jacents aux variables observées."),
                    
                    h5("Rotation"),
                    p("Une rotation", 
                      if(report$efa_available && !is.null(efa_result()$rotation)) {
                        paste0(efa_result()$rotation, " ")
                      } else {
                        "Promax "
                      },
                      "a été appliquée pour optimiser la structure factorielle et faciliter l'interprétation des facteurs."),
                    
                    h5("Critères de rétention des facteurs"),
                    p("Le nombre de facteurs a été déterminé en combinant plusieurs approches :"),
                    tags$ul(
                      tags$li("Tests statistiques : analyse parallèle, méthode Hull, et critère d'information empirique de Kaiser"),
                      tags$li("Scree plot : examen visuel du graphique des valeurs propres"),
                      tags$li("Interprétabilité : évaluation de la cohérence théorique de chaque facteur")
                    )
                  )
                )
              } else NULL,
              
              # Dans la partie qui génère le HTML pour le téléchargement, ajoutez après la section du modèle factoriel :
              
              # Section des résidus
              if(report$efa_available) {
                # Calculer les résidus pour le rapport
                residuals_data <- isolate(residual_matrix())
                
                div(
                  class = "card",
                  div(
                    class = "card-header",
                    h4("Analyse des résidus du modèle")
                  ),
                  div(
                    class = "card-body",
                    p("L'analyse des résidus permet d'évaluer l'adéquation du modèle factoriel aux données observées."),
                    
                    h5("Indicateurs d'ajustement basés sur les résidus"),
                    tags$table(
                      class = "table table-striped table-sm",
                      tags$thead(
                        tags$tr(
                          tags$th("Indicateur"),
                          tags$th("Valeur"),
                          tags$th("Interprétation")
                        )
                      ),
                      tags$tbody(
                        tags$tr(
                          tags$td("RMSR (Root Mean Square of Residuals)"),
                          tags$td(round(residuals_data$rmsr, 3)),
                          tags$td(
                            if(residuals_data$rmsr < 0.05) {
                              "Excellent (<0.05)"
                            } else if(residuals_data$rmsr < 0.08) {
                              "Bon (0.05-0.08)"
                            } else if(residuals_data$rmsr < 0.10) {
                              "Acceptable (0.08-0.10)"
                            } else {
                              "Médiocre (>0.10)"
                            }
                          )
                        ),
                        tags$tr(
                          tags$td("Résidu maximal (valeur absolue)"),
                          tags$td(round(residuals_data$max_residual, 3)),
                          tags$td(
                            if(residuals_data$max_residual < 0.10) {
                              "Excellent (<0.10)"
                            } else if(residuals_data$max_residual < 0.15) {
                              "Bon (0.10-0.15)"
                            } else {
                              "Préoccupant (>0.15)"
                            }
                          )
                        ),
                        tags$tr(
                          tags$td("Proportion de résidus >0.05 (valeur absolue)"),
                          tags$td(paste0(round(residuals_data$prop_large_residuals * 100, 1), "%")),
                          tags$td(
                            if(residuals_data$prop_large_residuals < 0.10) {
                              "Excellent (<10%)"
                            } else if(residuals_data$prop_large_residuals < 0.25) {
                              "Acceptable (10-25%)"
                            } else {
                              "Problématique (>25%)"
                            }
                          )
                        )
                      )
                    ),
                    
                    # Conclusion sur les résidus
                    if(residuals_data$rmsr < 0.08 && residuals_data$max_residual < 0.15 && residuals_data$prop_large_residuals < 0.25) {
                      div(
                        class = "alert alert-success",
                        tags$i(class = "fas fa-check-circle", style = "margin-right: 10px;"),
                        strong("Conclusion : "),
                        "Le modèle à ", strong(report$nfactors), " facteurs reproduit bien la structure de corrélation des données observées."
                      )
                    } else if(residuals_data$rmsr < 0.10 && residuals_data$max_residual < 0.20 && residuals_data$prop_large_residuals < 0.40) {
                      div(
                        class = "alert alert-warning",
                        tags$i(class = "fas fa-exclamation-triangle", style = "margin-right: 10px;"),
                        strong("Conclusion : "),
                        "Le modèle reproduit de façon acceptable la structure de corrélation, mais certaines relations entre variables ne sont pas parfaitement capturées."
                      )
                    } else {
                      div(
                        class = "alert alert-danger",
                        tags$i(class = "fas fa-times-circle", style = "margin-right: 10px;"),
                        strong("Conclusion : "),
                        "Le modèle ne reproduit pas adéquatement la structure de corrélation. Considérez d'augmenter le nombre de facteurs ou de modifier la spécification du modèle."
                      )
                    },
                    
                    # Visualisation optionnelle - description textuelle des résidus importants
                    h5("Résidus importants"),
                    if(residuals_data$max_residual > 0.10) {
                      # Identifier les paires de variables avec des résidus importants
                      residuals_mat <- residuals_data$residuals
                      var_names <- rownames(residuals_mat)
                      if(is.null(var_names)) {
                        var_names <- paste0("Variable ", 1:nrow(residuals_mat))
                      }
                      
                      # Trouver les résidus élevés (>0.10 en valeur absolue)
                      high_residuals <- which(abs(residuals_mat) > 0.10 & lower.tri(residuals_mat), arr.ind = TRUE)
                      
                      if(nrow(high_residuals) > 0) {
                        # Créer une table des résidus importants
                        residuals_table <- data.frame(
                          Variable1 = var_names[high_residuals[, 1]],
                          Variable2 = var_names[high_residuals[, 2]],
                          Residual = residuals_mat[high_residuals]
                        )
                        
                        # Trier par valeur absolue des résidus
                        residuals_table <- residuals_table[order(-abs(residuals_table$Residual)), ]
                        
                        # Limiter à 10 résidus maximum pour la lisibilité
                        if(nrow(residuals_table) > 10) {
                          residuals_table <- residuals_table[1:10, ]
                        }
                        
                        # Créer un tableau HTML
                        tags$table(
                          class = "table table-striped table-sm",
                          tags$thead(
                            tags$tr(
                              tags$th("Variable 1"),
                              tags$th("Variable 2"),
                              tags$th("Résidu")
                            )
                          ),
                          tags$tbody(
                            lapply(1:nrow(residuals_table), function(i) {
                              tags$tr(
                                tags$td(residuals_table$Variable1[i]),
                                tags$td(residuals_table$Variable2[i]),
                                tags$td(round(residuals_table$Residual[i], 3))
                              )
                            })
                          )
                        )
                      } else {
                        p("Aucun résidu important (>0.10) n'a été identifié.")
                      }
                    } else {
                      p("Tous les résidus sont inférieurs à 0.10 en valeur absolue, indiquant une bonne reproduction des corrélations observées.")
                    },
                    
                    # Interprétation et recommandations
                    h5("Implications pour le modèle"),
                    p("L'analyse des résidus suggère que :", 
                      if(residuals_data$rmsr < 0.08) {
                        "le nombre de facteurs retenus semble approprié pour représenter la structure des données."
                      } else if(residuals_data$prop_large_residuals > 0.25) {
                        "certaines relations entre variables ne sont pas bien expliquées par le modèle actuel. Un modèle avec plus de facteurs pourrait être envisagé."
                      } else {
                        "le modèle pourrait être amélioré, mais offre une représentation acceptable de la structure des données."
                      }
                    )
                  )
                )
              } else NULL,
              
              # Ajouter ici les sections restantes du rapport...
              # (sections conditionnelles, conclusion, etc.)
              
              # Conclusion générale (toujours incluse)
              div(
                class = "card",
                div(
                  class = "card-header",
                  h4("Conclusion générale")
                ),
                div(
                  class = "card-body",
                  p("L'analyse factorielle exploratoire a identifié une structure à", strong(report$nfactors), "facteurs."),
                  
                  p("Les résultats de cette analyse peuvent être utilisés pour :", 
                    tags$ul(
                      tags$li("Développer ou raffiner des instruments de mesure"),
                      tags$li("Explorer la structure sous-jacente des données"),
                      tags$li("Guider des analyses confirmatoires ultérieures")
                    )
                  ),
                  
                  h5("Recommandations"),
                  tags$ul(
                    tags$li("Valider cette structure factorielle sur un échantillon indépendant"),
                    tags$li("Comparer cette solution à des modèles théoriques existants dans la littérature"),
                    tags$li("Envisager une analyse factorielle confirmatoire pour tester plus rigoureusement cette structure")
                  )
                )
              )
            )
          }
        })
      )
    )
    
    # Écrire le fichier HTML complet
    html_str <- as.character(html_content)
    writeLines(html_str, file)
  }
)

screeplot_img <- reactive({
  req(dimtests_results())
  
  # Créer un fichier temporaire pour l'image
  temp_file <- tempfile(fileext = ".png")
  
  # Sauvegarder le graphique dans le fichier temporaire
  png(temp_file, width = 800, height = 500, res = 100)
  
  # Récupérer les valeurs propres
  results <- dimtests_results()
  
  # Utiliser les valeurs propres appropriées
  if (!is.null(results$parallel_data) && !is.null(results$parallel_data$observed)) {
    eigenvalues <- results$parallel_data$observed
  } else {
    eigenvalues <- results$eigenvalues
  }
  
  # Créer le screeplot de base
  n_factors <- length(eigenvalues)
  factor_numbers <- 1:n_factors
  
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  # Tracer le screeplot
  plot(factor_numbers, eigenvalues, type = "b", pch = 19, col = "blue",
       xlab = "Numéro de facteur", ylab = "Valeur propre",
       main = "Scree Plot", cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.1)
  
  # Ajouter la ligne de Kaiser
  abline(h = 1, col = "darkgreen", lty = 2, lwd = 2)
  text(n_factors * 0.8, 1.1, "Critère de Kaiser (λ = 1)", col = "darkgreen", pos = 3)
  
  # Ajouter les données du test parallèle si disponibles
  if (!is.null(results$parallel_data) && !is.null(results$parallel_data$Percentile)) {
    percentile_values <- results$parallel_data$Percentile
    lines(1:length(percentile_values), percentile_values, col = "red", lty = 2, type = "b", pch = 17)
    legend("topright", 
           legend = c("Valeurs propres observées", "Test parallèle (95%)"),
           col = c("blue", "red"), lty = c(1, 2), pch = c(19, 17),
           cex = 0.8, bg = "white")
  }
  
  # Fermer le périphérique graphique
  dev.off()
  
  # Lire le fichier en base64
  base64_img <- base64enc::base64encode(temp_file)
  
  # Construire la balise img
  img_tag <- tags$img(
    src = paste0("data:image/png;base64,", base64_img),
    style = "width: 100%; max-width: 800px; margin: 0 auto; display: block;"
  )
  
  return(img_tag)
})


output$rapport_contenu <- renderUI({
  req(rapport_data())
  report <- rapport_data()
  
  div(
    h2(report$title),
    p(class = "text-muted", "Généré le ", format(report$timestamp, "%d/%m/%Y à %H:%M:%S")),
    hr(),
    
    # Résumé exécutif 
    card(
      card_header(
        div(class = "d-flex align-items-center",
            span(icon("info-circle"), style = "margin-right: 10px;"),
            h4("Résumé exécutif", class = "m-0"))
      ),
      card_body(
        p("Cette analyse factorielle exploratoire a été réalisée sur un ensemble de données simulées comportant", 
          strong(report$data_summary$n_variables), "variables et", 
          strong(report$data_summary$n_observations), "observations."),
        p("Le nombre optimal de facteurs retenu est de", strong(report$nfactors), 
          "sur la base des tests de dimensionnalité et de considérations théoriques."),
        p("La solution factorielle présente une structure", 
          if(report$structure_type == "simple") {
            "relativement simple, avec peu de charges croisées significatives entre les facteurs."
          } else {
            "complexe, avec plusieurs variables présentant des charges significatives sur plus d'un facteur."
          })
      )
    ),
    
    # Méthodologie (si demandée)
    if(report$include_methodology) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("cogs"), style = "margin-right: 10px;"),
              h4("Méthodologie", class = "m-0"))
        ),
        card_body(
          h5("Extraction des facteurs"),
          p("La méthode d'extraction utilisée est la factorisation en axes principaux (Principal Axis Factoring). Cette méthode est adaptée pour identifier les facteurs latents sous-jacents aux variables observées."),
          
          h5("Rotation"),
          p("Une rotation", 
            if(report$efa_available && !is.null(efa_result()$rotation)) {
              paste0(efa_result()$rotation, " ")
            } else {
              "Promax "
            },
            "a été appliquée pour optimiser la structure factorielle et faciliter l'interprétation des facteurs."),
          
          h5("Critères de rétention des facteurs"),
          p("Le nombre de facteurs a été déterminé en combinant plusieurs approches :"),
          tags$ul(
            tags$li("Tests statistiques : analyse parallèle, méthode Hull, et critère d'information empirique de Kaiser"),
            tags$li("Scree plot : examen visuel du graphique des valeurs propres"),
            tags$li("Interprétabilité : évaluation de la cohérence théorique de chaque facteur")
          )
        )
      )
    } else NULL,
    
    # Tests de dimensionnalité
    if(report$dimtests_available) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("ruler"), style = "margin-right: 10px;"),
              h4("Tests de dimensionnalité", class = "m-0"))
        ),
        card_body(
          p("Les tests de dimensionnalité ont été réalisés pour déterminer le nombre optimal de facteurs à extraire."),
          
          # Scree Plot (si disponible)
          if(!is.null(screeplot_img())) {
            div(
              h5("Scree Plot"),
              div(
                style = "margin: 20px 0;",
                screeplot_img()
              )
            )
          } else NULL,
          
          # Tableau récapitulatif des tests
          h5("Résultats des tests statistiques"),
          tags$table(
            class = "table table-striped table-bordered",
            tags$thead(
              tags$tr(
                tags$th("Test"),
                tags$th("Nombre de facteurs suggéré")
              )
            ),
            tags$tbody(
              lapply(names(dimtests_results()$final), function(test_name) {
                if (test_name != "NOTE") {
                  tags$tr(
                    tags$td(test_name),
                    tags$td(dimtests_results()$final[[test_name]])
                  )
                }
              })
            )
          ),
          
          p("Sur la base de ces résultats, une solution à", strong(report$nfactors), "facteurs a été retenue pour l'analyse."),
          
          # Explication du choix si différent des résultats des tests
          if (!is.null(report$user_selected_nfactors) && report$user_selected_nfactors != report$nfactors) {
            div(
              class = "alert alert-info",
              icon("info-circle"),
              strong("Note : "),
              "Bien que les tests suggèrent majoritairement", strong(report$nfactors), "facteurs, l'utilisateur a spécifié", 
              strong(report$user_selected_nfactors), "facteurs dans la configuration du rapport."
            )
          } else NULL
        )
      )
    } else NULL,
    
    # Structure factorielle
    if(report$efa_available) {
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("table"), style = "margin-right: 10px;"),
              h4("Structure factorielle", class = "m-0"))
        ),
        card_body(
          p("La solution à", strong(report$nfactors), "facteurs présente la structure suivante :"),
          
          # Matrice des charges factorielles
          renderTable({
            # Récupérer les charges factorielles
            loadings_table <- efa_result()$loadings
            
            # Formater les charges factorielles (arrondir et afficher seulement les charges significatives)
            loadings_table <- round(loadings_table, 3)
            loadings_table[abs(loadings_table) < 0.3] <- ""
            
            # Renommer les colonnes
            colnames(loadings_table) <- paste0("Facteur ", 1:ncol(loadings_table))
            
            # Ajouter les communalités si demandé
            if(report$include_communalities && !is.null(efa_result()$communalities)) {
              loadings_table <- cbind(loadings_table, h2 = round(efa_result()$communalities, 3))
            }
            
            # Retourner le tableau formaté
            return(loadings_table)
          }, rownames = TRUE, align = "c", width = "100%"),
          
          # Notes sur les charges factorielles
          tags$small(
            tags$em("Note : "),
            "Seules les charges factorielles ≥ 0.3 sont affichées. ",
            if(report$include_communalities) {
              "La colonne h² représente les communalités."
            } else {
              ""
            }
          ),
          
          # Interprétation de la structure factorielle
          h5("Interprétation des facteurs"),
          p("Sur la base des charges factorielles, les facteurs peuvent être interprétés comme suit :"),
          
          tags$ul(
            lapply(1:report$nfactors, function(i) {
              # Récupérer les charges factorielles pour ce facteur
              factor_loadings <- efa_result()$loadings[, i]
              
              # Trouver les variables avec des charges importantes
              important_vars <- which(abs(factor_loadings) >= 0.4)
              
              # Si aucune variable importante, retourner un message
              if(length(important_vars) == 0) {
                return(tags$li(
                  tags$strong(paste0("Facteur ", i, ": ")),
                  "Aucune variable avec une charge ≥ 0.4 sur ce facteur."
                ))
              }
              
              # Déterminer le signe dominant
              dominant_sign <- sign(mean(factor_loadings[important_vars]))
              
              # Créer un message d'interprétation
              tags$li(
                tags$strong(paste0("Facteur ", i, ": ")),
                "Caractérisé principalement par ",
                if(length(important_vars) <= 3) {
                  paste0(
                    paste(
                      ifelse(factor_loadings[important_vars] * dominant_sign < 0, "faible ", "forte "),
                      rownames(efa_result()$loadings)[important_vars],
                      collapse = ", "
                    ),
                    "."
                  )
                } else {
                  paste0(
                    length(important_vars),
                    " variables avec des charges ≥ 0.4, dont les plus importantes sont ",
                    paste(
                      rownames(efa_result()$loadings)[important_vars[order(abs(factor_loadings[important_vars]), decreasing = TRUE)][1:min(3, length(important_vars))]], 
                      collapse = ", "
                    ),
                    "."
                  )
                }
              )
            })
          ),
          
          # Corrélations entre facteurs (si rotation oblique)
          if(!is.null(efa_result()$phi) && !identical(efa_result()$phi, diag(report$nfactors))) {
            div(
              h5("Corrélations entre facteurs"),
              renderTable({
                phi_matrix <- efa_result()$phi
                rownames(phi_matrix) <- colnames(phi_matrix) <- paste0("F", 1:ncol(phi_matrix))
                round(phi_matrix, 3)
              }, rownames = TRUE, align = "c", width = "100%"),
              p("Ces corrélations entre facteurs indiquent", 
                if(max(abs(efa_result()$phi[lower.tri(efa_result()$phi)])) > 0.3) {
                  "une structure factorielle où certaines dimensions sont substantiellement corrélées, justifiant l'utilisation d'une rotation oblique."
                } else {
                  "une structure factorielle où les dimensions sont relativement indépendantes."
                }
              )
            )
          } else NULL
        )
      )
    } else NULL,
    
    # Analyse des résidus (si l'EFA est disponible)
    if(report$efa_available) {
      # Calculer les résidus pour le rapport
      residuals_data <- residual_matrix()
      
      card(
        card_header(
          div(class = "d-flex align-items-center",
              span(icon("check-double"), style = "margin-right: 10px;"),
              h4("Analyse des résidus du modèle", class = "m-0"))
        ),
        card_body(
          p("L'analyse des résidus permet d'évaluer l'adéquation du modèle factoriel aux données observées."),
          
          h5("Indicateurs d'ajustement basés sur les résidus"),
          tags$table(
            class = "table table-striped table-sm",
            tags$thead(
              tags$tr(
                tags$th("Indicateur"),
                tags$th("Valeur"),
                tags$th("Interprétation")
              )
            ),
            tags$tbody(
              tags$tr(
                tags$td("RMSR (Root Mean Square of Residuals)"),
                tags$td(round(residuals_data$rmsr, 3)),
                tags$td(
                  if(residuals_data$rmsr < 0.05) {
                    "Excellent (<0.05)"
                  } else if(residuals_data$rmsr < 0.08) {
                    "Bon (0.05-0.08)"
                  } else if(residuals_data$rmsr < 0.10) {
                    "Acceptable (0.08-0.10)"
                  } else {
                    "Médiocre (>0.10)"
                  }
                )
              ),
              tags$tr(
                tags$td("Résidu maximal (valeur absolue)"),
                tags$td(round(residuals_data$max_residual, 3)),
                tags$td(
                  if(residuals_data$max_residual < 0.10) {
                    "Excellent (<0.10)"
                  } else if(residuals_data$max_residual < 0.15) {
                    "Bon (0.10-0.15)"
                  } else {
                    "Préoccupant (>0.15)"
                  }
                )
              ),
              tags$tr(
                tags$td("Proportion de résidus >0.05 (valeur absolue)"),
                tags$td(paste0(round(residuals_data$prop_large_residuals * 100, 1), "%")),
                tags$td(
                  if(residuals_data$prop_large_residuals < 0.10) {
                    "Excellent (<10%)"
                  } else if(residuals_data$prop_large_residuals < 0.25) {
                    "Acceptable (10-25%)"
                  } else {
                    "Problématique (>25%)"
                  }
                )
              )
            )
          ),
          
          # Conclusion sur les résidus
          if(residuals_data$rmsr < 0.08 && residuals_data$max_residual < 0.15 && residuals_data$prop_large_residuals < 0.25) {
            div(
              class = "alert alert-success",
              tags$i(class = "fas fa-check-circle", style = "margin-right: 10px;"),
              strong("Conclusion : "),
              "Le modèle à ", strong(report$nfactors), " facteurs reproduit bien la structure de corrélation des données observées."
            )
          } else if(residuals_data$rmsr < 0.10 && residuals_data$max_residual < 0.20 && residuals_data$prop_large_residuals < 0.40) {
            div(
              class = "alert alert-warning",
              tags$i(class = "fas fa-exclamation-triangle", style = "margin-right: 10px;"),
              strong("Conclusion : "),
              "Le modèle reproduit de façon acceptable la structure de corrélation, mais certaines relations entre variables ne sont pas parfaitement capturées."
            )
          } else {
            div(
              class = "alert alert-danger",
              tags$i(class = "fas fa-times-circle", style = "margin-right: 10px;"),
              strong("Conclusion : "),
              "Le modèle ne reproduit pas adéquatement la structure de corrélation. Considérez d'augmenter le nombre de facteurs ou de modifier la spécification du modèle."
            )
          },
          
          # Visualisation optionnelle - description textuelle des résidus importants
          h5("Résidus importants"),
          if(residuals_data$max_residual > 0.10) {
            # Identifier les paires de variables avec des résidus importants
            residuals_mat <- residuals_data$residuals
            var_names <- rownames(residuals_mat)
            if(is.null(var_names)) {
              var_names <- paste0("Variable ", 1:nrow(residuals_mat))
            }
            
            # Trouver les résidus élevés (>0.10 en valeur absolue)
            high_residuals <- which(abs(residuals_mat) > 0.10 & lower.tri(residuals_mat), arr.ind = TRUE)
            
            if(nrow(high_residuals) > 0) {
              # Créer une table des résidus importants
              residuals_table <- data.frame(
                Variable1 = var_names[high_residuals[, 1]],
                Variable2 = var_names[high_residuals[, 2]],
                Residual = residuals_mat[high_residuals]
              )
              
              # Trier par valeur absolue des résidus
              residuals_table <- residuals_table[order(-abs(residuals_table$Residual)), ]
              
              # Limiter à 10 résidus maximum pour la lisibilité
              if(nrow(residuals_table) > 10) {
                residuals_table <- residuals_table[1:10, ]
              }
              
              # Créer un tableau HTML
              tags$table(
                class = "table table-striped table-sm",
                tags$thead(
                  tags$tr(
                    tags$th("Variable 1"),
                    tags$th("Variable 2"),
                    tags$th("Résidu")
                  )
                ),
                tags$tbody(
                  lapply(1:nrow(residuals_table), function(i) {
                    tags$tr(
                      tags$td(residuals_table$Variable1[i]),
                      tags$td(residuals_table$Variable2[i]),
                      tags$td(round(residuals_table$Residual[i], 3),
                              style = if(abs(residuals_table$Residual[i]) > 0.15) {
                                "color: #dc3545; font-weight: bold;"
                              } else {
                                "color: #495057;"
                              })
                    )
                  })
                )
              )
            } else {
              p("Aucun résidu important (>0.10) n'a été identifié.")
            }
          } else {
            p("Tous les résidus sont inférieurs à 0.10 en valeur absolue, indiquant une bonne reproduction des corrélations observées.")
          },
          
          # Interprétation et recommandations
          h5("Implications pour le modèle"),
          p("L'analyse des résidus suggère que :", 
            if(residuals_data$rmsr < 0.08) {
              "le nombre de facteurs retenus semble approprié pour représenter la structure des données."
            } else if(residuals_data$prop_large_residuals > 0.25) {
              "certaines relations entre variables ne sont pas bien expliquées par le modèle actuel. Un modèle avec plus de facteurs pourrait être envisagé."
            } else {
              "le modèle pourrait être amélioré, mais offre une représentation acceptable de la structure des données."
            }
          ),
          p("Les résidus élevés indiquent des relations entre variables qui ne sont pas capturées par le modèle factoriel. Ces relations pourraient suggérer des facteurs spécifiques, des relations dyadiques entre variables, ou des problèmes de mesure.")
        )
      )
    } else NULL,
    
    # Conclusion générale (toujours incluse)
    div(
      class = "card",
      div(
        class = "card-header",
        h4("Conclusion générale")
      ),
      div(
        class = "card-body",
        p("L'analyse factorielle exploratoire a identifié une structure à", strong(report$nfactors), "facteurs."),
        
        p("Les résultats de cette analyse peuvent être utilisés pour :", 
          tags$ul(
            tags$li("Développer ou raffiner des instruments de mesure"),
            tags$li("Explorer la structure sous-jacente des données"),
            tags$li("Guider des analyses confirmatoires ultérieures")
          )
        ),
        
        h5("Recommandations"),
        tags$ul(
          tags$li("Valider cette structure factorielle sur un échantillon indépendant"),
          tags$li("Comparer cette solution à des modèles théoriques existants dans la littérature"),
          tags$li("Envisager une analyse factorielle confirmatoire pour tester plus rigoureusement cette structure")
        )
      )
    )
  )
})

# Observateur pour le bouton de calcul des résidus
observeEvent(input$calculate_residuals, {
  req(efa_result())
  withProgress(message = 'Calcul des résidus...', {
    # Forcer le calcul des résidus
    res <- residual_matrix()
    if (!is.null(res)) {
      showNotification("Résidus calculés avec succès !", type = "message")
    } else {
      showNotification("Erreur lors du calcul des résidus. Vérifiez l'analyse factorielle.", type = "error")
    }
  })
})

# Affichage des statistiques des résidus
output$residual_stats <- renderUI({
  req(residual_matrix())
  
  res_data <- residual_matrix()
  
  tags$table(
    class = "table table-striped table-sm",
    tags$thead(
      tags$tr(
        tags$th("Indicateur"),
        tags$th("Valeur"),
        tags$th("Interprétation")
      )
    ),
    tags$tbody(
      tags$tr(
        tags$td("RMSR (Root Mean Square of Residuals)"),
        tags$td(round(res_data$rmsr, 3)),
        tags$td(
          if(res_data$rmsr < 0.05) {
            span(style = "color: #28a745;", "Excellent (<0.05)")
          } else if(res_data$rmsr < 0.08) {
            span(style = "color: #28a745;", "Bon (0.05-0.08)")
          } else if(res_data$rmsr < 0.10) {
            span(style = "color: #ffc107;", "Acceptable (0.08-0.10)")
          } else {
            span(style = "color: #dc3545;", "Médiocre (>0.10)")
          }
        )
      ),
      tags$tr(
        tags$td("Résidu maximal (valeur absolue)"),
        tags$td(round(res_data$max_residual, 3)),
        tags$td(
          if(res_data$max_residual < 0.10) {
            span(style = "color: #28a745;", "Excellent (<0.10)")
          } else if(res_data$max_residual < 0.15) {
            span(style = "color: #28a745;", "Bon (0.10-0.15)")
          } else {
            span(style = "color: #dc3545;", "Préoccupant (>0.15)")
          }
        )
      ),
      tags$tr(
        tags$td("Proportion de résidus >0.05 (valeur absolue)"),
        tags$td(paste0(round(res_data$prop_large_residuals * 100, 1), "%")),
        tags$td(
          if(res_data$prop_large_residuals < 0.10) {
            span(style = "color: #28a745;", "Excellent (<10%)")
          } else if(res_data$prop_large_residuals < 0.25) {
            span(style = "color: #28a745;", "Acceptable (10-25%)")
          } else {
            span(style = "color: #dc3545;", "Problématique (>25%)")
          }
        )
      )
    )
  )
})

# Interprétation des résidus
output$residual_interpretation <- renderUI({
  req(residual_matrix())
  
  res_data <- residual_matrix()
  
  # Déterminer la qualité globale de l'ajustement
  if(res_data$rmsr < 0.08 && res_data$max_residual < 0.15 && res_data$prop_large_residuals < 0.25) {
    div(
      class = "alert alert-success",
      tags$i(class = "fas fa-check-circle", style = "margin-right: 10px;"),
      strong("Bon ajustement : "),
      "Le modèle reproduit bien la structure de corrélation des données."
    )
  } else if(res_data$rmsr < 0.10 && res_data$max_residual < 0.20 && res_data$prop_large_residuals < 0.40) {
    div(
      class = "alert alert-warning",
      tags$i(class = "fas fa-exclamation-triangle", style = "margin-right: 10px;"),
      strong("Ajustement acceptable : "),
      "Le modèle reproduit de façon acceptable la structure de corrélation, mais certaines relations entre variables ne sont pas parfaitement capturées."
    )
  } else {
    div(
      class = "alert alert-danger",
      tags$i(class = "fas fa-times-circle", style = "margin-right: 10px;"),
      strong("Ajustement médiocre : "),
      "Le modèle ne reproduit pas adéquatement la structure de corrélation. Considérez d'augmenter le nombre de facteurs ou de modifier la spécification du modèle."
    )
  }
})

# Affichage de la matrice des résidus avec unicités sur la diagonale et légende des couleurs
output$residual_matrix_table <- renderDT({
  req(residual_matrix())
  req(efa_result())
  
  res_data <- residual_matrix()
  
  # Vérifier que res_data n'est pas NULL
  if (is.null(res_data)) {
    return(datatable(
      data.frame(Message = "Impossible de calculer les résidus. Assurez-vous que l'analyse factorielle est complète."),
      options = list(dom = 't')
    ))
  }
  
  # Récupérer la matrice des résidus
  residuals_mat <- res_data$residuals
  
  # Créer une matrice triangulaire supérieure (les valeurs en dessous de la diagonale sont NA)
  triangular_mat <- residuals_mat
  triangular_mat[lower.tri(triangular_mat, diag = FALSE)] <- NA  # On garde la diagonale cette fois
  
  # Arrondir les valeurs pour l'affichage
  triangular_mat <- round(triangular_mat, 3)
  
  # Convertir en dataframe pour DT
  df <- as.data.frame(triangular_mat)
  
  # Créer la table avec formatage conditionnel et options améliorées
  dt <- datatable(
    df,
    rownames = TRUE,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      autoWidth = TRUE,
      dom = 'tip',
      compact = TRUE,  # Style compact
      columnDefs = list(
        # Augmenter la largeur de la première colonne (rownames)
        list(width = '60px', targets = 0),  # La colonne d'index 0 est la colonne des rownames
        list(className = 'dt-center', targets = '_all')
      )
    ),
    caption = HTML("Matrice des résidus (triangulaire supérieure) - Les valeurs sur la diagonale représentent les unicités (1 - communalité)<br>
                   <div style='margin-top: 8px;'>
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #d73027; margin-right: 5px;'></span> < -0.15
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #fc8d59; margin-right: 5px; margin-left: 10px;'></span> -0.15 à -0.10
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #fee090; margin-right: 5px; margin-left: 10px;'></span> -0.10 à -0.05
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #ffffbf; margin-right: 5px; margin-left: 10px;'></span> -0.05 à 0.05
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #e0f3f8; margin-right: 5px; margin-left: 10px;'></span> 0.05 à 0.10
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #91bfdb; margin-right: 5px; margin-left: 10px;'></span> 0.10 à 0.15
                     <span style='display: inline-block; width: 12px; height: 12px; background-color: #4575b4; margin-right: 5px; margin-left: 10px;'></span> > 0.15
                   </div>")
  ) %>%
    formatStyle(
      columns = names(df),
      backgroundColor = styleInterval(
        c(-0.15, -0.10, -0.05, 0.05, 0.10, 0.15),
        c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4")
      ),
      color = styleInterval(
        c(-0.15, -0.10, 0.10, 0.15),
        c("#ffffff", "#000000", "#000000", "#000000", "#ffffff")
      ),
      fontWeight = styleInterval(
        c(-0.15, -0.10, 0.10, 0.15),
        c("bold", "normal", "normal", "normal", "bold")
      )
    ) %>%
    # Formatter les cellules NA pour qu'elles soient vides
    formatStyle(
      columns = names(df),
      valueColumns = names(df),
      target = 'cell',
      backgroundColor = styleEqual(NA, '#f9f9f9'),  # Couleur de fond légèrement grisée pour les NA
      color = styleEqual(NA, '#f9f9f9')  # Rendre le texte invisible pour les NA
    )
  
  return(dt)
})

# Fonction réactive pour calculer la matrice des résidus
residual_matrix <- reactive({
  message("Calculating residual matrix...")
  req(efa_result())
  req(sim_data())  # S'assurer que les données sont disponibles
  
  efa_res <- efa_result()
  
  # Vérifier si les charges factorielles sont disponibles
  loadings <- NULL
  
  # Chercher les charges factorielles dans l'ordre de priorité
  if (!is.null(efa_res$pattern) && is.matrix(efa_res$pattern)) {
    loadings <- efa_res$pattern
    message("Utilisation de pattern")
  } else if (!is.null(efa_res$loadingsROT) && is.matrix(efa_res$loadingsROT)) {
    loadings <- efa_res$loadingsROT
    message("Utilisation de loadingsROT")
  } else if (!is.null(efa_res$loadingsNOROT) && is.matrix(efa_res$loadingsNOROT)) {
    loadings <- efa_res$loadingsNOROT
    message("Utilisation de loadingsNOROT")
  }
  
  # Si aucune matrice de charges n'est disponible, arrêter le calcul
  if (is.null(loadings)) {
    message("Aucune matrice de charges factorielles trouvée")
    return(NULL)
  }
  
  # Récupérer les données brutes
  response_data <- sim_data()$response
  
  # Vérifier que response_data est bien une matrice ou un data.frame
  if (!is.matrix(response_data) && !is.data.frame(response_data)) {
    message("Les données de réponse ne sont pas dans un format matriciel")
    return(NULL)
  }
  
  # Calculer la matrice de corrélation observée directement à partir des données brutes
  message("Calcul de la matrice de corrélation à partir des données brutes")
  observed_corr <- tryCatch({
    cor(response_data, use = "pairwise.complete.obs")
  }, error = function(e) {
    message("Erreur lors du calcul de la matrice de corrélation:", e$message)
    return(NULL)
  })
  
  # Vérifier si observed_corr est bien calculée
  if (is.null(observed_corr)) {
    message("Matrice de corrélation observée non disponible")
    return(NULL)
  }
  
  # S'assurer que loadings est bien une matrice
  loadings <- as.matrix(loadings)
  
  # Calculer la matrice reproduite selon la rotation
  tryCatch({
    # Vérifier si les dimensions sont compatibles
    if (nrow(loadings) != nrow(observed_corr) || nrow(loadings) != ncol(observed_corr)) {
      message("Dimensions incompatibles: loadings:", paste(dim(loadings), collapse="x"), 
              ", observed_corr:", paste(dim(observed_corr), collapse="x"))
      return(NULL)
    }
    
    # Calculer la matrice reproduite en fonction du type de rotation
    if (!is.null(efa_res$phi) && is.matrix(efa_res$phi) && ncol(efa_res$phi) == ncol(loadings)) {
      # Vérifier que la matrice phi n'est pas identique à la matrice identité
      is_identity <- all(diag(efa_res$phi) == 1) && 
        sum(abs(efa_res$phi - diag(ncol(efa_res$phi)))) < 1e-10
      
      if (!is_identity) {
        # Rotation oblique
        message("Calcul avec rotation oblique")
        reproduced_corr <- loadings %*% efa_res$phi %*% t(loadings)
      } else {
        # Matrice phi est l'identité (rotation orthogonale)
        message("Calcul avec rotation orthogonale (phi est identité)")
        reproduced_corr <- loadings %*% t(loadings)
      }
    } else {
      # Rotation orthogonale ou par défaut
      message("Calcul avec rotation orthogonale")
      reproduced_corr <- loadings %*% t(loadings)
    }
    
    # Diagonal de la matrice reproduite = communalités
    # CORRECTION POUR GÉRER LA MATRICE DE COMMUNALITÉS
    if (!is.null(efa_res$communalities)) {
      if (is.matrix(efa_res$communalities) && ncol(efa_res$communalities) > 1 && 
          "Extraction" %in% colnames(efa_res$communalities)) {
        # Si c'est une matrice avec une colonne "Extraction"
        diag(reproduced_corr) <- efa_res$communalities[, "Extraction"]
        message("Utilisation des communalités (colonne Extraction) pour la diagonale")
      } else if (is.numeric(efa_res$communalities) && length(efa_res$communalities) == nrow(loadings)) {
        # Si c'est un vecteur numérique de la bonne longueur
        diag(reproduced_corr) <- efa_res$communalities
        message("Utilisation des communalités pour la diagonale")
      } else {
        # Format non reconnu
        diag(reproduced_corr) <- 1
        message("Format des communalités non reconnu, diagonale fixée à 1")
      }
    } else if (!is.null(efa_res$h2) && is.numeric(efa_res$h2) && 
               length(efa_res$h2) == nrow(loadings)) {
      diag(reproduced_corr) <- efa_res$h2
      message("Utilisation de h2 pour la diagonale")
    } else {
      # Si pas de communalités disponibles, mettre 1 sur la diagonale
      diag(reproduced_corr) <- 1
      message("Diagonale fixée à 1 (pas de communalités disponibles)")
    }
    
    # Calculer la matrice des résidus
    residuals <- observed_corr - reproduced_corr
    
    # Calculer les indicateurs d'ajustement
    rmsr <- sqrt(mean(residuals[lower.tri(residuals)]^2))  # Root Mean Square Residual
    max_residual <- max(abs(residuals[lower.tri(residuals)]))
    prop_large_residuals <- mean(abs(residuals[lower.tri(residuals)]) > 0.05)  # Proportion de résidus > 0.05
    
    message("Calcul des résidus terminé avec succès. RMSR =", rmsr)
    
    # Stocker les résultats
    return(list(
      observed = observed_corr,
      reproduced = reproduced_corr,
      residuals = residuals,
      rmsr = rmsr,
      max_residual = max_residual,
      prop_large_residuals = prop_large_residuals
    ))
  }, error = function(e) {
    # En cas d'erreur dans les calculs matriciels
    message("Erreur dans le calcul des résidus: ", e$message)
    return(NULL)
  })
})

# Ajouter un indicateur de statut pour le calcul des résidus
output$residuals_status <- renderText({
  if (is.null(residual_matrix())) {
    "Cliquez sur 'Calculer les résidus' pour lancer l'analyse"
  } else {
    "Résidus calculés avec succès"
  }
})

# Modification de l'observateur pour le bouton de calcul des résidus
observeEvent(input$calculate_residuals, {
  req(efa_result())
  showNotification("Calcul des résidus en cours...", type = "message", duration = NULL, id = "residuals_calc")
  
  # Forcer l'invalidation du cache
  invalidateLater(100)
  
  withProgress(message = 'Calcul des résidus...', {
    # Forcer le calcul des résidus
    res <- isolate({
      # Utiliser isolate pour éviter des boucles de calcul
      tryCatch({
        rm <- residual_matrix()
        if (is.null(rm)) "NULL" else "OK"
      }, error = function(e) {
        paste("Erreur:", e$message)
      })
    })
    
    # Afficher un message plus explicite basé sur le résultat
    if (res == "OK") {
      showNotification("Résidus calculés avec succès !", type = "message", id = "residuals_calc", duration = 4)
    } else {
      showNotification(paste("Problème lors du calcul des résidus:", res), 
                       type = "error", id = "residuals_calc", duration = 8)
    }
  })
})

}

shinyApp(ui, server)