library(shiny)
library(bslib)
library(DT)
library(plotly)
library(readxl)
library(writexl)
library(haven)
library(stats)
library(psych) 
library(corrplot)
library(lavaan) 
library(EFA.dimensions)
library(EFAfactors)
library(EFAtools)
library(EFAutilities)

# Fonction pour convertir un modèle EFA en syntaxe CFA
convert_efa_to_cfa <- function(efa_result, loading_threshold = 0.3) {
  # Extraire les charges factorielles 
  loadings_matrix <- efa_result$loadings
  
  # Identifier les facteurs
  factors <- colnames(loadings_matrix)
  
  # Initialiser la syntaxe CFA
  cfa_syntax <- c()
  
  # Parcourir chaque facteur et attribuer les variables indicatrices
  for (factor_idx in 1:length(factors)) {
    factor <- factors[factor_idx]
    
    # Obtenir les variables qui chargent sur ce facteur au-dessus du seuil
    factor_items <- rownames(loadings_matrix)[abs(loadings_matrix[, factor_idx]) >= loading_threshold]
    
    if (length(factor_items) > 0) {
      # Construire la syntaxe du modèle CFA pour ce facteur
      cfa_syntax <- c(cfa_syntax, paste(factor, "=~", paste(factor_items, collapse=" + ")))
    }
  }
  
  # Combiner en une seule chaîne
  cfa_model <- paste(cfa_syntax, collapse = "\n")
  
  return(cfa_model)
}

# Fonction pour générer des interprétations automatiques des facteurs
generate_factor_interpretation <- function(loadings, threshold = 0.3) {
  # Obtenir le nombre de facteurs
  n_factors <- ncol(loadings)
  
  # Structure pour stocker les interprétations
  interpretations <- list()
  
  # Pour chaque facteur
  for (f in 1:n_factors) {
    # Trouver les variables qui chargent significativement sur ce facteur
    significant_loadings <- loadings[, f]
    
    # Trier les charges par ordre décroissant de valeur absolue
    sorted_idx <- order(abs(significant_loadings), decreasing = TRUE)
    sorted_loadings <- significant_loadings[sorted_idx]
    var_names <- rownames(loadings)[sorted_idx]
    
    # Filtrer les charges au-dessus du seuil
    above_threshold <- abs(sorted_loadings) >= threshold
    strong_vars <- var_names[above_threshold]
    strong_loadings <- sorted_loadings[above_threshold]
    
    # Séparer les charges positives et négatives
    positive_vars <- strong_vars[strong_loadings > 0]
    negative_vars <- strong_vars[strong_loadings < 0]
    
    # Générer l'interprétation pour ce facteur
    factor_name <- paste0("Facteur ", f)
    
    # Analyser la force et le nombre de charges
    n_significant <- sum(above_threshold)
    avg_loading <- mean(abs(strong_loadings))
    
    # Évaluer la clarté du facteur
    if (n_significant == 0) {
      clarity <- "Ce facteur n'a pas de charges significatives au seuil défini."
      factor_quality <- "insuffisant"
    } else if (n_significant < 3) {
      clarity <- paste0("Ce facteur a seulement ", n_significant, " variable(s) avec des charges significatives, ce qui est généralement insuffisant pour une interprétation robuste.")
      factor_quality <- "faible"
    } else if (avg_loading < 0.4) {
      clarity <- "Ce facteur a des charges modérées, rendant son interprétation moins claire."
      factor_quality <- "modéré"
    } else if (avg_loading >= 0.6) {
      clarity <- "Ce facteur a des charges fortes et clairement définies."
      factor_quality <- "excellent"
    } else {
      clarity <- "Ce facteur a des charges de force moyenne à bonne."
      factor_quality <- "bon"
    }
    
    # Résumé des variables principales
    if (length(positive_vars) > 0) {
      pos_summary <- paste0("Variables avec charges positives : ", paste(positive_vars, collapse=", "))
    } else {
      pos_summary <- "Aucune variable avec charge positive significative"
    }
    
    if (length(negative_vars) > 0) {
      neg_summary <- paste0("Variables avec charges négatives : ", paste(negative_vars, collapse=", "))
    } else {
      neg_summary <- "Aucune variable avec charge négative significative"
    }
    
    # Suggestion d'interprétation
    suggestion <- ""
    if (n_significant >= 3) {
      suggestion <- paste0("Suggestion : Examinez les ", min(5, n_significant), " variables ayant les charges les plus élevées pour identifier le concept latent commun. ",
                           "Envisagez de nommer ce facteur en fonction du contenu théorique ou sémantique partagé par ces variables.")
    }
    
    # Regrouper toutes les informations
    interpretations[[f]] <- list(
      name = factor_name,
      quality = factor_quality,
      clarity = clarity,
      positive_vars = positive_vars,
      negative_vars = negative_vars,
      pos_summary = pos_summary,
      neg_summary = neg_summary,
      suggestion = suggestion,
      n_vars = n_significant,
      mean_loading = avg_loading
    )
  }
  
  return(interpretations)
}

# Palette de couleurs centralisée de l'application
pal_efa <- list(
  # Distributions univariées (histogrammes, boxplots, QQ-plots)
  hist_fill        = "skyblue",
  box_fill         = "lightgreen",
  qq_points        = "darkblue",
  density_obs      = "red",
  density_theo     = "blue",

  # Graphique des éboulis / analyse parallèle
  eigen_real       = "blue",
  eigen_parallel   = "red",
  kaiser_line      = "darkgreen",
  suggested_line   = "orange",
  grid_line        = "lightgray",

  # Biplot des charges factorielles (statut des variables)
  var_active       = "blue",
  var_inactive     = "gray",
  axis_line        = "grey",

  # Variance expliquée
  bar_variance     = "steelblue",
  line_cumulative  = "red",

  # Tableau des charges factorielles
  loading_pos_bg   = "#a8d08d",
  loading_neg_bg   = "#e15759",
  loading_pos_text = "#118dff",
  loading_neg_text = "#e15759",
  loading_neutral  = "white",

  # Communalités
  comm_bar         = "steelblue",
  comm_bar_modal   = "lightsteelblue",

  # Échelle à 4 niveaux (fiabilité, stabilité bootstrap)
  scale4_poor      = "#ffcccc",
  scale4_medium    = "#ffff9c",
  scale4_good      = "#a8d08d",
  scale4_excellent = "#32CD32",

  # Indices de modification (CFA)
  mi_low           = "white",
  mi_medium        = "#ffff9c",
  mi_high          = "#ffcccc",

  # Bootstrap
  boot_bg          = "#f5f5f5",
  boot_marker      = "black",

  # Qualité d'interprétation automatique des facteurs
  quality = list(
    "excellent"   = "#28a745",
    "bon"         = "#5cb85c",
    "modéré"      = "#ffc107",
    "faible"      = "#dc3545",
    "insuffisant" = "#6c757d"
  )
)

# Bulle d'aide contextuelle : icône (i) à placer à côté d'un titre de carte,
# affichant une courte explication au survol.
info_tip <- function(text, placement = "right") {
  tooltip(
    icon("info-circle", style = "color:#6c757d; font-size:0.85em; margin-left:6px; cursor:help;"),
    text,
    placement = placement
  )
}

# Interface utilisateur
ui <- page_sidebar(
  title = "Analyse Factorielle Exploratoire",
  theme = bs_theme(bootswatch = "flatly"),
  
  # Barre latérale avec options d'importation et paramètres d'analyse
  sidebar = sidebar(
    width = "350px",
    
    # Section d'importation des données
    accordion(
      id = "sidebar_accordion",
      accordion_panel(
        "Importation des données",
        radioButtons("import_method", "Méthode d'importation :",
                     choices = c(
                       "Charger un fichier" = "file",
                       "Données simulées" = "example",
                       "Données classiques" = "classic" 
                     ),
                     selected = "file",
                     inline = TRUE),
        
        conditionalPanel(
          condition = "input.import_method == 'file'",
          div(
            fileInput("data_file", NULL,
                      accept = c(".csv", ".sav", ".dta", ".xlsx", ".xls"),
                      buttonLabel = "Parcourir...",
                      placeholder = "Aucun fichier sélectionné"),
            
            conditionalPanel(
              condition = "input.data_file != null && input.data_file.name.toLowerCase().endsWith('.csv')",
              layout_columns(
                col_widths = c(6, 6),
                checkboxInput("header", "Noms des variables en 1ère ligne", TRUE),
                radioButtons("sep", "Séparateur (CSV)",
                             choices = c(Virgule = ",",
                                         "Point-virgule" = ";",
                                         Tabulation = "\\t"),
                             selected = ",",
                             inline = TRUE)
              ),
              layout_columns(
                col_widths = c(12),
                radioButtons("dec", "Point décimal (CSV)",
                             choices = c(Point = ".", Virgule = ","),
                             selected = ".",
                             inline = TRUE)
              )
            ),
            
            div(
              class = "border rounded p-3 mt-3",
              style = "background-color: #f8f9fa; width: 100%; min-height: 90px;",
              uiOutput("file_summary_compact")
            )
          )
        ),
        conditionalPanel(
          condition = "input.import_method == 'example'",
          h4("Configuration des données simulées"),
          numericInput("sim_nobs", "Nombre d'observations", value = 200, min = 50, max = 1000),
          numericInput("sim_nfactors", "Nombre de facteurs latents", value = 2, min = 1, max = 5),
          numericInput("sim_nvars", "Nombre de variables par facteur", value = 3, min = 2, max = 10),
          sliderInput("sim_loading", "Force des charges principales", min = 0.5, max = 0.9, value = 0.55, step = 0.05),
          sliderInput("sim_cross_loading", "Force des charges secondaires", min = 0, max = 0.4, value = 0.25, step = 0.05),
          sliderInput("sim_error", "Bruit aléatoire", min = 0.1, max = 1.0, value = 0.75, step = 0.05),
          checkboxInput("sim_correlated", "Facteurs corrélés", value = TRUE),
          conditionalPanel(
            condition = "input.sim_correlated == true",
            sliderInput("sim_factor_cor", "Corrélation entre facteurs", min = 0.1, max = 0.9, value = 0.3, step = 0.05)
          ),
          layout_columns(
            col_widths = c(6),
            actionButton("generate_example", "Générer", 
                         icon = icon("random"),
                         class = "btn-success",
                         width = "100%")
          )
        ),
        conditionalPanel(
          condition = "input.import_method == 'classic'",
          div(
            h4("Sélection de jeux de données classiques"),
            selectInput("classic_dataset", "Choisir un jeu de données :",
                        choices = c(
                          "-- Sélectionner --" = "",
                          "Packages standards" = c(
                            "Holzinger & Swineford (1939)" = "HolzingerSwineford1939",
                            "BFI Personality (25 items)" = "bfi_dataset"
                          ),
                          "EFA.dimensions" = c(
                            "NEO-PI-R (Big Five) - Neurotiscism" = "neo_n",
                            "NEO-PI-R (Big Five) - Extraversion" = "neo_e",
                            "NEO-PI-R (Big Five) - Ouverture" = "neo_o",
                            "Intelligence Émotionnelle (MSCEIT)" = "msceit"
                          ),
                          "EFAfactors" = c(
                            "BFI-2 (15 facets)" = "bfi2",
                            "SAPA Personality Inventory" = "spi"
                          ),
                          "EFAtools" = c(
                            "Holzinger (1939) - 9 tests cognitifs" = "holzinger",
                            "WIT data (14 variables)" = "WIT" 
                          ),
                          "EFAutilities" = c(
                            "Thurstone (1947) - Box dataset" = "box_data",
                            "Jöreskog (2003) - Données WAIS" = "wais_data"
                          )
                        )),
            tags$div(
              id = "dataset_description",
              style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",
              uiOutput("dataset_info")
            ),
            layout_columns(
              col_widths = c(6),
              actionButton("load_classic", "Charger", 
                           icon = icon("database"),
                           class = "btn-primary",
                           width = "100%")
            )
          )
        )
      ),
      
      # Section de filtrage des données
      accordion_panel(
        "Filtrage des données",
        conditionalPanel(
          condition = "output.data_loaded",
          selectInput("filter_var", "Variable de filtre",
                      choices = NULL),
          uiOutput("filter_controls"),
          actionButton("apply_filter", "Appliquer le filtre", 
                       icon = icon("filter"),
                       class = "btn-info",
                       width = "100%"),
          actionButton("reset_filter", "Réinitialiser", 
                       icon = icon("undo"),
                       width = "100%")
        )
      ),
      
      # Nouvelle section pour les paramètres d'analyse
      accordion_panel(
        "Paramètres de l'analyse",
        conditionalPanel(
          condition = "output.data_loaded",
          selectInput("efa_variables", "Variables à inclure", 
                      choices = NULL, multiple = TRUE),
          selectInput("efa_method", 
                      label = tagList("Méthode d'extraction", info_tip(
                        "Détermine comment les facteurs sont estimés à partir des corrélations observées. Le maximum de vraisemblance (ML) est le plus courant et permet des tests d'ajustement statistiques."
                      )),
                      choices = c(
                        "Courantes" = c(
                          "Maximum de vraisemblance" = "ML",
                          "Maximum de vraisemblance robuste" = "MLR",
                          "Moindres carrés non pondérés" = "ULS"
                        ),
                        "Avancées" = c(
                          "Maximum de vraisemblance avec correction de moyenne" = "MLM",
                          "Maximum de vraisemblance avec correction de moyenne et variance" = "MLMV",
                          "Maximum de vraisemblance avec correction de moyenne, variance et kurtosis" = "MLMVS",
                          "Maximum de vraisemblance avec correction de Fang" = "MLF",
                          "Moindres carrés généralisés" = "GLS",
                          "Moindres carrés pondérés" = "WLS",
                          "Moindres carrés pondérés diagonalement" = "DWLS",
                          "Moindres carrés à pondération distributionnelle" = "DLS"
                        )
                      ),
                      selected = "ML"),

          layout_column_wrap(
            width = 1,
            card(
              card_header("Rotation dynamique"),
              card_body(
                selectInput("dynamic_rotation", 
                            label = tagList("Changer la méthode de rotation", info_tip(
                              "La rotation ne change pas la qualité globale du modèle, mais redistribue les charges entre facteurs pour faciliter leur interprétation. Orthogonale (ex. Varimax) : facteurs non corrélés. Oblique (ex. Oblimin, Promax) : facteurs corrélés, souvent plus réaliste en sciences sociales."
                            )),
                            choices = c(
                              "Recommandées" = c(
                                "Aucune" = "none",
                                "Varimax (orthogonale)" = "varimax",
                                "Oblimin (oblique)" = "oblimin",
                                "Promax (oblique)" = "promax",
                                "Quartimax (orthogonale)" = "quartimax"
                              ),
                              "Avancées" = c(
                                "Orthomax" = "orthomax",
                                "Quartimin" = "quartimin",
                                "Geomin" = "geomin",
                                "Entropy" = "entropy",
                                "McCammon" = "mccammon",
                                "Infomax" = "infomax",
                                "Tandem1" = "tandem1",
                                "Tandem2" = "tandem2",
                                "Oblimax" = "oblimax",
                                "Bentler" = "bentler",
                                "Simplimax" = "simplimax",
                                "Target" = "target",
                                "PST" = "pst",
                                "CF" = "cf",
                                "Crawford-Ferguson" = "crawford-ferguson",
                                "CF-Quartimax" = "cf-quartimax",
                                "CF-Varimax" = "cf-varimax",
                                "CF-Equamax" = "cf-equamax",
                                "CF-Parsimax" = "cf-parsimax",
                                "CF-Facparsim" = "cf-facparsim"
                              )
                            ),
                            selected = "none"),
                actionButton("apply_rotation", "Appliquer la rotation", 
                             icon = icon("sync"),
                             class = "btn-primary",
                             width = "100%")
              )
            )
          ),
          numericInput("efa_nfactors", "Nombre de facteurs", 
                       value = 2, min = 1, max = 10),
          checkboxInput("efa_correlation", "Afficher la matrice de corrélation", TRUE),
          actionButton("run_efa", "Exécuter l'analyse", 
                       icon = icon("play"),
                       class = "btn-success",
                       width = "100%"),

        )
      )
    ),
    hr(),
    actionButton("help_button", "Aide", icon = icon("question-circle"), width = "100%")
  ),
  
  # Contenu principal avec système d'onglets
  navset_card_tab(
    title = "Analyse des données",
    id = "main_tabs",
    # Onglet d'introduction (ajouter au début des onglets principaux)
    nav_panel(
      title = "Introduction",
      icon = icon("info-circle"),
      
      layout_column_wrap(
        width = 1,
        card(
          card_header("Analyse Factorielle Exploratoire"),
          card_body(
            tags$h4("Bienvenue dans cette application d'analyse factorielle exploratoire"),
            tags$p("Cette application vous permet de réaliser des analyses factorielles exploratoires (AFE) 
               sur vos données de manière simple et interactive, tout en offrant des fonctionnalités avancées."),
            
            tags$h5("Fonctionnalités principales :"),
            tags$ul(
              tags$li(tags$b("Importation de données"), " - Chargez vos fichiers CSV, Excel, SPSS ou Stata, 
                  ou générez des données simulées pour tester l'application"),
              tags$li(tags$b("Tests préliminaires"), " - Vérifiez si vos données sont adaptées à l'AFE 
                  (KMO, Bartlett, déterminant)"),
              tags$li(tags$b("Analyse factorielle"), " - Plusieurs méthodes d'extraction et de rotation disponibles"),
              tags$li(tags$b("Visualisation"), " - Visualisez les résultats avec des graphiques interactifs"),
              tags$li(tags$b("Analyse confirmatoire"), " - Passez directement de l'AFE à l'AFC"),
              tags$li(tags$b("Bootstrap"), " - Évaluez la stabilité de votre solution factorielle")
            ),
            
            tags$h5("Comment utiliser cette application :"),
            tags$ol(
              tags$li("Commencez par importer vos données via la barre latérale à gauche"),
              tags$li("Explorez vos données dans l'onglet 'Aperçu des données'"),
              tags$li("Exécutez les tests préliminaires pour vérifier l'adéquation de vos données"),
              tags$li("Lancez l'analyse factorielle et interprétez les résultats"),
              tags$li("Utilisez les fonctionnalités de visualisation interactive pour affiner votre analyse")
            )
          )
        ),
        
        card(
          card_header("Ressources et références"),
          card_body(
            tags$h5("Principes de l'analyse factorielle exploratoire"),
            tags$p("L'AFE est une technique statistique permettant d'identifier des facteurs latents 
               (non directement observables) qui expliquent les corrélations entre variables observées."),
            
            tags$h5("Étapes clés de l'AFE :"),
            tags$ol(
              tags$li("Vérification de l'adéquation des données"),
              tags$li("Détermination du nombre de facteurs à extraire"),
              tags$li("Choix d'une méthode d'extraction des facteurs"),
              tags$li("Sélection d'une méthode de rotation"),
              tags$li("Interprétation des facteurs")
            ),
            
            tags$h5("Références utiles :"),
            tags$ul(
              tags$li(tags$a(href="https://www.jstor.org/stable/40247698", "Fabrigar, L. R., Wegener, D. T., MacCallum, R. C., & 
                         Strahan, E. J. (1999). Evaluating the use of exploratory factor analysis in psychological research. 
                         Psychological Methods, 4(3), 272-299.")),
              tags$li(tags$a(href="https://doi.org/10.1080/00273171.2017.1344538", "Schmitt, T. A. (2017). Current methodological 
                         considerations in exploratory and confirmatory factor analysis. 
                         Journal of Psychoeducational Assessment, 29(4), 304-321.")),
              tags$li(tags$a(href="https://doi.org/10.1177/0013164417719308", "Goretzko, D., Pham, T. T. H., & Bühner, M. (2021).
                         Exploratory factor analysis: Current use, methodological developments and recommendations for good practice.
                         Current Psychology, 40, 3510-3521."))
            ),
            
            tags$div(
              style = "margin-top: 20px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
              tags$p(style = "font-style: italic;", 
                     "Si vous rencontrez des difficultés dans l'utilisation de cette application, 
                 cliquez sur le bouton 'Aide' dans la barre latérale pour obtenir des informations supplémentaires.")
            )
          )
        )
      )
    ),
    # Onglet de visualisation des données
    nav_panel(
      title = "Aperçu des données",
      icon = icon("table"),
      
      # Utiliser des sous-onglets pour organiser l'aperçu des données
      navset_card_tab(
        id = "data_overview_tabs",
        
        # Sous-onglet pour l'affichage des données
        nav_panel(
          title = "Données",
          icon = icon("table"),
          
          layout_column_wrap(
            width = 1,
            card(
              card_header("Aperçu des données"),
              card_body(
                conditionalPanel(
                  condition = "output.data_loaded",
                  layout_columns(
                    col_widths = c(6, 6),
                    conditionalPanel(
                      condition = "input.import_method == 'example'",
                      div(
                        style = "display: flex; align-items: center; gap: 15px;",
                        downloadButton("download_simulated_data", "Télécharger", 
                                       icon = icon("download"),
                                       class = "btn-info"),
                        div(
                          style = "margin-top: 5px;",
                          radioButtons("download_format", "Format :",
                                       choices = c("CSV" = "csv",
                                                   "SPSS" = "sav",
                                                   "Stata" = "dta",
                                                   "Excel" = "xlsx"),
                                       selected = "csv",
                                       inline = TRUE)
                        )
                      )
                    )
                  ),
                  DTOutput("data_preview")
                )
              )
            )
          )
        ),
        
        # Sous-onglet pour les statistiques descriptives
        nav_panel(
          title = "Statistiques descriptives",
          icon = icon("chart-bar"),
          
          # Utiliser des sous-onglets pour organiser les statistiques descriptives
          navset_card_tab(
            id = "stats_tabs",
            
            # Sous-onglet pour la matrice de corrélation (optimisé pour utiliser tout l'espace)
            nav_panel(
              title = "Matrice de corrélation",
              icon = icon("table"),
              
              # Utiliser tout l'espace disponible pour la matrice
              card(
                card_header("Matrice de corrélation"),
                card_body(
                  radioButtons("corr_display_type", "Affichage :",
                               choices = c("Graphique" = "plot", "Texte" = "text"),
                               selected = "plot",
                               inline = TRUE),
                  uiOutput("correlation_output")
                )
              )
            ),
            
            # Sous-onglet pour le résumé statistique et tests de normalité (automatiques)
            nav_panel(
              title = "Résumé et normalité",
              icon = icon("info-circle"),
              
              layout_column_wrap(
                width = 1,
                card(
                  card_header("Résumé statistique"),
                  card_body(verbatimTextOutput("data_summary"))
                )
              ),
              
              layout_column_wrap(
                width = 1,
                card(
                  card_header("Tests de normalité"),
                  card_body(
                    DTOutput("normality_results_table")
                  )
                )
              )
            ),
            
            # Sous-onglet pour la matrice de visualisations
            nav_panel(
              title = "Visualisation",
              icon = icon("chart-column"),
              
              layout_column_wrap(
                width = 1,
                card(
                  card_header("Visualisation des variables"),
                  card_body(
                    fluidRow(
                      column(6,
                             selectInput("viz_variables", "Sélectionner des variables (max 8) :", 
                                         choices = NULL, multiple = TRUE, selected = NULL)
                      ),
                      column(3,
                             radioButtons("matrix_plot_type", "Type de graphique :",
                                          choices = c("Histogrammes" = "hist", 
                                                      "Boxplots" = "box",
                                                      "QQ-Plots" = "qq"),
                                          selected = "hist",
                                          inline = TRUE)
                      ),
                      column(3,
                             checkboxInput("matrix_density", "Ajouter courbes de densité", value = TRUE)
                      )
                    ),
                    plotOutput("matrix_plot", height = "600px")
                  )
                )
              )
            )
          )
        )
      )
    ),
    
    # Onglet pour l'analyse factorielle (caché par défaut)
    nav_panel(
      title = "Analyse factorielle",
      icon = icon("chart-line"),
      
      # Utiliser des sous-onglets pour organiser le contenu
      navset_card_tab(
        id = "efa_tabs",
        
        # Sous-onglet des tests préliminaires
        nav_panel(
          title = "Tests préliminaires",
          icon = icon("check-circle"),
          
          # Layout en deux colonnes, avec la gauche divisée en deux rangées
          layout_columns(
            col_widths = c(1/3, 2/3),
            
            # Colonne de gauche avec deux cartes empilées
            div(
              # Test de Bartlett en haut
              card(
                card_header("Test de Bartlett", info_tip(
                  "Teste si la matrice de corrélation diffère significativement d'une matrice identité (variables indépendantes). Un p < 0,05 indique que les variables sont assez corrélées entre elles pour justifier une analyse factorielle."
                )),
                card_body(
                  verbatimTextOutput("bartlett_test")
                ),
                height = "240px"
              ),
              # Déterminant de la matrice en bas
              card(
                card_header("Déterminant de la matrice", info_tip(
                  "Indique le degré de redondance entre les variables. Proche de 0 : certaines variables sont presque des combinaisons linéaires d'autres (multicolinéarité), ce qui peut déstabiliser l'estimation du modèle."
                )),
                card_body(
                  verbatimTextOutput("determinant_test")
                ),
                height = "240px"
              )
            ),
            
            # Colonne de droite avec mesure KMO
            card(
              card_header("Mesure de Kaiser-Meyer-Olkin (KMO)", info_tip(
                "Mesure la proportion de variance commune entre les variables par rapport au bruit. Plus le KMO est proche de 1, plus les corrélations observées ont des chances de refléter de véritables facteurs communs plutôt que du hasard."
              )),
              card_body(
                verbatimTextOutput("kmo_test")
              ),
              height = "480px"
            )
          ),
          
          # Graphique des éboulis en dessous
          layout_column_wrap(
            width = 1,
            card(
              card_header("Choix du nombre de facteurs", info_tip(
                "Chaque point est la variance expliquée (valeur propre) d'un facteur. La ligne rouge simule des données aléatoires (analyse parallèle) : ne retenez que les facteurs dont la valeur propre réelle dépasse la valeur simulée correspondante."
              )),
              card_body(
                plotOutput("scree_plot", height = "350px")
              )
            )
          )
        ),
        
        # Sous-onglet des résultats principaux
        nav_panel(
          title = "Résultats principaux",
          icon = icon("list"),
          
          layout_column_wrap(
            width = 1,
            card(
              card_header("Résultats de l'analyse factorielle"),
              card_body(
                verbatimTextOutput("efa_results")
              )
            )
          )
        ),
        
        # Sous-onglet de visualisation
        nav_panel(
          title = "Visualisation",
          icon = icon("chart-bar"),
          
          layout_columns(
            col_widths = c(3, 9),
            
            # Colonne de gauche avec tous les inputs regroupés dans une card
            card(
              card_header("Contrôles de visualisation"),
              card_body(
                selectInput("plot_factor_x", "Facteur en X:", 
                            choices = NULL,
                            selected = 1,
                            width = "150px"),
                selectInput("plot_factor_y", "Facteur en Y:", 
                            choices = NULL,
                            selected = 2,
                            width = "150px"),
                checkboxInput("auto_recalculate", "Recalcul auto", value = FALSE),
                actionButton("manual_recalculate", "Recalcul man", 
                             icon = icon("refresh"), 
                             class = "btn-info",
                             width = "50%"),
                selectizeInput("disabled_variables", "Variables désactivées :", 
                               choices = NULL, multiple = TRUE, 
                               options = list(plugins = list("remove_button"))),
                # Ajout d'une boîte pour les informations après recalcul
                tags$hr(),
                tags$h5("Informations sur l'analyse recalculée"),
                htmlOutput("recalculated_info")
              )
            ),
            
            # Colonne de droite avec le graphique
            card(
              card_header("Visualisation des facteurs"),
              card_body(
                plotlyOutput("factor_plot", height = "500px")
              )
            )
          )
        ),
        
        # Sous-onglet de comparaison des modèles
        nav_panel(
          title = "Comparaison des modèles",
          icon = icon("balance-scale"),
          
          layout_column_wrap(
            width = 1,
            card(
              card_header("Comparaison des modèles selon le nombre de facteurs", info_tip(
                "CFI/TLI proches de 1 et RMSEA/SRMR proches de 0 indiquent un bon ajustement. AIC/BIC servent à comparer les modèles entre eux : le plus petit est préférable."
              )),
              card_body(
                verbatimTextOutput("model_comparison_output")
              )
            )
          )
        ),
        # Sous-onglet pour l'analyse de fiabilité
        nav_panel(
          title = "Fiabilité des facteurs",
          icon = icon("chart-line"),
          
          layout_columns(
            col_widths = c(6, 6),  # Répartition 50-50 sur une ligne
            
            # Card 1 - Tableau des indices de fiabilité
            card(
              card_header("Indices de fiabilité", info_tip(
                "L'alpha de Cronbach et l'oméga de McDonald mesurent la cohérence interne des variables associées à un facteur — dans quelle mesure elles semblent mesurer la même chose. L'oméga est généralement préféré car il ne suppose pas que toutes les charges soient égales."
              )),
              card_body(
                selectInput("reliability_method", "Méthodes de fiabilité à afficher :",
                            choices = c("Alpha de Cronbach" = "alpha",
                                        "Oméga de McDonald" = "omega",
                                        "Les deux" = "both"),
                            selected = "both"),
                tags$div(
                  style = "font-size: 0.85em; color: #6c757d; margin-bottom: 8px;",
                  tags$span(style = "color: #ffcccc; font-weight: bold;", "■ "), "< 0,6 insuffisant   ",
                  tags$span(style = "color: #ffff9c; font-weight: bold;", "■ "), "0,6-0,7 discutable   ",
                  tags$span(style = "color: #a8d08d; font-weight: bold;", "■ "), "0,7-0,8 acceptable   ",
                  tags$span(style = "color: #32CD32; font-weight: bold;", "■ "), "≥ 0,8 bon"
                ),
                DTOutput("reliability_table")
              )
            ),
            
            # Card 2 - Analyse détaillée de fiabilité
            card(
              card_header("Analyse détaillée de fiabilité"),
              card_body(
                selectInput("factor_for_analysis", "Sélectionner un facteur :", 
                            choices = NULL),
                verbatimTextOutput("detailed_reliability")
              )
            )
          )
        )
      )
    ),

    
    # Onglet pour l'interprétation des résultats (caché par défaut)
    nav_panel(
      title = "Interprétation",
      icon = icon("search"),
      
      layout_column_wrap(
        width = 1,
        card(
          card_header("Tableaux des charges factorielles", info_tip(
            "Une charge indique la force du lien entre une variable et un facteur (de -1 à 1). Au-delà du seuil choisi, la variable est considérée comme représentative de ce facteur. Repère usuel : ≥ 0,30 acceptable, ≥ 0,50 solide."
          )),
          card_body(
            sliderInput("loadings_threshold", "Seuil de charge factorielle", 
                        min = 0, max = 1, value = 0.3, step = 0.05),
            tags$div(
              style = "font-size: 0.85em; color: #6c757d; margin-bottom: 8px;",
              tags$span(style = "color: #118dff; font-weight: bold;", "■ "), "charge positive au-delà du seuil   ",
              tags$span(style = "color: #e15759; font-weight: bold;", "■ "), "charge négative au-delà du seuil   ",
              tags$span(style = "color: black;", "■ "), "sous le seuil (peu informative pour ce facteur)"
            ),
            DTOutput("loadings_table")
          )
        )
      ),
      
      layout_column_wrap(
        width = 1/2,
        card(
          card_header("Variance expliquée", info_tip(
            "Proportion de la variance totale des variables captée par chaque facteur. La courbe cumulée indique combien de facteurs sont nécessaires pour expliquer une part satisfaisante de l'information (souvent visé : 50-60 % et plus)."
          )),
          card_body(
            plotlyOutput("variance_plot")
          )
        ),
        card(
          card_header("Communalités", info_tip(
            "Part de la variance d'une variable expliquée par l'ensemble des facteurs retenus (communalité + unicité = 1). Une communalité faible (< 0,30) signifie que la variable est mal représentée par la solution factorielle."
          )),
          card_body(
            DTOutput("communalities_table")
          )
        ),
        layout_column_wrap(
          width = 1,
          card(
            card_header("Suggestion d'interprétation automatique"),
            card_body(
              p("Ces suggestions sont générées automatiquement à partir des charges factorielles et peuvent vous aider à interpréter les facteurs."),
              uiOutput("auto_interpretation")
            )
          )
        )
      )
    ),
    
# Analyse bootstrap
nav_panel(
  title = "Bootstrap",
  icon = icon("random"),
  
  layout_columns(
    col_widths = c(7, 5),
    
    # Tableau des intervalles de confiance pour les charges
    card(
      card_header("Intervalles de confiance des charges factorielles"),
      card_body(
        layout_columns(
          col_widths = c(6, 6),
          numericInput("bootstrap_samples", "Nombre d'échantillons", 
                       value = 100, min = 50, max = 1000, step = 50),
          numericInput("bootstrap_ci", "Niveau de confiance (%)", 
                       value = 95, min = 80, max = 99, step = 1)
        ),
        actionButton("run_bootstrap", "Exécuter l'analyse bootstrap", 
                     icon = icon("play"),
                     class = "btn-success",
                     width = "100%"),
        tags$div(
          class = "mt-3",
          style = "font-size: 0.9em; color: #6c757d;",
          "Le bootstrap permet d'estimer la stabilité des charges factorielles en répliquant l'analyse sur de nouveaux échantillons créés par tirage avec remise."
        ),
        tags$hr(),
        selectInput("bootstrap_display_factor", "Afficher le facteur :", choices = NULL),
        DTOutput("bootstrap_loadings_table")
      )
    ),
    
    # Visualisation de la stabilité
    card(
      card_header("Stabilité des charges", info_tip(
        "Chaque violon montre la distribution des charges obtenues sur les échantillons bootstrap pour la variable sélectionnée. L'étoile noire est la charge originale ; plus le violon est étroit, plus l'estimation est stable."
      )),
      card_body(
        selectInput("bootstrap_stability_var", "Sélectionner une variable :", choices = NULL),
        plotlyOutput("bootstrap_stability_plot")
      )
    )
  )
),

    # Ajouter cet onglet après l'onglet "Interprétation"
    nav_panel(
      title = "Analyse Confirmatoire",
      icon = icon("check-circle"),
      
      # Utiliser layout_columns au lieu de layout_column_wrap pour plus de flexibilité
      layout_columns(
        col_widths = c(4, 8), # 1/3 - 2/3 ratio approximatif
        
        # Colonne de gauche: génération du modèle et syntaxe superposées verticalement
        card(
          card_header("Configurations CFA"),
          card_body(
            # Contrôles pour générer le modèle CFA
            card(
              card_header("Génération du modèle CFA"),
              card_body(
                sliderInput("cfa_threshold", "Seuil de charge factorielle", 
                            min = 0, max = 1, value = 0.3, step = 0.05),
                actionButton("generate_cfa", "Générer modèle CFA", 
                             icon = icon("cogs"), 
                             class = "btn-primary", 
                             width = "100%"),
                tags$p(class = "mt-2", style = "font-size: 0.9em; color: #6c757d;", 
                       "Cette fonction génère automatiquement un modèle CFA à partir des résultats de l'analyse exploratoire.")
              )
            ),
            
            # Syntaxe du modèle CFA
            card(
              card_header("Syntaxe du modèle CFA"),
              card_body(
                selectInput("cfa_method", "Méthode d'estimation",
                            choices = c("Maximum de vraisemblance" = "ML", 
                                        "Maximum de vraisemblance robuste" = "MLR",
                                        "Maximum de vraisemblance avec correction de moyenne" = "MLM",
                                        "Maximum de vraisemblance avec correction de moyenne et variance" = "MLMV",
                                        "Moindres carrés non pondérés" = "ULS",
                                        "Moindres carrés pondérés diagonalement" = "DWLS",
                                        "Moindres carrés pondérés" = "WLS"),
                            selected = "MLR"),
                textAreaInput("cfa_syntax", NULL, 
                              rows = 8, 
                              width = "100%", 
                              placeholder = "Le modèle CFA apparaîtra ici après génération..."),
                actionButton("run_cfa", "Exécuter l'analyse confirmatoire", 
                             icon = icon("play"), 
                             class = "btn-success", 
                             width = "100%")
              )
            )
          )
        ),
        
        # Colonne de droite: résultats de l'analyse
        card(
          card_header("Résultats de l'analyse confirmatoire", info_tip(
            "Indices d'ajustement : CFI/TLI proches de 1 et RMSEA/SRMR proches de 0 indiquent un bon ajustement. Diagramme : les flèches sont les charges standardisées, les cercles les facteurs latents (non observés), les rectangles les variables mesurées."
          )),
          card_body(
            navset_card_tab(
              nav_panel(
                title = "Résumé du modèle",
                verbatimTextOutput("cfa_summary")
              ),
              nav_panel(
                title = "Indices d'ajustement",
                verbatimTextOutput("cfa_fit_indices")
              ),
              nav_panel(
                title = "Diagramme du modèle",
                plotOutput("cfa_diagram", height = "600px")
              ),
              nav_panel(
                title = "Modification d'indices",
                DTOutput("modification_indices")
              )
            )
          )
        )
      )
    )
  )
)

# Serveur
server <- function(input, output, session) {
  
  # Définir l'onglet Introduction comme onglet actif au démarrage
  updateTabsetPanel(session, "main_tabs", selected = "Introduction")
  
  # Réactives pour stocker les données
  raw_data <- reactiveVal(NULL)
  filtered_data <- reactiveVal(NULL)
  efa_result <- reactiveVal(NULL)
  model_comparison <- reactiveVal(NULL)
  unrotated_solution <- reactiveVal(NULL)
  
  # Fonction pour nettoyer les objets réactifs
  cleanup_reactive_objects <- function() {
    # Réinitialiser les objets réactifs principaux
    efa_result(NULL)
    model_comparison(NULL)
    unrotated_solution(NULL)
    bootstrap_results(NULL)
    cfa_result(NULL)
    normality_results_data(NULL)
    
    # Libérer la mémoire explicitement
    gc()
    
    # Message de débogage en console (pas visible par l'utilisateur)
    message("Nettoyage de la mémoire effectué")
  }
  
  # Indicateur de chargement des données
  output$data_loaded <- reactive({
    return(!is.null(raw_data()))
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  # Cacher les onglets d'analyse au démarrage
  observe({
    hideTab(inputId = "main_tabs", target = "Analyse factorielle")
    hideTab(inputId = "main_tabs", target = "Interprétation")
    hideTab(inputId = "main_tabs", target = "Bootstrap")  # Ajoutez cette ligne
    hideTab(inputId = "main_tabs", target = "Analyse Confirmatoire")
  })
  
  # Chargement automatique lorsqu'un fichier est sélectionné
  observeEvent(input$data_file, {
    req(input$data_file)
    
    # Nettoyer les objets existants avant de charger de nouvelles données
    cleanup_reactive_objects()
    
    # Afficher une notification pendant le chargement
    id <- showNotification("Chargement des données en cours...", type = "message", duration = NULL)
    
    # Obtenir l'extension et le chemin du fichier
    file_ext <- tools::file_ext(input$data_file$name)
    file_path <- input$data_file$datapath
    
    # Importer les données selon le format
    tryCatch({
      data <- NULL
      
      if (file_ext == "csv") {
        # S'assurer que les contrôles sont disponibles avant de lire le fichier
        req(input$header, input$sep, input$dec)
        
        # Afficher un message de débogage
        message("Chargement CSV avec: sep=", input$sep, ", dec=", input$dec, ", header=", input$header)
        
        data <- read.csv(file_path, 
                         header = input$header, 
                         sep = input$sep, 
                         dec = input$dec,
                         stringsAsFactors = FALSE)
      } else if (file_ext %in% c("xls", "xlsx")) {
        data <- as.data.frame(readxl::read_excel(file_path))
      } else if (file_ext == "sav") {
        data <- as.data.frame(haven::read_spss(file_path))
      } else if (file_ext == "dta") {
        data <- as.data.frame(haven::read_stata(file_path))
      } else {
        stop("Format de fichier non pris en charge")
      }
      
      # Correction des noms de colonnes vides
      if (any(names(data) == "")) {
        empty_names <- which(names(data) == "")
        new_names <- names(data)
        new_names[empty_names] <- paste0("Colonne_", empty_names)
        names(data) <- new_names
      }
      
      # Conversion des variables labellées
      if (file_ext %in% c("sav", "dta")) {
        for (col in names(data)) {
          if (inherits(data[[col]], "haven_labelled")) {
            if (is.numeric(haven::zap_labels(data[[col]]))) {
              data[[col]] <- as.numeric(haven::zap_labels(data[[col]]))
            } else {
              data[[col]] <- as.factor(haven::as_factor(data[[col]]))
            }
          }
        }
      }
      
      # Stocker les données
      raw_data(data)
      filtered_data(data)
      
      # Mettre à jour les sélecteurs
      updateSelectInput(session, "filter_var", 
                        choices = c("Aucun" = "", names(data)))
      
      numeric_cols <- names(data)[sapply(data, is.numeric)]
      updateSelectInput(session, "efa_variables", 
                        choices = numeric_cols,
                        selected = numeric_cols)
      
      # Afficher l'onglet d'analyse
      showTab(inputId = "main_tabs", target = "Analyse factorielle")
      
      # Fermer "Importation des données" et ouvrir "Paramètres de l'analyse"
      accordion_panel_set(id = "sidebar_accordion", values = "Paramètres de l'analyse")
      
      removeNotification(id)
      showNotification("Données importées avec succès", type = "message")
    })
  })
  
  # Créer une réactive pour stocker les résultats des tests de normalité
  # (alimentée automatiquement par l'observe() de la section "Aperçu des données")
  normality_results_data <- reactiveVal(NULL)
  
  observeEvent(input$apply_rotation, {
    req(unrotated_solution())
    
    solution <- isolate(unrotated_solution())
    
    withProgress(message = "Application de la nouvelle rotation...", {
      tryCatch({
        # Code pour appliquer la rotation
        # Vérifier d'abord si data_complete existe dans solution
        if (!("data_complete" %in% names(solution))) {
          # Si data_complete n'existe pas, utiliser les variables sélectionnées dans filtered_data
          data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
          data_complete <- na.omit(data_subset)
        } else {
          # Sinon, utiliser la data_complete déjà stockée
          data_complete <- solution$data_complete
        }
        
        new_efa_model <- lavaan::efa(
          data = data_complete,
          nfactors = solution$nfactors,
          rotation = input$dynamic_rotation,
          estimator = solution$method,
          control = list(
            optim.method = "BFGS",
            optim.control = list(maxit = 1000),
            check.gradient = FALSE
          )
        )
        model_summary <- summary(new_efa_model, nd = 3, cutoff = 0.001)
        
        if (inherits(new_efa_model, "efaList")) {
          loadings_matrix <- model_summary$loading
          
          if (is.null(loadings_matrix) || all(is.na(loadings_matrix))) {
            model <- new_efa_model[[1]]
            params <- lavaan::parameterEstimates(model)
            vars <- unique(params$rhs[params$op == "=~"])
            factors <- unique(params$lhs[params$op == "=~"])
            
            loadings_matrix <- matrix(0, nrow = length(vars), ncol = length(factors))
            rownames(loadings_matrix) <- vars
            colnames(loadings_matrix) <- factors
            
            for (i in seq_len(nrow(params))) {
              if (params$op[i] == "=~") {
                f_idx <- match(params$lhs[i], factors)
                v_idx <- match(params$rhs[i], vars)
                loadings_matrix[v_idx, f_idx] <- params$est[i]
              }
            }
          }
          
          colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
          
          communalities <- model_summary$communalities
          if (is.null(communalities)) {
            communalities <- rowSums(loadings_matrix^2)
          }
          
          fit_indices <- lavaan::fitMeasures(new_efa_model[[1]])
          
        } else {
          loadings_matrix <- lavaan::inspect(new_efa_model, "est")$lambda
          colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
          
          communalities <- rowSums(loadings_matrix^2)
          
          fit_indices <- lavaan::fitMeasures(new_efa_model)
        }
        
        uniqueness <- 1 - communalities
        names(communalities) <- rownames(loadings_matrix)
        names(uniqueness) <- rownames(loadings_matrix)
        
        variance_explained <- colSums(loadings_matrix^2) / nrow(loadings_matrix)
        cumulative_variance <- cumsum(variance_explained)
        
        # Créer un objet de résultat complet
        efa_result(list(
          efa = new_efa_model,
          loadings = loadings_matrix,
          uniquenesses = uniqueness,
          communality = communalities,
          method = solution$method,
          rotation = input$dynamic_rotation,
          nfactors = solution$nfactors,
          Vaccounted = data.frame(
            SS_loadings = colSums(loadings_matrix^2),
            Proportion = variance_explained,
            Cumulative = cumulative_variance
          ),
          fit_indices = fit_indices,
          summary = model_summary,  # Ajouter le résumé complet
          reliability = calculate_reliability(data_complete, loadings_matrix, input$loadings_threshold)
        ))
        
        showNotification(paste("Rotation", input$dynamic_rotation, "appliquée avec succès"), 
                         type = "message")
        
      }, error = function(e) {
        showNotification(
          paste("Erreur lors de l'application de la rotation:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
  })

  
  # Visualisation des facteurs
  output$factor_plot <- renderPlotly({
    req(efa_result())
    efa <- efa_result()
    
    # Définir les indices des facteurs à afficher
    factor_idx_x <- 1
    factor_idx_y <- 2
    
    # Si plus de 2 facteurs et des sélecteurs sont utilisés
    if (efa$nfactors > 2 && !is.null(input$plot_factor_x) && !is.null(input$plot_factor_y)) {
      factor_idx_x <- as.numeric(input$plot_factor_x)
      factor_idx_y <- as.numeric(input$plot_factor_y)
    }
    
    # Extraire les charges factorielles pour les facteurs sélectionnés
    loadings <- efa$loadings
    
    if (ncol(loadings) < 2) {
      # Si un seul facteur est extrait
      var_names <- rownames(loadings)
      factor1 <- loadings[, 1]
      factor2 <- rep(0, nrow(loadings))
    } else {
      var_names <- rownames(loadings)
      factor1 <- loadings[, factor_idx_x]
      factor2 <- loadings[, factor_idx_y]
    }
    
    # Vérifier s'il y a des variables désactivées
    disabled_vars <- input$disabled_variables
    if (is.null(disabled_vars)) {
      disabled_vars <- c()
    }
    
    # Créer un data frame pour le graphique
    plot_data <- data.frame(
      Variable = var_names,
      Factor1 = factor1,
      Factor2 = factor2,
      Status = ifelse(var_names %in% disabled_vars, "Désactivée", "Active")
    )
    
    # Obtenir les pourcentages de variance pour l'étiquetage
    variance_percent <- round(efa$Vaccounted$Proportion * 100, 1)
    
    # Créer le graphique
    p <- plot_ly(plot_data, x = ~Factor1, y = ~Factor2, 
                 text = ~Variable,
                 type = "scatter", 
                 mode = "markers+text", 
                 textposition = "top center",
                 color = ~Status,
                 colors = c("Active" = pal_efa$var_active, "Désactivée" = pal_efa$var_inactive),
                 marker = list(size = 10),
                 customdata = ~Variable) %>%
      layout(
        title = "Visualisation des charges factorielles<br><span style='font-size:12px'>(Cliquez sur un point pour activer/désactiver la variable)</span>",
        xaxis = list(title = paste0("Facteur ", factor_idx_x, " (", variance_percent[factor_idx_x], "% de variance)")),
        yaxis = list(title = paste0("Facteur ", factor_idx_y, " (", variance_percent[factor_idx_y], "% de variance)")),
        shapes = list(
          list(
            type = "line", line = list(color = pal_efa$axis_line, dash = "dash"), 
            x0 = 0, x1 = 0, y0 = -1, y1 = 1
          ),
          list(
            type = "line", line = list(color = pal_efa$axis_line, dash = "dash"), 
            y0 = 0, y1 = 0, x0 = -1, x1 = 1
          )
        ),
        legend = list(title = list(text = "Statut de la variable"))
      ) %>%
      event_register("plotly_click")
    
    return(p)
  })
      
      # Tests préliminaires
      output$bartlett_test <- renderPrint({
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        
        data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
        data_complete <- na.omit(data_subset)
        
        bartlett_res <- psych::cortest.bartlett(cor(data_complete), n = nrow(data_complete))
        
        cat("Test de sphéricité de Bartlett\n")
        cat("----------------------------\n")
        cat("Chi-carré:", round(bartlett_res$chisq, 2), "\n")
        cat("Degrés de liberté:", bartlett_res$df, "\n")
        cat("p-value:", format.pval(bartlett_res$p.value, digits = 4), "\n\n")
        
        if (bartlett_res$p.value < 0.05) {
          cat("Interprétation: Les variables sont corrélées (p < 0.05).\n")
          cat("L'analyse factorielle est appropriée.")
        } else {
          cat("Interprétation: Les variables ne sont pas suffisamment corrélées (p > 0.05).\n")
          cat("L'analyse factorielle peut ne pas être appropriée.")
        }
      })
      
      # Modifier la sortie du test KMO dans l'onglet "Tests préliminaires"
      output$kmo_test <- renderPrint({
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        
        data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
        data_complete <- na.omit(data_subset)
        
        kmo_res <- psych::KMO(cor(data_complete))
        
        cat("Mesure d'adéquation de Kaiser-Meyer-Olkin (KMO)\n")
        cat("-------------------------------------------\n")
        cat("KMO général:", round(kmo_res$MSA, 3), "\n\n")
        
        # Interprétation du KMO général
        if (kmo_res$MSA >= 0.8) {
          cat("Interprétation générale: Excellent (≥ 0.8)\n\n")
        } else if (kmo_res$MSA >= 0.7) {
          cat("Interprétation générale: Bon (≥ 0.7)\n\n")
        } else if (kmo_res$MSA >= 0.6) {
          cat("Interprétation générale: Moyen (≥ 0.6)\n\n")
        } else if (kmo_res$MSA >= 0.5) {
          cat("Interprétation générale: Médiocre (≥ 0.5)\n\n")
        } else {
          cat("Interprétation générale: Inacceptable (< 0.5)\n")
          cat("L'analyse factorielle n'est pas recommandée.\n\n")
        }
        
        # Ajout des MSA individuelles pour chaque variable
        cat("MSA pour chaque variable:\n")
        cat("-----------------------\n")
        
        # Extraire les MSA pour chaque variable
        individual_msa <- kmo_res$MSAi
        
        # Créer une table avec deux colonnes: Variable et MSA
        msa_table <- data.frame(
          Variable = names(individual_msa),
          MSA = round(individual_msa, 3)
        )
        
        # Ajouter une colonne d'interprétation
        msa_table$Interprétation <- cut(
          msa_table$MSA,
          breaks = c(0, 0.5, 0.6, 0.7, 0.8, 1),
          labels = c("Inacceptable", "Médiocre", "Moyen", "Bon", "Excellent"),
          include.lowest = TRUE
        )
        
        # Trier par ordre croissant de MSA pour identifier les variables problématiques
        msa_table <- msa_table[order(msa_table$MSA), ]
        
        # Afficher le tableau
        print(msa_table, row.names = FALSE)
        
        # Recommandation si certaines variables ont un MSA faible
        if (any(msa_table$MSA < 0.5)) {
          cat("\nRecommandation: Considérez retirer les variables avec MSA < 0.5 ")
          cat("avant de procéder à l'analyse factorielle.\n")
        }
      })
      
      # Scree Plot avec analyse parallèle
      output$scree_plot <- renderPlot({
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        
        withProgress(message = "Génération du graphique des éboulis...", {
          tryCatch({
            # Sélection des variables pour l'analyse
            data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
            data_complete <- na.omit(data_subset)
            
            # Vérifier qu'il y a suffisamment d'observations
            if (nrow(data_complete) < 3) {
              plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
              text(0, 0, "Données insuffisantes pour générer le graphique", cex = 1.2)
              return()
            }
            
            # Calculer la matrice de corrélation manuellement
            R <- cor(data_complete, use = "pairwise.complete.obs")
            
            # Calculer les valeurs propres manuellement
            eigen_values <- eigen(R)$values
            
            # Calculer l'analyse parallèle
            # Utiliser psych::fa.parallel pour calculer l'analyse parallèle
            incProgress(0.3, detail = "Calcul de l'analyse parallèle...")
            
            # Désactiver les sorties de fa.parallel
            op <- options(console.quiet = TRUE)
            on.exit(options(op))
            
            # Analyse parallèle
            fa_parallel <- psych::fa.parallel(data_complete, 
                                              fa = "both", # Calculer FA et PCA
                                              fm = "minres", # Méthode d'extraction
                                              n.iter = 100, # Nombre de simulations
                                              show.legend = FALSE,
                                              plot = FALSE) # Ne pas afficher le graphique
            
            # Extraire les valeurs simulées
            sim_eigen_values <- fa_parallel$fa.values
            sim_eigen_values <- sim_eigen_values[1:length(eigen_values)]
            
            # Définir les marges et initialiser un nouveau dispositif graphique
            par(mar = c(4, 4, 3, 2))
            
            # Définir les limites du graphique
            y_max <- max(max(eigen_values), max(sim_eigen_values)) * 1.1
            
            # Tracer les valeurs propres
            plot(eigen_values, type = "b", pch = 19, col = pal_efa$eigen_real,
                 xlab = "Nombre de facteurs", ylab = "Valeur propre",
                 main = "Graphique des éboulis avec analyse parallèle",
                 xlim = c(1, length(eigen_values)), 
                 ylim = c(0, y_max),
                 xaxt = "n")  # Supprimer l'axe des x par défaut
            
            # Ajouter l'analyse parallèle
            lines(sim_eigen_values, type = "b", pch = 17, col = pal_efa$eigen_parallel, lty = 2)
            
            # Ajouter un axe des x personnalisé avec les numéros de facteurs
            axis(1, at = 1:length(eigen_values), labels = 1:length(eigen_values))
            
            # Ajouter une ligne horizontale à la valeur 1 (critère de Kaiser)
            abline(h = 1, lty = 3, col = pal_efa$kaiser_line)
            
            # Trouver où les valeurs propres réelles croisent les valeurs simulées
            cross_point <- which(eigen_values < sim_eigen_values)[1]
            if (!is.na(cross_point) && !is.null(cross_point)) {
              suggested_factors <- cross_point - 1
              # Ajouter une ligne verticale au point de croisement suggéré
              abline(v = suggested_factors, lty = 3, col = pal_efa$suggested_line)
              # Ajouter du texte pour indiquer le nombre recommandé de facteurs
              text(suggested_factors, y_max * 0.9, paste("Recommandé:", suggested_factors), 
                   col = pal_efa$suggested_line, pos = 4, cex = 0.9)
            }
            
            # Ajouter une grille pour faciliter la lecture
            grid(nx = NULL, ny = NULL, lty = "dotted", col = pal_efa$grid_line)
            
            # Ajouter une légende
            legend("topright", 
                   legend = c("Valeurs propres réelles", 
                              "Analyse parallèle (simulations)", 
                              "Critère de Kaiser (VP = 1)"),
                   col = c(pal_efa$eigen_real, pal_efa$eigen_parallel, pal_efa$kaiser_line), 
                   lty = c(1, 2, 3), 
                   pch = c(19, 17, NA),
                   cex = 0.8,
                   bg = "white")
            
          }, error = function(e) {
            # En cas d'erreur, afficher un message
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, paste("Erreur:", e$message), cex = 1.2)
          })
        })
      })
      
      # EFA Results
      output$efa_results <- renderPrint({
        req(efa_result())
        efa <- efa_result()
        
        cat("ANALYSE FACTORIELLE EXPLORATOIRE\n")
        cat("================================\n\n")
        cat("Méthode d'extraction: ", efa$method, "\n")
        cat("Méthode de rotation: ", efa$rotation, "\n")
        cat("Nombre de facteurs: ", efa$nfactors, "\n\n")
        
        # Vérifier si l'objet efa est de classe efaList ou un modèle lavaan
        if (inherits(efa$efa, "efaList")) {
          cat("RÉSUMÉ DU MODÈLE\n")
          cat("--------------\n")
          print(summary(efa$efa, nd = 3, cutoff = 0.2, dot.cutoff = 0.05))
          
          cat("\nINDICES D'AJUSTEMENT DU MODÈLE\n")
          cat("----------------------------\n")
          # Utiliser fit.measures = "all" pour obtenir tous les indices d'ajustement
          fit_measures <- fitMeasures(efa$efa, fit.measures = "all")
          print(fit_measures)
        } else {
          # Pour les anciennes versions de lavaan
          cat("INDICES D'AJUSTEMENT DU MODÈLE\n")
          cat("----------------------------\n")
          # Utiliser fit.measures = "all" pour obtenir tous les indices d'ajustement
          fit_indices <- fitMeasures(efa$efa, fit.measures = "all")
          print(fit_indices)
          
          cat("\nMATRICE DES CHARGES FACTORIELLES\n")
          cat("------------------------------\n")
          loadings <- round(efa$loadings, 3)
          print(loadings)
        }
        
        cat("\nVARIANCE EXPLIQUÉE PAR FACTEUR\n")
        cat("---------------------------\n")
        vacc <- efa$Vaccounted
        vacc_rounded <- round(vacc * 100, 2)
        print(vacc_rounded)
      })
      outputOptions(output, "efa_results", suspendWhenHidden = FALSE)
      
      # Loadings Table
      output$loadings_table <- renderDT({
        req(efa_result())
        loadings <- efa_result()$loadings
        threshold <- input$loadings_threshold
        
        # Appliquer le seuil
        loadings_thresholded <- loadings
        loadings_thresholded[abs(loadings_thresholded) < threshold] <- NA
        
        # Formater pour l'affichage
        loadings_formatted <- round(loadings_thresholded, 3)
        
        # Créer un tableau avec coloration selon la valeur
        datatable(
          loadings_formatted,
          options = list(pageLength = 25, scrollX = TRUE),
          rownames = TRUE
        ) %>%
          formatStyle(
            columns = names(loadings_formatted),
            backgroundColor = styleInterval(
              c(-threshold, threshold),
              c(pal_efa$loading_neutral, pal_efa$loading_neutral, pal_efa$loading_pos_bg)
            ),
            color = styleInterval(
              c(-threshold, threshold),
              c(pal_efa$loading_neg_text, "black", pal_efa$loading_pos_text)
            )
          )
      })
      
      # Variance Plot
      output$variance_plot <- renderPlotly({
        req(efa_result())
        efa <- efa_result()
        
        # Préparer les données pour le graphique
        var_data <- data.frame(
          Factor = paste0("F", 1:efa$nfactors),
          Variance = efa$Vaccounted$Proportion * 100,
          Cumulative = efa$Vaccounted$Cumulative * 100
        )
        
        # Créer le graphique
        plot_ly() %>%
          add_bars(x = ~var_data$Factor, y = ~var_data$Variance, name = "Variance expliquée (%)",
                   marker = list(color = pal_efa$bar_variance)) %>%
          add_trace(x = ~var_data$Factor, y = ~var_data$Cumulative, name = "Cumulative (%)",
                    type = "scatter", mode = "lines+markers", yaxis = "y2",
                    line = list(color = pal_efa$line_cumulative, width = 2),
                    marker = list(color = pal_efa$line_cumulative, size = 8)) %>%
          layout(
            title = "Variance expliquée par facteur",
            xaxis = list(title = "Facteurs"),
            yaxis = list(title = "Variance expliquée (%)"),
            yaxis2 = list(
              title = "Variance cumulée (%)",
              overlaying = "y",
              side = "right",
              range = c(0, 100)
            ),
            legend = list(x = 0.6, y = 0.9),
            barmode = "group"
          )
      })
      
      # Table des communalités
      output$communalities_table <- renderDT({
        req(efa_result())
        comm <- efa_result()$communality
        uniq <- efa_result()$uniquenesses
        
        # Préparer les données pour le tableau
        comm_data <- data.frame(
          Variable = names(comm),
          Communalité = round(comm, 3),
          Unicité = round(uniq, 3)
        )
        
        # Créer le tableau
        datatable(
          comm_data,
          options = list(pageLength = 25, scrollX = TRUE),
          rownames = FALSE
        ) %>%
          formatStyle(
            "Communalité",
            background = styleColorBar(c(0, 1), pal_efa$comm_bar),
            backgroundSize = "100% 90%",
            backgroundRepeat = "no-repeat",
            backgroundPosition = "center"
          )
      })
      
      # Résumé statistique
      output$data_summary <- renderPrint({
        req(filtered_data())
        summary(filtered_data())
      })
      outputOptions(output, "data_summary", suspendWhenHidden = FALSE)
      
      # Application du filtre
      observeEvent(input$apply_filter, {
        req(input$filter_var, raw_data())
        
        if (input$filter_var == "") {
          filtered_data(raw_data())
          return()
        }
        
        data <- raw_data()
        var <- data[[input$filter_var]]
        
        if (is.numeric(var) || is.integer(var)) {
          req(input$filter_range)
          filtered <- data[var >= input$filter_range[1] & var <= input$filter_range[2], ]
        } else if (is.factor(var) || is.character(var)) {
          req(input$filter_categories)
          filtered <- data[var %in% input$filter_categories, ]
        } else {
          filtered <- data
        }
        
        filtered_data(filtered)
        showNotification(paste0(nrow(filtered), " observations après filtrage"), type = "message")
      })
      
      # Réinitialisation du filtre
      observeEvent(input$reset_filter, {
        filtered_data(raw_data())
        updateSelectInput(session, "filter_var", selected = "")
        showNotification("Filtres réinitialisés", type = "message")
      })
      
      # Aide
      observeEvent(input$help_button, {
        showModal(modalDialog(
          title = "Aide - Analyse Factorielle Exploratoire",
          tags$div(
            tags$h4("Comment utiliser cette application"),
            tags$ol(
              tags$li("Importez vos données à l'aide de la barre latérale"),
              tags$li("Examinez vos données dans l'onglet 'Aperçu des données'"),
              tags$li("Configurez et lancez votre analyse factorielle"),
              tags$li("Interprétez les résultats dans l'onglet 'Interprétation'")
            ),
            tags$h4("À propos de l'analyse factorielle"),
            tags$p("L'analyse factorielle exploratoire (AFE) permet d'identifier les facteurs latents à partir d'un ensemble de variables observées."),
            tags$p("Quelques conseils :"),
            tags$ul(
              tags$li("Vérifiez que les tests préliminaires (KMO, Bartlett) sont satisfaits"),
              tags$li("Utilisez le graphique des éboulis pour déterminer le nombre optimal de facteurs"),
              tags$li("Interprétez les facteurs en examinant les charges factorielles (> 0.3 ou 0.4)")
            )
          ),
          easyClose = TRUE,
          footer = modalButton("Fermer")
        ))
      })
      
      # Exécution de l'analyse factorielle avec lavaan
      observeEvent(input$run_efa, {
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        
        withProgress(message = "Calcul de l'analyse factorielle...", {
          tryCatch({
            # Sélection des variables pour l'analyse
            data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
            
            # Vérifier les NA et les gérer
            data_complete <- na.omit(data_subset)
            if (nrow(data_complete) < nrow(data_subset)) {
              showNotification(
                paste0(nrow(data_subset) - nrow(data_complete), " observations avec valeurs manquantes ont été supprimées"),
                type = "warning"
              )
            }
            
            # Augmenter la précision numérique
            op <- options(digits = 15)
            on.exit(options(op))
            
            # Vérifier la factorabilité de la matrice de corrélation
            data_complete <- na.omit(data_subset)
            cor_matrix <- cor(data_complete)
            
            # Vérifier les problèmes de multicollinéarité
            eigen_values <- eigen(cor_matrix)$values
            condition_number <- max(eigen_values) / min(eigen_values)
            
            if (condition_number > 1000) {
              showNotification(
                "Attention: Forte multicollinéarité détectée dans les données. Cela peut causer des problèmes d'estimation.",
                type = "warning",
                duration = 10
              )
            }
            
            # Analyse factorielle avec lavaan::efa()
            efa_model <- lavaan::efa(
              data = data_complete,
              nfactors = input$efa_nfactors,
              rotation = input$dynamic_rotation,
              estimator = input$efa_method,
              control = list(
                optim.method = "BFGS",
                optim.control = list(maxit = 1000),
                check.gradient = FALSE
              )
            )
            
            # Vérifier si des variances négatives ont été produites
            model_for_check <- if (inherits(efa_model, "efaList")) efa_model[[1]] else efa_model
            fit_measures_check <- fitMeasures(model_for_check)
            if (!is.null(attr(fit_measures_check, "warnings"))) {
              warning_messages <- attr(fit_measures_check, "warnings")
              
              # Si des avertissements de variance négative sont présents
              if (any(grepl("variance", warning_messages, ignore.case = TRUE))) {
                showNotification(
                  HTML("Attention: Des variances négatives ont été détectées.<br/>Envisagez de:<br/>- Réduire le nombre de facteurs<br/>- Exclure certaines variables<br/>- Augmenter la taille de l'échantillon"),
                  type = "warning",
                  duration = 15
                )
              }
            }
            
            # Utiliser summary() pour obtenir les résultats complets
            model_summary <- summary(efa_model, nd = 3, cutoff = 0.001)
            
            # Extraire les informations importantes du résumé
            if (inherits(efa_model, "efaList")) {
              # Pour les versions récentes de lavaan
              loadings_matrix <- model_summary$loading
              
              # Si la matrice de chargement n'est pas disponible directement
              if (is.null(loadings_matrix) || all(is.na(loadings_matrix))) {
                # Extraire le premier modèle de la liste
                model <- efa_model[[1]]
                params <- lavaan::parameterEstimates(model)
                vars <- unique(params$rhs[params$op == "=~"])
                factors <- unique(params$lhs[params$op == "=~"])
                
                loadings_matrix <- matrix(0, nrow = length(vars), ncol = length(factors))
                rownames(loadings_matrix) <- vars
                colnames(loadings_matrix) <- factors
                
                for (i in seq_len(nrow(params))) {
                  if (params$op[i] == "=~") {
                    f_idx <- match(params$lhs[i], factors)
                    v_idx <- match(params$rhs[i], vars)
                    loadings_matrix[v_idx, f_idx] <- params$est[i]
                  }
                }
              }
              
              # Renommer les colonnes pour la cohérence
              colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
              
              # Obtenir les communalités
              communalities <- model_summary$communalities
              if (is.null(communalities)) {
                communalities <- rowSums(loadings_matrix^2)
              }
              
              # Obtenir les mesures d'ajustement
              fit_indices <- lavaan::fitMeasures(efa_model[[1]])
              
            } else {
              # Pour les versions antérieures ou différentes de lavaan
              loadings_matrix <- lavaan::inspect(efa_model, "est")$lambda
              colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
              
              # Calculer les communalités
              communalities <- rowSums(loadings_matrix^2)
              
              # Obtenir les mesures d'ajustement
              fit_indices <- lavaan::fitMeasures(efa_model)
            }
            
            # Calculer les unicités
            uniqueness <- 1 - communalities
            
            # S'assurer que les noms sont corrects
            names(communalities) <- rownames(loadings_matrix)
            names(uniqueness) <- rownames(loadings_matrix)
            
            # Calculer la variance expliquée
            variance_explained <- colSums(loadings_matrix^2) / nrow(loadings_matrix)
            cumulative_variance <- cumsum(variance_explained)
            
            # Créer un objet de résultat complet
            efa_result(list(
              efa = efa_model,
              loadings = loadings_matrix,
              uniquenesses = uniqueness,
              communality = communalities,
              method = input$efa_method,
              rotation = input$dynamic_rotation,
              nfactors = input$efa_nfactors,
              Vaccounted = data.frame(
                SS_loadings = colSums(loadings_matrix^2),
                Proportion = variance_explained,
                Cumulative = cumulative_variance
              ),
              fit_indices = fit_indices,
              summary = model_summary  # Ajouter le résumé complet
            ))
            # Exécuter automatiquement la comparaison des modèles
            incProgress(0.6, message = "Comparaison des modèles...")
            
            # Définir la plage de facteurs à tester - limité par le nombre de variables
            n_variables <- length(input$efa_variables)
            # Borne de Ledermann : nombre maximal de facteurs identifiables pour n_variables indicateurs
            # (au-delà, lavaan::efa() lève une erreur car le modèle n'est plus identifié)
            max_factors <- floor((2 * n_variables + 1 - sqrt(8 * n_variables + 1)) / 2)
            max_factors <- min(max_factors, 5)  # plafond pratique pour la durée de calcul
            if (max_factors < 1) max_factors <- 1
            factors_to_test <- 1:max_factors
            
            # Exécuter l'EFA avec un vecteur de nfactors
            efa_models <- lavaan::efa(
              data = data_complete,
              nfactors = factors_to_test,
              rotation = input$dynamic_rotation,
              estimator = input$efa_method,
              control = list(
                optim.method = "BFGS",
                optim.control = list(maxit = 1000),
                check.gradient = FALSE
              )
            )
            
            # Obtenir les mesures d'ajustement pour tous les modèles
            fit_measures <- fitMeasures(efa_models, fit.measures = "all")
            
            # Convertir en dataframe avec des noms de colonnes appropriés
            if (is.matrix(fit_measures)) {
              fit_df <- as.data.frame(fit_measures)
              colnames(fit_df) <- paste0("nfct=", factors_to_test)
            } else {
              # Si seulement un modèle est testé
              fit_df <- data.frame(fit_measures)
              colnames(fit_df) <- paste0("nfct=", factors_to_test[1])
            }
            
            # Stocker les résultats de la comparaison
            model_comparison(list(
              fit_measures = fit_df,
              efa_models = efa_models
            ))
            # Sauvegarde de la solution sans rotation

            if (exists("data_complete") && !is.null(data_complete) && nrow(data_complete) > 0) {
              unrotated_solution(list(
                data_complete = data_complete,
                method = input$efa_method,
                nfactors = input$efa_nfactors
              ))
            } else {
              showNotification(
                "Erreur lors de la sauvegarde de la solution non rotée: données incomplètes",
                type = "error",
                duration = 5
              )
            }
            
            # Synchroniser le sélecteur de rotation dynamique 
            updateSelectInput(session, "dynamic_rotation", selected = input$dynamic_rotation)
            
            # Afficher l'onglet d'interprétation
            showTab(inputId = "main_tabs", target = "Interprétation")
            showTab(inputId = "main_tabs", target = "Bootstrap")
            showNotification("Analyse factorielle terminée", type = "message")
            
          }, error = function(e) {
            showNotification(
              paste("Erreur lors de l'analyse factorielle:", e$message),
              type = "error",
              duration = 10
            )
          })
        })
      })

      observeEvent(input$import_method, {
        if(input$import_method == "example") {
          # Mettre à jour les valeurs par défaut
          updateNumericInput(session, "sim_nobs", value = 200)
          updateNumericInput(session, "sim_nfactors", value = 2)
          updateNumericInput(session, "sim_nvars", value = 3)
          updateSliderInput(session, "sim_loading", value = 0.55)
          updateSliderInput(session, "sim_cross_loading", value = 0.25)
          updateSliderInput(session, "sim_error", value = 0.75)
          updateCheckboxInput(session, "sim_correlated", value = TRUE)
        }
      })
      
      # Gestion des données d'exemple
      observeEvent(input$generate_example, {
        # Nettoyer les objets existants avant de générer de nouvelles données
        cleanup_reactive_objects()
        withProgress(message = "Génération des données simulées...", {
          # Récupérer les paramètres 
          n <- input$sim_nobs
          num_factors <- input$sim_nfactors
          vars_per_factor <- input$sim_nvars
          loading_strength <- input$sim_loading
          cross_loading <- input$sim_cross_loading
          error_var <- input$sim_error
          correlated <- input$sim_correlated
          factor_cor <- input$sim_factor_cor
          
          # Initialiser un dataframe vide
          total_vars <- num_factors * vars_per_factor
          example_data <- matrix(0, nrow = n, ncol = total_vars)
          
          # Générer les facteurs latents
          current_time <- as.numeric(Sys.time())
          seed_value <- round(current_time %% 10000)
          set.seed(seed_value)
          
          # Générer les facteurs corrélés ou non
          if (correlated) {
            # Créer une matrice de corrélation pour les facteurs
            cor_matrix <- matrix(factor_cor, nrow = num_factors, ncol = num_factors)
            diag(cor_matrix) <- 1
            
            # Décomposition de Cholesky pour générer des facteurs corrélés
            L <- chol(cor_matrix)
            uncorrelated_factors <- matrix(rnorm(n * num_factors), nrow = n, ncol = num_factors)
            factors <- uncorrelated_factors %*% L
          } else {
            # Facteurs indépendants
            factors <- matrix(rnorm(n * num_factors), nrow = n, ncol = num_factors)
          }
          
          # Standardiser les facteurs pour assurer une variance de 1
          for (f in 1:num_factors) {
            factors[, f] <- scale(factors[, f])
          }
          
          # Générer les variables observées pour chaque facteur
          for (f in 1:num_factors) {
            # Indices des variables appartenant à ce facteur
            var_indices <- ((f-1) * vars_per_factor + 1):(f * vars_per_factor)
            
            # Pour chaque variable de ce facteur
            for (v in 1:vars_per_factor) {
              var_idx <- var_indices[v]
              
              # Initialiser avec zéro
              example_data[, var_idx] <- 0
              
              # Charge principale du facteur correspondant
              example_data[, var_idx] <- example_data[, var_idx] + loading_strength * factors[, f]
              
              # Ajouter des charges croisées (cross-loadings) des autres facteurs
              if (num_factors > 1) {
                other_factors <- setdiff(1:num_factors, f)
                for (of in other_factors) {
                  example_data[, var_idx] <- example_data[, var_idx] + cross_loading * factors[, of]
                }
              }
              
              # Ajouter un bruit aléatoire - important: le bruit doit être proportionnel à la variance non expliquée
              # Calculer la variance expliquée par les facteurs
              explained_var <- loading_strength^2 + (num_factors - 1) * cross_loading^2
              # Le bruit doit représenter (1 - explained_var) pour que la variance totale soit ~1
              # Mais nous voulons l'ajuster par le facteur error_var de l'utilisateur
              remaining_var <- 1 - explained_var
              noise_sd <- sqrt(remaining_var) * error_var
              
              # Ajouter le bruit avec la bonne échelle
              example_data[, var_idx] <- example_data[, var_idx] + rnorm(n, 0, noise_sd)
              
              # Standardiser chaque variable pour avoir une variance de 1
              example_data[, var_idx] <- scale(example_data[, var_idx])
            }
          }
          
          # Convertir en dataframe avec des noms de variables simples X1, X2, ...
          example_data <- as.data.frame(example_data)
          colnames(example_data) <- paste0("X", 1:total_vars)
          
          # Stocker les données
          raw_data(example_data)
          filtered_data(example_data)
          
          # Mettre à jour les sélecteurs
          updateSelectInput(session, "filter_var",
                            choices = c("Aucun" = "", names(example_data)))
          
          updateSelectInput(session, "efa_variables",
                            choices = names(example_data),
                            selected = names(example_data))
          
          # Mettre à jour le nombre de facteurs recommandé
          updateNumericInput(session, "efa_nfactors", 
                             value = num_factors)
          
          # Afficher l'onglet d'analyse
          showTab(inputId = "main_tabs", target = "Analyse factorielle")
          
          # Replier "Importation des données" et ouvrir "Paramètres de l'analyse"
          accordion_panel_set(id = "sidebar_accordion", values = "Paramètres de l'analyse")
          
          # Ajouter une notification avec un résumé des données simulées
          showNotification(
            HTML(paste0(
              "<b>Données simulées générées</b><br/>",
              "Observations: ", n, "<br/>",
              "Structure: ", num_factors, " facteurs avec ", vars_per_factor, " variables chacun<br/>",
              "Charges principales: ", loading_strength, "<br/>",
              "Facteurs ", ifelse(correlated, "corrélés", "non corrélés")
            )),
            type = "message",
            duration = 5
          )
        })
      })
      
      
      # Résumé du fichier - Chaque élément sur une ligne séparée
      output$file_summary_compact <- renderUI({
        req(input$data_file)
        
        # Déterminer le type de fichier
        file_ext <- tools::file_ext(input$data_file$name)
        file_size <- round(input$data_file$size / 1024, 1)
        
        # Style CSS pour les labels et les valeurs
        label_style <- "font-weight: bold; margin-right: 5px;"
        
        # Création d'un résumé avec chaque élément sur une ligne séparée
        div(
          div(
            span(style = label_style, "Fichier:"), 
            span(input$data_file$name)
          ),
          div(
            span(style = label_style, "Type:"), 
            span(toupper(file_ext))
          ),
          div(
            span(style = label_style, "Taille:"), 
            span(paste0(file_size, " ko"))
          )
        )
      })
      
      # Affichage de l'aperçu des données - limité à 4 décimales
      output$data_preview <- renderDT({
        if(is.null(filtered_data())) {
          return(datatable(data.frame(Message = "Aucune donnée chargée. Veuillez importer un fichier."),
                           options = list(dom = 't'),
                           rownames = FALSE))
        }
        
        # Créer une copie des données pour formater les nombres
        formatted_data <- as.data.frame(filtered_data())
        
        # Formater les colonnes numériques pour afficher au maximum 4 décimales
        for(col in names(formatted_data)) {
          if(is.numeric(formatted_data[[col]])) {
            formatted_data[[col]] <- round(formatted_data[[col]], 4)
          }
        }
        
        datatable(
          formatted_data,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 'ftip',
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/French.json')
          ),
          rownames = FALSE
        )
      })
      
      # Contrôles de filtre dynamiques
      output$filter_controls <- renderUI({
        req(input$filter_var, raw_data())
        
        if (input$filter_var == "") {
          return(NULL)
        }
        
        var <- raw_data()[[input$filter_var]]
        
        if (is.numeric(var) || is.integer(var)) {
          min_val <- min(var, na.rm = TRUE)
          max_val <- max(var, na.rm = TRUE)
          
          sliderInput("filter_range", 
                      label = paste("Plage de valeurs pour", input$filter_var),
                      min = min_val, 
                      max = max_val,
                      value = c(min_val, max_val),
                      step = (max_val - min_val) / 100)
        } else if (is.factor(var) || is.character(var)) {
          unique_vals <- unique(var)
          checkboxGroupInput("filter_categories",
                             label = paste("Catégories pour", input$filter_var),
                             choices = unique_vals,
                             selected = unique_vals)
        } else {
          return(tags$div("Type de variable non supporté pour le filtrage"))
        }
      })
      
      # Output pour afficher la comparaison des modèles
      output$model_comparison_output <- renderPrint({
        req(model_comparison())
        
        cat("COMPARAISON DES MODÈLES SELON LE NOMBRE DE FACTEURS\n")
        cat("=================================================\n\n")
        
        # Extraire et formater les mesures d'ajustement
        fit_df <- model_comparison()$fit_measures
        
        # Sélectionner seulement les indices importants
        important_indices <- c("npar", "fmin", "chisq", "df", "pvalue", 
                               "cfi", "tli", "nnfi", "rfi", "nfi", 
                               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", 
                               "srmr", "aic", "bic")
        
        if (all(important_indices %in% rownames(fit_df))) {
          fit_df <- fit_df[important_indices, ]
        }
        
        print(round(fit_df, 3))
      })
      outputOptions(output, "model_comparison_output", suspendWhenHidden = FALSE)
      
      output$download_simulated_data <- downloadHandler(
        filename = function() {
          ext <- input$download_format
          paste("donnees_simulees_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext, sep = "")
        },
        content = function(file) {
          req(raw_data())
          data <- raw_data()
          
          # Exporter selon le format sélectionné
          if (input$download_format == "csv") {
            write.csv(data, file, row.names = FALSE, quote = TRUE)
          } else if (input$download_format == "xlsx") {
            writexl::write_xlsx(data, path = file)
          } else if (input$download_format == "sav") {
            haven::write_sav(data, file)
          } else if (input$download_format == "dta") {
            haven::write_dta(data, file)
          }
        }
      )
      
      output$correlation_output <- renderUI({
        req(filtered_data())
        
        if(input$corr_display_type == "plot") {
          plotOutput("correlation_plot", height = "700px", width = "100%")  # Augmenter la hauteur
        } else {
          verbatimTextOutput("correlation_text")
        }
      })
      outputOptions(output, "correlation_output", suspendWhenHidden = FALSE)
      
      # Modifiez cette section dans le code serveur pour la matrice de corrélation graphique
      output$correlation_plot <- renderPlot({
        req(filtered_data())
        
        # Sélectionner uniquement les variables numériques
        data <- filtered_data()
        numeric_cols <- sapply(data, is.numeric)
        
        if(sum(numeric_cols) < 2) {
          # Si moins de 2 variables numériques, afficher un message
          plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(0, 0, "Pas assez de variables numériques\npour afficher une matrice de corrélation", cex = 1.2)
        } else {
          # Limiter à 30 variables numériques pour la lisibilité
          numeric_data <- data[, numeric_cols, drop = FALSE]
          if(ncol(numeric_data) > 30) {
            numeric_data <- numeric_data[, 1:30]
            title_cor <- "Matrice de corrélation (30 premières variables numériques)"
          } else {
            title_cor <- "Matrice de corrélation"
          }
          
          # Calculer la matrice de corrélation
          cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
          
          # Afficher avec corrplot.mixed
          corrplot.mixed(cor_matrix, 
                         lower = "number", 
                         upper = "color",
                         tl.col = "black",
                         tl.pos = "lt",
                         diag = "l",
                         number.cex = 0.9,
                         title = title_cor,
                         mar = c(0, 0, 2, 0))
        }
      }, height = 800, width = 900)  # Augmenter la hauteur et la largeur du graphique
      outputOptions(output, "correlation_plot", suspendWhenHidden = FALSE)
      
      # Ajouter cette nouvelle fonction dans la partie serveur
      output$correlation_text <- renderPrint({
        req(filtered_data())
        
        # Sélectionner uniquement les variables numériques
        data <- filtered_data()
        numeric_cols <- sapply(data, is.numeric)
        
        if(sum(numeric_cols) < 2) {
          cat("Pas assez de variables numériques pour afficher une matrice de corrélation")
        } else {
          # Utiliser toutes les variables numériques sans limitation pour l'affichage texte
          numeric_data <- data[, numeric_cols, drop = FALSE]
          
          # Calculer la matrice de corrélation
          cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
          
          # Afficher un message informatif
          if(ncol(numeric_data) > 20) {
            cat("Matrice de corrélation complète (", ncol(numeric_data), " variables)\n\n")
          } else {
            cat("Matrice de corrélation\n\n")
          }
          
          # Afficher la matrice arrondie à 2 décimales
          print(round(cor_matrix, 2))
        }
      })
      outputOptions(output, "correlation_text", suspendWhenHidden = FALSE)
      
      # Mettre à jour automatiquement les tests de normalité quand les données changent
      observe({
        req(filtered_data())
        
        # Appliquer les tests de normalité automatiquement à toutes les variables numériques
        data <- filtered_data()
        numeric_vars <- names(data)[sapply(data, is.numeric)]
        
        # Limiter à 20 variables pour des raisons de performance
        if (length(numeric_vars) > 20) {
          numeric_vars <- numeric_vars[1:20]
        }
        
        # Créer un tableau pour stocker les résultats
        results <- data.frame(
          Variable = character(),
          Shapiro_Stat = numeric(),
          Shapiro_p = numeric(),
          Normal = character(),
          stringsAsFactors = FALSE
        )
        
        # Exécuter les tests pour chaque variable
        for (var_name in numeric_vars) {
          var_data <- data[[var_name]]
          
          # Vérifier que la variable a suffisamment de données
          if (sum(!is.na(var_data)) >= 3) {
            # Limiter à 5000 observations pour le test de Shapiro (limitation de la méthode)
            if (length(var_data) > 5000) {
              var_data <- sample(var_data[!is.na(var_data)], 5000)
            } else {
              var_data <- var_data[!is.na(var_data)]
            }
            
            # Test de Shapiro-Wilk
            tryCatch({
              shapiro_test <- shapiro.test(var_data)
              
              # Ajouter les résultats
              results <- rbind(results, data.frame(
                Variable = var_name,
                Shapiro_Stat = round(shapiro_test$statistic, 4),
                Shapiro_p = round(shapiro_test$p.value, 4),
                Normal = ifelse(shapiro_test$p.value > 0.05, "Oui", "Non"),
                stringsAsFactors = FALSE
              ))
            }, error = function(e) {
              # En cas d'erreur, ajouter une ligne avec message d'erreur
              results <- rbind(results, data.frame(
                Variable = var_name,
                Shapiro_Stat = NA,
                Shapiro_p = NA,
                Normal = "Erreur de test",
                stringsAsFactors = FALSE
              ))
            })
          }
        }
        
        # Stocker les résultats dans la réactive
        normality_results_data(results)
      })
      
      # Afficher les résultats des tests de normalité dans un tableau interactif
      output$normality_results_table <- renderDT({
        req(normality_results_data())
        results <- normality_results_data()
        
        if (nrow(results) == 0) {
          return(datatable(data.frame(Message = "Aucune variable numérique trouvée pour les tests de normalité"),
                           options = list(dom = 't'),
                           rownames = FALSE))
        }
        
        # Préparer un tableau formaté pour l'affichage
        results_table <- data.frame(
          Variable = results$Variable,
          "Stat. Shapiro" = results$Shapiro_Stat,
          "p-value" = results$Shapiro_p,
          "Distribution normale" = ifelse(results$Normal == "Oui", 
                                          "Oui (p > 0.05)", 
                                          ifelse(results$Normal == "Non", 
                                                 "Non (p < 0.05)", 
                                                 results$Normal)),
          check.names = FALSE
        )
        
        # Créer le tableau interactif avec coloration des p-values
        datatable(
          results_table,
          options = list(
            pageLength = 15,
            dom = 'ftip',
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/French.json')
          ),
          rownames = FALSE
        ) %>%
          formatStyle(
            "p-value",
            color = styleInterval(0.05, c("red", "green")),
            fontWeight = styleInterval(0.05, c("bold", "normal"))
          )
      })
      outputOptions(output, "normality_results_table", suspendWhenHidden = FALSE)
      
      # Mettre à jour le selectInput pour la matrice de visualisations
      observe({
        req(filtered_data())
        data <- filtered_data()
        
        # Mettre à jour les sélecteurs avec les variables numériques
        numeric_vars <- names(data)[sapply(data, is.numeric)]
        updateSelectInput(session, "viz_variables", choices = numeric_vars)
      })
      
      # Créer la matrice de visualisations
      output$matrix_plot <- renderPlot({
        req(filtered_data(), input$viz_variables)
        
        # Vérifier qu'il y a des variables sélectionnées
        if (length(input$viz_variables) == 0) {
          plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(0, 0, "Sélectionnez au moins une variable à visualiser", cex = 1.5)
          return()
        }
        
        # Limiter à 8 variables maximum
        vars_to_plot <- input$viz_variables
        if (length(vars_to_plot) > 8) {
          vars_to_plot <- vars_to_plot[1:8]
        }
        
        # Récupérer les données pour ces variables
        data <- filtered_data()
        plot_data <- data[, vars_to_plot, drop = FALSE]
        
        # Déterminer la disposition de la matrice (nombre de lignes et colonnes)
        n_vars <- length(vars_to_plot)
        n_cols <- min(4, n_vars)  # Maximum 4 colonnes
        n_rows <- ceiling(n_vars / n_cols)
        
        # Configurer une mise en page en matrice
        par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
        
        # Créer un graphique pour chaque variable
        for (var in vars_to_plot) {
          var_data <- data[[var]]
          var_data <- var_data[!is.na(var_data)]
          
          if (length(var_data) == 0) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, "Pas de données", cex = 1.2)
            next
          }
          
          if (input$matrix_plot_type == "hist") {
            # Histogramme
            h <- hist(var_data, plot = FALSE, breaks = min(30, max(10, round(length(var_data)/30))))
            
            hist(var_data, 
                 breaks = h$breaks,
                 main = var,
                 xlab = "",
                 ylab = "Fréquence",
                 col = pal_efa$hist_fill,
                 border = "white",
                 prob = input$matrix_density)
            
            if (input$matrix_density) {
              d <- density(var_data, na.rm = TRUE)
              lines(d, col = pal_efa$density_obs, lwd = 2)
              x <- seq(min(var_data), max(var_data), length.out = 100)
              lines(x, dnorm(x, mean = mean(var_data), sd = sd(var_data)), 
                    col = pal_efa$density_theo, lty = 2, lwd = 2)
            }
            
          } else if (input$matrix_plot_type == "box") {
            # Boxplot
            boxplot(var_data, 
                    main = var,
                    xlab = "",
                    ylab = "",
                    col = pal_efa$box_fill,
                    horizontal = TRUE)
            
          } else if (input$matrix_plot_type == "qq") {
            # Q-Q plot
            qqnorm(var_data, 
                   main = var,
                   col = pal_efa$qq_points,
                   pch = 19)
            qqline(var_data, col = pal_efa$density_obs, lwd = 2)
          }
        }
        
        # Ajouter une légende commune si c'est un histogramme avec densité
        if (input$matrix_plot_type == "hist" && input$matrix_density) {
          # Ajouter une légende dans la marge externe supérieure
          par(fig = c(0, 1, 0, 1), oma = c(0, 0, 3, 0), mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          legend("top", 
                 legend = c("Densité observée", "Distribution normale"), 
                 col = c(pal_efa$density_obs, pal_efa$density_theo), 
                 lty = c(1, 2), 
                 lwd = 2,
                 horiz = TRUE,
                 bty = "n")
        }
      })
      outputOptions(output, "matrix_plot", suspendWhenHidden = FALSE)
      
      # Réactive pour stocker les résultats de la CFA
      cfa_result <- reactiveVal(NULL)
      
      # Fonction pour générer le modèle CFA à partir de l'EFA
      observeEvent(input$generate_cfa, {
        req(efa_result())
        
        # Générer la syntaxe CFA
        cfa_syntax <- convert_efa_to_cfa(efa_result(), input$cfa_threshold)
        
        # Mettre à jour le champ de texte avec la syntaxe
        updateTextAreaInput(session, "cfa_syntax", value = cfa_syntax)
      })
      
      # Fonction pour exécuter l'analyse confirmatoire
      observeEvent(input$run_cfa, {
        req(input$cfa_syntax, filtered_data())
        
        withProgress(message = "Exécution de l'analyse confirmatoire...", {
          tryCatch({
            # Obtenir la syntaxe et les données
            cfa_model_syntax <- input$cfa_syntax
            data <- filtered_data()
            
            # Sélectionner seulement les variables utilisées dans le modèle
            used_vars <- lavaan::lavParseModelString(cfa_model_syntax)
            used_vars <- unique(c(used_vars$lhs, used_vars$rhs))
            used_vars <- intersect(used_vars, names(data))
            
            # Vérifier la factoriabilité avant l'analyse
            if (length(used_vars) > 1) {
              data_subset <- data[, used_vars, drop = FALSE]
              data_subset <- na.omit(data_subset)
              
              # Calculer le déterminant de la matrice de corrélation
              cor_matrix <- cor(data_subset)
              det_cor <- det(cor_matrix)
              
              # Avertir si le déterminant est trop proche de zéro
              if (det_cor < 0.00001) {
                showNotification(
                  HTML(paste0(
                    "<b>Attention: Forte dépendance linéaire détectée</b><br/>",
                    "Le déterminant de la matrice de corrélation est très proche de zéro (", round(det_cor, 8), ").<br/>",
                    "Cela peut indiquer que certaines variables sont presque parfaitement corrélées, ",
                    "ce qui peut causer des problèmes d'estimation."
                  )),
                  type = "warning",
                  duration = 10
                )
              }
            }
            
            # Exécuter le modèle CFA avec la méthode sélectionnée
            cfa_model <- lavaan::cfa(
              model = cfa_model_syntax,
              data = data,
 #             missing = "fiml",
              std.lv = TRUE,
              estimator = input$cfa_method,  # Utiliser la méthode sélectionnée
              se = ifelse(input$cfa_method %in% c("MLR", "MLM", "MLMV"), "robust", "standard"),
              control = list(
                iter.max = 10000,
                optimizer = "NLMINB",
                ridge = 0.1
              )
            )
            
            # Vérifier si des problèmes existent après ajustement
            fit_measures <- fitMeasures(cfa_model)
            
            # Vérifier les variances négatives
            param_est <- parameterEstimates(cfa_model)
            var_params <- param_est[param_est$op == "~~" & param_est$lhs == param_est$rhs, ]
            
            # S'il y a des variances négatives
            if (any(var_params$est < 0)) {
              negative_vars <- var_params$lhs[var_params$est < 0]
              
              showNotification(
                HTML(paste0(
                  "<b>Variances négatives détectées</b><br/>",
                  "Variables: ", paste(negative_vars, collapse=", "), "<br/>",
                  "Cela peut indiquer un sur-ajustement ou un problème de spécification.<br/>",
                  "Suggestions:<br/>",
                  "- Simplifier votre modèle<br/>",
                  "- Essayer une méthode d'estimation différente (ex: ULS)<br/>",
                  "- Ajouter des contraintes de non-négativité"
                )),
                type = "warning",
                duration = 10
              )
            }
            
            # Stocker les résultats
            cfa_result(cfa_model)
            
            # Notification de succès
            showNotification(paste0("Analyse confirmatoire terminée (Méthode: ", input$cfa_method, ")"), type = "message")
            
          }, error = function(e) {
            # En cas d'erreur, essayer avec ULS
            tryCatch({
              showNotification(
                paste0("Échec avec la méthode ", input$cfa_method, ". Tentative avec ULS..."),
                type = "warning",
                duration = 5
              )
              
              # Obtenir la syntaxe et les données
              cfa_model_syntax <- input$cfa_syntax
              data <- filtered_data()
              
              # Essayer avec ULS qui est moins sensible aux problèmes numériques
              cfa_model <- lavaan::cfa(
                model = cfa_model_syntax,
                data = data,
 #               missing = "fiml",
                std.lv = TRUE,
                estimator = "ULS",
                control = list(
                  iter.max = 15000,
                  optimizer = "NLMINB",
                  ridge = 0.5
                )
              )
              
              # Stocker les résultats
              cfa_result(cfa_model)
              
              showNotification(
                "Analyse confirmatoire terminée avec méthode alternative (ULS)",
                type = "warning",
                duration = 7
              )
              
            }, error = function(e2) {
              # Si même l'ULS échoue
              showNotification(
                HTML(paste0(
                  "<b>Erreur lors de l'analyse confirmatoire</b><br/>",
                  e$message, "<br/><br/>",
                  "Votre modèle pourrait avoir des problèmes d'identification ou de spécification."
                )),
                type = "error",
                duration = 15
              )
            })
          })
        })
      })
      
      # Afficher le résumé du modèle CFA
      output$cfa_summary <- renderPrint({
        req(cfa_result())
        summary(cfa_result(), fit.measures = TRUE, standardized = TRUE)
      })
      outputOptions(output, "cfa_summary", suspendWhenHidden = FALSE)
      
      # Afficher les indices d'ajustement
      output$cfa_fit_indices <- renderPrint({
        req(cfa_result())
        cat("INDICES D'AJUSTEMENT DU MODÈLE CFA\n")
        cat("================================\n\n")
        fit_indices <- fitMeasures(cfa_result(), fit.measures = "all")
        print(fit_indices)
      })
      outputOptions(output, "cfa_fit_indices", suspendWhenHidden = FALSE)
      
      # Diagramme du modèle CFA
      output$cfa_diagram <- renderPlot({
        req(cfa_result())
        
        withProgress(message = "Génération du diagramme du modèle...", {
          # Générer un diagramme du modèle avec semPaths
          tryCatch({
            if (!requireNamespace("semPlot", quietly = TRUE)) {
              plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
              text(0, 0, "Package 'semPlot' requis pour le diagramme", cex = 1.2)
            } else {
              semPlot::semPaths(cfa_result(),
                                intercepts = FALSE,
                                whatLabels = "std",
                                rotation = 2,
                                sizeMan = 6,
                                sizeMan2 = 3,
                                sizeLat = 8,
                                style = "lisrel",
                                curvature = 1.8,
                                edge.label.cex = .8,
                                label.cex = .8,
                                residScale = 12)
            }
          }, error = function(e) {
            plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(0, 0, paste("Erreur lors de la génération du diagramme:", e$message), cex = 1.2)
          })
        })
      })
      outputOptions(output, "cfa_diagram", suspendWhenHidden = FALSE)
      
      # Afficher les indices de modification
      output$modification_indices <- renderDT({
        req(cfa_result())
        
        # Obtenir les indices de modification
        mod_indices <- modificationIndices(cfa_result(), sort = TRUE)
        
        # Si aucun index de modification, afficher un message
        if (nrow(mod_indices) == 0) {
          return(datatable(data.frame(Message = "Aucun indice de modification disponible"),
                           options = list(dom = 't'),
                           rownames = FALSE))
        }
        
        # Filtrer et trier les indices
        mod_indices <- mod_indices[mod_indices$mi > 3.84, ]  # Seuil chi² à df=1, p<0.05
        mod_indices <- mod_indices[order(-mod_indices$mi), ]
        
        # Formater pour l'affichage
        mod_indices_formatted <- mod_indices[, c("lhs", "op", "rhs", "mi", "epc")]
        colnames(mod_indices_formatted) <- c("De", "Opération", "Vers", "IM", "CPC")
        
        # Créer le tableau
        datatable(
          mod_indices_formatted,
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("IM", "CPC"), digits = 2) %>%
          formatStyle(
            "IM",
            background = styleInterval(c(3.84, 10), c(pal_efa$mi_low, pal_efa$mi_medium, pal_efa$mi_high)),
            fontWeight = styleInterval(10, c("normal", "bold"))
          )
      })
      outputOptions(output, "modification_indices", suspendWhenHidden = FALSE)
      
      # Mise à jour de l'état après l'exécution de l'EFA
      observeEvent(efa_result(), {
        # Si l'EFA est terminée, afficher l'onglet CFA
        if (!is.null(efa_result())) {
          showTab(inputId = "main_tabs", target = "Analyse Confirmatoire")
        }
      })
      
      # Sortie pour le déterminant de la matrice de corrélation
      output$determinant_test <- renderPrint({
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        
        data_subset <- filtered_data()[, input$efa_variables, drop = FALSE]
        data_complete <- na.omit(data_subset)
        
        # Calculer la matrice de corrélation et son déterminant
        cor_matrix <- cor(data_complete)
        det_cor <- det(cor_matrix)
        
        cat("Déterminant de la matrice de corrélation\n")
        cat("-------------------------------------\n")
        cat("Déterminant:", format(det_cor, digits = 5, scientific = TRUE), "\n\n")
        
        # Interprétation du déterminant
        if (det_cor > 0.00001) {
          cat("Interprétation: La matrice est bien conditionnée.\n")
          cat("Les variables ne sont pas trop fortement corrélées entre elles.\n")
        } else if (det_cor > 0) {
          cat("Interprétation: Déterminant très proche de zéro.\n")
          cat("Attention: Forte multicollinéarité potentielle.\n")
          cat("Certaines variables peuvent être trop fortement corrélées,\n")
          cat("ce qui peut causer des problèmes d'estimation.\n")
        } else if (det_cor == 0) {
          cat("Interprétation: Déterminant égal à zéro.\n")
          cat("Alerte: La matrice est singulière (non inversible).\n")
          cat("Il existe une dépendance linéaire parfaite entre certaines variables.\n")
          cat("L'analyse factorielle peut être impossible sans modification des données.\n")
        } else {
          cat("Erreur de calcul ou matrice non définie positive.\n")
        }
      })
      
      # Calculer les indices de fiabilité pour chaque facteur
      calculate_reliability <- function(data, loadings, threshold = 0.3) {
        req(data, loadings)
        
        # Nombre de facteurs
        n_factors <- ncol(loadings)
        
        # Initialiser le tableau de résultats
        reliability_results <- data.frame(
          Facteur = paste0("F", 1:n_factors),
          Alpha = numeric(n_factors),
          Omega = numeric(n_factors),
          Variables = character(n_factors)
        )
        
        # Pour chaque facteur
        for (f in 1:n_factors) {
          # Identifier les variables qui chargent sur ce facteur au-dessus du seuil
          items <- rownames(loadings)[abs(loadings[, f]) >= threshold]
          
          # Si au moins 2 variables sont associées au facteur
          if (length(items) >= 2) {
            # Extraire les données pour ces variables
            factor_data <- data[, items, drop = FALSE]
            
            # Calculer alpha de Cronbach
            alpha_result <- tryCatch({
              psych::alpha(factor_data)$total$raw_alpha
            }, error = function(e) {
              NA
            })
            
            # Calculer oméga de McDonald
            omega_result <- tryCatch({
              psych::omega(factor_data, nfactors = 1, plot = FALSE)$omega.tot
            }, error = function(e) {
              NA
            })
            
            # Stocker les résultats
            reliability_results$Alpha[f] <- alpha_result
            reliability_results$Omega[f] <- omega_result
            reliability_results$Variables[f] <- paste(items, collapse = ", ")
          } else {
            # Si moins de 2 variables, la fiabilité ne peut pas être calculée
            reliability_results$Alpha[f] <- NA
            reliability_results$Omega[f] <- NA
            reliability_results$Variables[f] <- ifelse(length(items) == 0, 
                                                       "Aucune variable (seuil trop élevé)",
                                                       paste(items, collapse = ", "))
          }
        }
        
        return(reliability_results)
      }
      
      # Mettre à jour le sélecteur de facteurs lorsque l'EFA est terminée
      observe({
        req(efa_result())
        efa <- efa_result()
        
        # Mettre à jour le sélecteur de facteurs avec les facteurs disponibles
        updateSelectInput(session, "factor_for_analysis",
                          choices = setNames(1:efa$nfactors, 
                                             paste0("Facteur ", 1:efa$nfactors)))
      })
      
      # Tableau des indices de fiabilité
      output$reliability_table <- renderDT({
        req(efa_result(), filtered_data())
        
        efa <- efa_result()
        data <- filtered_data()
        loadings <- efa$loadings
        threshold <- input$loadings_threshold
        
        # Calculer les indices de fiabilité
        reliability_data <- calculate_reliability(data, loadings, threshold)
        
        # Formatage pour l'affichage selon la méthode sélectionnée
        if (input$reliability_method == "alpha") {
          display_data <- reliability_data[, c("Facteur", "Alpha", "Variables")]
          colnames(display_data) <- c("Facteur", "Alpha de Cronbach", "Variables incluses")
          
          # Créer le tableau avec coloration selon la valeur - ALPHA SEULEMENT
          datatable(
            display_data,
            options = list(pageLength = 25, scrollX = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(
              columns = "Alpha de Cronbach",
              digits = 4
            ) %>%
            formatStyle(
              columns = "Alpha de Cronbach",
              backgroundColor = styleInterval(
                c(0.6, 0.7, 0.8),
                c(pal_efa$scale4_poor, pal_efa$scale4_medium, pal_efa$scale4_good, pal_efa$scale4_excellent)
              ),
              color = "black"
            )
          
        } else if (input$reliability_method == "omega") {
          display_data <- reliability_data[, c("Facteur", "Omega", "Variables")]
          colnames(display_data) <- c("Facteur", "Oméga de McDonald", "Variables incluses")
          
          # Créer le tableau avec coloration selon la valeur - OMEGA SEULEMENT
          datatable(
            display_data,
            options = list(pageLength = 25, scrollX = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(
              columns = "Oméga de McDonald",
              digits = 4
            ) %>%
            formatStyle(
              columns = "Oméga de McDonald",
              backgroundColor = styleInterval(
                c(0.6, 0.7, 0.8),
                c(pal_efa$scale4_poor, pal_efa$scale4_medium, pal_efa$scale4_good, pal_efa$scale4_excellent)
              ),
              color = "black"
            )
          
        } else { # both
          display_data <- reliability_data[, c("Facteur", "Alpha", "Omega", "Variables")]
          colnames(display_data) <- c("Facteur", "Alpha de Cronbach", "Oméga de McDonald", "Variables incluses")
          
          # Créer le tableau avec coloration selon la valeur - LES DEUX
          datatable(
            display_data,
            options = list(pageLength = 25, scrollX = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(
              columns = c("Alpha de Cronbach", "Oméga de McDonald"),
              digits = 4
            ) %>%
            formatStyle(
              columns = c("Alpha de Cronbach", "Oméga de McDonald"),
              backgroundColor = styleInterval(
                c(0.6, 0.7, 0.8),
                c(pal_efa$scale4_poor, pal_efa$scale4_medium, pal_efa$scale4_good, pal_efa$scale4_excellent)
              ),
              color = "black"
            )
        }
      })
      outputOptions(output, "reliability_table", suspendWhenHidden = FALSE)
      
      # Analyse détaillée de fiabilité pour un facteur spécifique
      output$detailed_reliability <- renderPrint({
        req(efa_result(), filtered_data(), input$factor_for_analysis)
        
        efa <- efa_result()
        data <- filtered_data()
        loadings <- efa$loadings
        threshold <- input$loadings_threshold
        factor_idx <- as.numeric(input$factor_for_analysis)
        
        # Identifier les variables qui chargent sur ce facteur
        items <- rownames(loadings)[abs(loadings[, factor_idx]) >= threshold]
        
        if (length(items) < 2) {
          cat("Fiabilité non calculable : le facteur doit avoir au moins 2 variables avec des charges\n")
          cat("supérieures au seuil (", threshold, ").\n\n")
          cat("Variables associées à ce facteur (si existantes) :\n")
          cat(paste(items, collapse = ", "))
          return()
        }
        
        # Extraire les données pour ces variables
        factor_data <- data[, items, drop = FALSE]
        
        # Afficher les charges factorielles pour ces variables
        cat("ANALYSE DE FIABILITÉ - FACTEUR ", factor_idx, "\n")
        cat("======================================\n\n")
        cat("Variables incluses (seuil >= ", threshold, "):\n")
        item_loadings <- data.frame(
          Variable = items,
          Charge = loadings[items, factor_idx]
        )
        print(item_loadings)
        cat("\n")
        
        # Calculer et afficher alpha de Cronbach avec détails
        cat("ALPHA DE CRONBACH\n")
        cat("----------------\n")
        tryCatch({
          alpha_result <- psych::alpha(factor_data)
          print(alpha_result$total)
          cat("\nStatistiques par item:\n")
          print(alpha_result$item.stats)
          cat("\nAlpha si l'item est supprimé:\n")
          alpha_drop <- alpha_result$alpha.drop
          if (!is.null(alpha_drop)) {
            print(alpha_drop[, c("raw_alpha", "std.alpha")])
          } else {
            cat("Non disponible\n")
          }
        }, error = function(e) {
          cat("Erreur lors du calcul de l'alpha :", e$message, "\n")
        })
        
        cat("\n\nOMÉGA DE MCDONALD\n")
        cat("----------------\n")
        tryCatch({
          omega_result <- psych::omega(factor_data, nfactors = 1, plot = FALSE)
          cat("Oméga total :", round(omega_result$omega.tot, 3), "\n")
          cat("Oméga hiérarchique :", round(omega_result$omega_h, 3), "\n")
          cat("Variance totale expliquée :", round(omega_result$total.var, 3), "\n")
        }, error = function(e) {
          cat("Erreur lors du calcul de l'oméga :", e$message, "\n")
        })
      })
      outputOptions(output, "detailed_reliability", suspendWhenHidden = FALSE)
      
      # Réactive pour stocker les résultats du bootstrap
      bootstrap_results <- reactiveVal(NULL)
      
      # Fonction pour exécuter l'analyse bootstrap
      observeEvent(input$run_bootstrap, {
        req(filtered_data(), efa_result(), input$efa_variables)
        
        # Récupérer les paramètres nécessaires
        data <- filtered_data()
        efa <- efa_result()
        selected_vars <- input$efa_variables
        nfactors <- efa$nfactors
        rotation <- efa$rotation
        method <- efa$method
        n_samples <- input$bootstrap_samples
        ci_level <- input$bootstrap_ci / 100
        
        # Vérifier que les données sont suffisantes
        if (length(selected_vars) < nfactors + 1) {
          showNotification(
            "Erreur: Trop peu de variables pour le nombre de facteurs spécifié",
            type = "error",
            duration = 5
          )
          return()
        }
        
        withProgress(message = "Exécution du bootstrap...", {
          tryCatch({
            # Créer une liste pour stocker les résultats
            all_loadings <- list()
            
            # Sélectionner les variables pour l'analyse
            data_subset <- data[, selected_vars, drop = FALSE]
            data_complete <- na.omit(data_subset)
            n_obs <- nrow(data_complete)
            
            incProgress(0.1, detail = "Initialisation...")
            
            # Boucle de bootstrap
            for (i in 1:n_samples) {
              # Mettre à jour la barre de progression tous les 10 échantillons
              if (i %% 10 == 0) {
                incProgress(0.8 * (10/n_samples), 
                            detail = paste0("Échantillon ", i, " sur ", n_samples))
              }
              
              # Créer un échantillon bootstrap (tirage avec remise)
              boot_indices <- sample(1:n_obs, n_obs, replace = TRUE)
              boot_data <- data_complete[boot_indices, ]
              
              # Essayer d'exécuter l'EFA sur l'échantillon bootstrap
              tryCatch({
                boot_efa <- lavaan::efa(
                  data = boot_data,
                  nfactors = nfactors,
                  rotation = rotation,
                  estimator = method,
                  control = list(
                    optim.method = "BFGS",
                    optim.control = list(maxit = 500),
                    check.gradient = FALSE
                  )
                )
                
                # Extraire les charges factorielles
                if (inherits(boot_efa, "efaList")) {
                  model_summary <- summary(boot_efa, nd = 3, cutoff = 0.001)
                  loadings_matrix <- model_summary$loading
                  
                  if (is.null(loadings_matrix) || all(is.na(loadings_matrix))) {
                    model <- boot_efa[[1]]
                    params <- lavaan::parameterEstimates(model)
                    vars <- unique(params$rhs[params$op == "=~"])
                    factors <- unique(params$lhs[params$op == "=~"])
                    
                    loadings_matrix <- matrix(0, nrow = length(vars), ncol = length(factors))
                    rownames(loadings_matrix) <- vars
                    colnames(loadings_matrix) <- factors
                    
                    for (j in seq_len(nrow(params))) {
                      if (params$op[j] == "=~") {
                        f_idx <- match(params$lhs[j], factors)
                        v_idx <- match(params$rhs[j], vars)
                        loadings_matrix[v_idx, f_idx] <- params$est[j]
                      }
                    }
                  }
                } else {
                  loadings_matrix <- lavaan::inspect(boot_efa, "est")$lambda
                }
                
                # S'assurer que les colonnes sont correctement nommées
                colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
                
                # Ajouter à la liste des résultats
                all_loadings[[i]] <- loadings_matrix
                
              }, error = function(e) {
                # Si une erreur se produit, ignorer cet échantillon
                message("Erreur dans l'échantillon bootstrap ", i, ": ", e$message)
              })
            }
            
            # Filtrer les échantillons qui ont échoué
            all_loadings <- all_loadings[!sapply(all_loadings, is.null)]
            n_successful <- length(all_loadings)
            
            if (n_successful < n_samples * 0.5) {
              showNotification(
                paste0("Attention: Seulement ", n_successful, " échantillons sur ", 
                       n_samples, " ont convergé. Les résultats peuvent être instables."),
                type = "warning",
                duration = 10
              )
            }
            
            incProgress(0.9, detail = "Calcul des statistiques...")
            
            # Initialiser des tableaux pour les statistiques
            var_names <- rownames(efa$loadings)
            factor_names <- colnames(efa$loadings)
            n_vars <- length(var_names)
            
            # Tableaux pour stocker les résultats
            mean_loadings <- matrix(0, nrow = n_vars, ncol = nfactors)
            sd_loadings <- matrix(0, nrow = n_vars, ncol = nfactors)
            ci_lower <- matrix(0, nrow = n_vars, ncol = nfactors)
            ci_upper <- matrix(0, nrow = n_vars, ncol = nfactors)
            
            rownames(mean_loadings) <- rownames(sd_loadings) <- 
              rownames(ci_lower) <- rownames(ci_upper) <- var_names
            colnames(mean_loadings) <- colnames(sd_loadings) <- 
              colnames(ci_lower) <- colnames(ci_upper) <- factor_names
            
            # Pour chaque variable et chaque facteur
            for (v in 1:n_vars) {
              for (f in 1:nfactors) {
                # Extraire toutes les charges factorielles pour cette variable/facteur
                loadings_vf <- sapply(all_loadings, function(l) l[v, f])
                
                # Calculer la moyenne et l'écart-type
                mean_loadings[v, f] <- mean(loadings_vf, na.rm = TRUE)
                sd_loadings[v, f] <- sd(loadings_vf, na.rm = TRUE)
                
                # Calculer les intervalles de confiance
                alpha <- (1 - ci_level) / 2
                ci_lower[v, f] <- quantile(loadings_vf, alpha, na.rm = TRUE)
                ci_upper[v, f] <- quantile(loadings_vf, 1 - alpha, na.rm = TRUE)
              }
            }
            
            # Stocker tous les résultats
            bootstrap_results(list(
              mean_loadings = mean_loadings,
              sd_loadings = sd_loadings,
              ci_lower = ci_lower,
              ci_upper = ci_upper,
              raw_loadings = all_loadings,
              n_samples = n_successful,
              ci_level = ci_level
            ))
            
            # Mettre à jour les sélecteurs de facteurs et variables
            updateSelectInput(session, "bootstrap_display_factor", 
                              choices = setNames(1:nfactors, paste0("Facteur ", 1:nfactors)),
                              selected = 1)
            
            updateSelectInput(session, "bootstrap_stability_var",
                              choices = var_names,
                              selected = var_names[1])
            
            showNotification(paste0("Analyse bootstrap terminée (", n_successful, " échantillons)"), 
                             type = "message")
            
          }, error = function(e) {
            showNotification(
              paste("Erreur lors de l'analyse bootstrap:", e$message),
              type = "error",
              duration = 10
            )
          })
        })
      })
      
      # Tableau des intervalles de confiance des charges
      output$bootstrap_loadings_table <- renderDT({
        req(bootstrap_results(), input$bootstrap_display_factor)
        
        results <- bootstrap_results()
        factor_idx <- as.numeric(input$bootstrap_display_factor)
        
        # Extraire les résultats pour le facteur sélectionné
        original_loadings <- efa_result()$loadings[, factor_idx]
        mean_loadings <- results$mean_loadings[, factor_idx]
        sd_loadings <- results$sd_loadings[, factor_idx]
        ci_lower <- results$ci_lower[, factor_idx]
        ci_upper <- results$ci_upper[, factor_idx]
        
        # Créer un tableau des résultats
        bootstrap_table <- data.frame(
          Variable = rownames(results$mean_loadings),
          Original = round(original_loadings, 4),
          Moyenne = round(mean_loadings, 4),
          "Écart-type" = round(sd_loadings, 4),
          "IC Inférieur" = round(ci_lower, 4),
          "IC Supérieur" = round(ci_upper, 4),
          check.names = FALSE
        )
        
        # Créer le tableau avec formatage
        datatable(
          bootstrap_table,
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            dom = 'ftip'
          ),
          rownames = FALSE
        ) %>%
          formatStyle(
            columns = "Original",
            fontWeight = "bold"
          ) %>%
          formatStyle(
            columns = c("IC Inférieur", "IC Supérieur"),
            backgroundColor = pal_efa$boot_bg
          ) %>%
          formatStyle(
            columns = "Moyenne",
            backgroundColor = styleInterval(
              c(0.3, 0.5, 0.7),
              c(pal_efa$loading_neutral, pal_efa$scale4_medium, pal_efa$scale4_good, pal_efa$scale4_excellent)
            ),
            color = "black"
          )
      })
      
      # Graphique de la stabilité des charges
      output$bootstrap_stability_plot <- renderPlotly({
        req(bootstrap_results(), input$bootstrap_stability_var)
        
        results <- bootstrap_results()
        selected_var <- input$bootstrap_stability_var
        var_idx <- match(selected_var, rownames(results$mean_loadings))
        
        # S'assurer que la variable existe
        if (is.na(var_idx)) {
          plot_ly() %>%
            add_annotations(
              text = "Variable non trouvée dans les résultats",
              showarrow = FALSE,
              font = list(size = 14)
            )
        }
        
        # Créer un dataframe avec les distributions des charges pour cette variable
        factor_names <- colnames(results$mean_loadings)
        n_factors <- length(factor_names)
        
        # Extraire toutes les charges bootstrappées pour cette variable
        all_samples <- results$raw_loadings
        n_samples <- length(all_samples)
        
        # Créer un dataframe pour stocker toutes les charges bootstrappées
        plot_data <- data.frame(
          Facteur = character(),
          Charge = numeric(),
          stringsAsFactors = FALSE
        )
        
        # Remplir le dataframe
        for (f in 1:n_factors) {
          factor_loadings <- sapply(all_samples, function(l) l[var_idx, f])
          plot_data <- rbind(
            plot_data,
            data.frame(
              Facteur = factor_names[f],
              Charge = factor_loadings
            )
          )
        }
        
        # Ajouter les valeurs originales pour référence
        original_loadings <- efa_result()$loadings[var_idx, ]
        
        # Créer le graphique violin plot
        p <- plot_ly() %>%
          add_trace(
            data = plot_data,
            x = ~Facteur,
            y = ~Charge,
            type = "violin",
            box = list(visible = TRUE),
            meanline = list(visible = TRUE),
            color = ~Facteur,
            points = "outliers"
          ) %>%
          add_trace(
            x = factor_names,
            y = original_loadings,
            type = "scatter",
            mode = "markers",
            marker = list(size = 10, symbol = "star", color = pal_efa$boot_marker),
            name = "Valeur originale"
          ) %>%
          layout(
            title = paste("Distribution des charges factorielles pour", selected_var),
            xaxis = list(title = "Facteur"),
            yaxis = list(title = "Charge factorielle", range = c(-1, 1)),
            showlegend = FALSE
          )
        
        return(p)
      })
      
      # Mettre à jour les contrôles du bootstrap quand l'analyse factorielle est exécutée
      observe({
        req(efa_result())
        
        # Obtenir les noms des variables et facteurs
        var_names <- rownames(efa_result()$loadings)
        n_factors <- ncol(efa_result()$loadings)
        
        # Mettre à jour les sélecteurs
        updateSelectInput(
          session,
          "bootstrap_stability_var",
          choices = var_names,
          selected = if(length(var_names) > 0) var_names[1] else NULL
        )
        
        updateSelectInput(
          session,
          "bootstrap_display_factor",
          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
          selected = if(n_factors > 0) 1 else NULL
        )
      })
      
      # Mettre à jour les sélecteurs après le bootstrap également
      observeEvent(bootstrap_results(), {
        req(bootstrap_results(), efa_result())
        
        results <- bootstrap_results()
        
        # Mettre à jour les sélecteurs de facteurs et variables
        var_names <- rownames(efa_result()$loadings)
        n_factors <- ncol(efa_result()$loadings)
        
        updateSelectInput(
          session,
          "bootstrap_display_factor", 
          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
          selected = 1
        )
        
        updateSelectInput(
          session,
          "bootstrap_stability_var",
          choices = var_names,
          selected = var_names[1]
        )
      })
      
      observeEvent(event_data("plotly_click"), {
        req(efa_result())
        click_data <- event_data("plotly_click")
        
        # Si on a des données de clic
        if (!is.null(click_data)) {
          # Récupérer la variable cliquée
          clicked_var <- click_data$customdata
          
          # Récupérer la liste actuelle des variables désactivées
          current_disabled <- input$disabled_variables
          if (is.null(current_disabled)) {
            current_disabled <- c()
          }
          
          # Basculer l'état de la variable (activer/désactiver)
          if (clicked_var %in% current_disabled) {
            # Si déjà désactivée, réactiver
            new_disabled <- setdiff(current_disabled, clicked_var)
          } else {
            # Si active, désactiver
            new_disabled <- c(current_disabled, clicked_var)
          }
          
          # Mettre à jour la liste des variables désactivées
          updateSelectizeInput(session, "disabled_variables", selected = new_disabled)
          
          # Si le recalcul automatique est activé, relancer l'analyse factorielle
          if (input$auto_recalculate) {
            # Déclencher le recalcul
            recalculate_efa_with_active_vars()
          }
        }
      })
      
      # Handler pour le bouton de recalcul manuel
      observeEvent(input$manual_recalculate, {
        recalculate_efa_with_active_vars()
      })
      
      # Ajouter cette fonction dans la partie serveur
      recalculate_efa_with_active_vars <- function() {
        req(filtered_data(), input$efa_variables, length(input$efa_variables) >= 2)
        
        # Obtenir les variables actives (celles qui ne sont pas désactivées)
        all_vars <- input$efa_variables
        disabled_vars <- input$disabled_variables
        if (is.null(disabled_vars)) {
          disabled_vars <- c()
        }
        active_vars <- setdiff(all_vars, disabled_vars)
        
        # Vérifier qu'il reste assez de variables pour l'analyse
        if (length(active_vars) < 2) {
          showNotification(
            "Il faut au moins 2 variables pour l'analyse factorielle",
            type = "error"
          )
          return()
        }
        
        withProgress(message = "Recalcul de l'analyse factorielle...", {
          tryCatch({
            # Sélection des variables pour l'analyse
            data_subset <- filtered_data()[, active_vars, drop = FALSE]
            
            # Vérifier les NA et les gérer
            data_complete <- na.omit(data_subset)
            
            # Exécuter l'analyse factorielle
            new_efa_model <- lavaan::efa(
              data = data_complete,
              nfactors = input$efa_nfactors,
              rotation = input$dynamic_rotation,
              estimator = input$efa_method,
              control = list(
                optim.method = "BFGS",
                optim.control = list(maxit = 1000),
                check.gradient = FALSE
              )
            )
            
            # Traiter les résultats comme dans la fonction run_efa
            model_summary <- summary(new_efa_model, nd = 3, cutoff = 0.001)
            
            if (inherits(new_efa_model, "efaList")) {
              loadings_matrix <- model_summary$loading
              
              if (is.null(loadings_matrix) || all(is.na(loadings_matrix))) {
                model <- new_efa_model[[1]]
                params <- lavaan::parameterEstimates(model)
                vars <- unique(params$rhs[params$op == "=~"])
                factors <- unique(params$lhs[params$op == "=~"])
                
                loadings_matrix <- matrix(0, nrow = length(vars), ncol = length(factors))
                rownames(loadings_matrix) <- vars
                colnames(loadings_matrix) <- factors
                
                for (i in seq_len(nrow(params))) {
                  if (params$op[i] == "=~") {
                    f_idx <- match(params$lhs[i], factors)
                    v_idx <- match(params$rhs[i], vars)
                    loadings_matrix[v_idx, f_idx] <- params$est[i]
                  }
                }
              }
              
              colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
              
              communalities <- model_summary$communalities
              if (is.null(communalities)) {
                communalities <- rowSums(loadings_matrix^2)
              }
              
              fit_indices <- lavaan::fitMeasures(new_efa_model[[1]])
              
            } else {
              loadings_matrix <- lavaan::inspect(new_efa_model, "est")$lambda
              colnames(loadings_matrix) <- paste0("F", 1:ncol(loadings_matrix))
              
              communalities <- rowSums(loadings_matrix^2)
              
              fit_indices <- lavaan::fitMeasures(new_efa_model)
            }
            
            uniqueness <- 1 - communalities
            names(communalities) <- rownames(loadings_matrix)
            names(uniqueness) <- rownames(loadings_matrix)
            
            variance_explained <- colSums(loadings_matrix^2) / nrow(loadings_matrix)
            cumulative_variance <- cumsum(variance_explained)
            
            
            # Mettre à jour l'objet des résultats
            efa_result(list(
              efa = new_efa_model,
              loadings = loadings_matrix,
              uniquenesses = uniqueness,
              communality = communalities,
              method = input$efa_method,
              rotation = input$dynamic_rotation,
              nfactors = input$efa_nfactors,
              Vaccounted = data.frame(
                SS_loadings = colSums(loadings_matrix^2),
                Proportion = variance_explained,
                Cumulative = cumulative_variance
              ),
              fit_indices = fit_indices,
              summary = model_summary
            ))
            
            showNotification("Analyse factorielle recalculée avec les variables actives", type = "message")
            
          }, error = function(e) {
            showNotification(
              paste("Erreur lors du recalcul de l'analyse factorielle:", e$message),
              type = "error",
              duration = 10
            )
          })
        })
      }
      
      # Observer qui met à jour les sélecteurs de facteurs lorsque l'analyse factorielle est exécutée
      observe({
        req(efa_result())
        efa <- efa_result()
        n_factors <- efa$nfactors
        
        # Mettre à jour les sélecteurs de facteurs
        updateSelectInput(session, "plot_factor_x", 
                          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
                          selected = 1)
        
        updateSelectInput(session, "plot_factor_y", 
                          choices = setNames(1:n_factors, paste0("Facteur ", 1:n_factors)),
                          selected = min(2, n_factors))
        
        # Mettre à jour les variables désactivées
        updateSelectizeInput(session, "disabled_variables", 
                             choices = rownames(efa$loadings), 
                             selected = character(0))
      })
      
      # Ajouter cette fonction pour générer une sortie HTML avec les informations sur l'analyse recalculée
      output$recalculated_info <- renderUI({
        req(efa_result())
        efa <- efa_result()
        
        # Créer un affichage HTML des informations clés
        html_content <- tags$div(
          tags$p(tags$b("Méthode d'extraction:"), efa$method),
          tags$p(tags$b("Méthode de rotation:"), efa$rotation),
          tags$p(tags$b("Nombre de facteurs:"), efa$nfactors),
          tags$p(tags$b("Variables actives:"), length(rownames(efa$loadings))),
          tags$hr(),
          tags$p(tags$b("Variance expliquée (%):"), style = "margin-bottom: 5px;")
        )
        
        # Ajouter un tableau pour la variance expliquée
        variance_table <- tags$table(
          style = "width: 100%; border-collapse: collapse;",
          tags$thead(
            tags$tr(
              tags$th(style = "border: 1px solid #ddd; padding: 4px; text-align: left;", "Facteur"),
              tags$th(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", "% Var"),
              tags$th(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", "% Cum")
            )
          ),
          tags$tbody(
            lapply(1:efa$nfactors, function(i) {
              tags$tr(
                tags$td(style = "border: 1px solid #ddd; padding: 4px;", paste0("F", i)),
                tags$td(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", 
                        sprintf("%.1f", efa$Vaccounted$Proportion[i] * 100)),
                tags$td(style = "border: 1px solid #ddd; padding: 4px; text-align: right;", 
                        sprintf("%.1f", efa$Vaccounted$Cumulative[i] * 100))
              )
            })
          )
        )
        
        # Ajouter un bouton pour afficher plus de détails
        detail_button <- tags$button(
          id = "show_detail_modal",
          class = "btn btn-sm btn-secondary mt-2",
          style = "width: 100%",
          onclick = "Shiny.setInputValue('show_detailed_results', Math.random(), {priority: 'event'});",
          "Afficher les détails complets"
        )
        
        return(tags$div(html_content, variance_table, detail_button))
      })
      
      # Ajouter ce gestionnaire d'événements pour afficher une fenêtre modale avec les résultats complets
      observeEvent(input$show_detailed_results, {
        req(efa_result())
        
        showModal(modalDialog(
          title = "Résultats détaillés de l'analyse factorielle",
          size = "l",
          easyClose = TRUE,
          
          tabsetPanel(
            tabPanel("Charges factorielles", 
                     DTOutput("modal_loadings_table")),
            tabPanel("Communalités", 
                     DTOutput("modal_communalities_table")),
            tabPanel("Indices d'ajustement", 
                     verbatimTextOutput("modal_fit_indices"))
          ),
          
          footer = modalButton("Fermer")
        ))
      })
      
      # Tableaux pour la fenêtre modale
      output$modal_loadings_table <- renderDT({
        req(efa_result(), input$loadings_threshold)
        loadings <- efa_result()$loadings
        threshold <- input$loadings_threshold
        
        # Formater pour l'affichage
        loadings_formatted <- round(loadings, 3)
        
        # Créer un tableau avec coloration (même convention que l'onglet Interprétation)
        datatable(
          loadings_formatted,
          options = list(pageLength = 15, scrollX = TRUE),
          rownames = TRUE
        ) %>%
          formatStyle(
            columns = names(loadings_formatted),
            backgroundColor = styleInterval(
              c(-threshold, threshold),
              c(pal_efa$loading_neutral, pal_efa$loading_neutral, pal_efa$loading_pos_bg)
            ),
            color = styleInterval(
              c(-threshold, threshold),
              c(pal_efa$loading_neg_text, "black", pal_efa$loading_pos_text)
            )
          )
      })
      
      output$modal_communalities_table <- renderDT({
        req(efa_result())
        
        comm <- efa_result()$communality
        uniq <- efa_result()$uniquenesses
        
        # Préparer les données pour le tableau
        comm_data <- data.frame(
          Variable = names(comm),
          Communalité = round(comm, 3),
          Unicité = round(uniq, 3)
        )
        
        # Créer le tableau
        datatable(
          comm_data,
          options = list(pageLength = 15, scrollX = TRUE),
          rownames = FALSE
        ) %>%
          formatStyle(
            "Communalité",
            background = styleColorBar(c(0, 1), pal_efa$comm_bar_modal),
            backgroundSize = "100% 90%",
            backgroundRepeat = "no-repeat",
            backgroundPosition = "center"
          )
      })
      
      output$modal_fit_indices <- renderPrint({
        req(efa_result())
        
        fit_indices <- efa_result()$fit_indices
        
        cat("INDICES D'AJUSTEMENT DU MODÈLE RECALCULÉ\n")
        cat("=====================================\n\n")
        
        # Sélectionner et afficher les indices les plus importants
        important_indices <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
        
        if (!is.null(fit_indices)) {
          for (idx in important_indices) {
            if (idx %in% names(fit_indices)) {
              cat(sprintf("%-10s: %8.4f\n", idx, fit_indices[idx]))
            }
          }
        } else {
          cat("Aucun indice d'ajustement disponible.")
        }
      })
      
      onSessionEnded(function() {
        # Nettoyer la mémoire lorsque l'utilisateur ferme l'application
        cleanup_reactive_objects()
        message("Session terminée - Mémoire libérée")
      })
      
      # Rendu des interprétations automatiques
      # Dans la partie rendu des interprétations automatiques, modifiez la section qui crée les cartes :
      
      output$auto_interpretation <- renderUI({
        req(efa_result())
        
        # Obtenir les charges factorielles et le seuil actuel
        loadings <- efa_result()$loadings
        threshold <- input$loadings_threshold
        
        # Générer les interprétations
        interpretations <- generate_factor_interpretation(loadings, threshold)
        
        # Créer des éléments d'interface pour chaque facteur
        factor_cards <- lapply(interpretations, function(interp) {
          # Obtenir la couleur appropriée depuis la palette centralisée
          quality_color <- pal_efa$quality[[interp$quality]]
          
          div(
            class = "card mb-3",
            div(
              class = "card-header",
              style = paste0("background-color: ", quality_color, "; color: white;"),
              h4(interp$name, style = "margin: 0;")
            ),
            div(
              class = "card-body",
              div(
                class = "alert", 
                style = paste0("background-color: ", quality_color, "20; color: black;"), # Ajout de "color: black"
                p(strong("Qualité de l'interprétation : "), style = "text-transform: uppercase;", interp$quality),
                p(interp$clarity)
              ),
              h5("Variables associées"),
              p(interp$pos_summary),
              if (length(interp$negative_vars) > 0) p(interp$neg_summary) else NULL,
              hr(),
              if (interp$n_vars >= 3) div(
                class = "alert alert-info",
                style = "color: black;", # Ajout de cette ligne pour garantir le texte en noir
                h5("Suggestion d'interprétation"),
                p(interp$suggestion),
                if (interp$n_vars >= 5) {
                  top_vars <- interp$positive_vars[1:min(5, length(interp$positive_vars))]
                  p("Variables clés à examiner : ", strong(paste(top_vars, collapse=", ")))
                }
              ) else div(
                class = "alert alert-warning",
                style = "color: black;", # Ajout de cette ligne pour garantir le texte en noir
                p("Nombre insuffisant de variables avec des charges significatives pour une interprétation fiable.")
              )
            )
          )
        })
        
        # Renvoyer tous les éléments
        div(factor_cards)
      })
      
      # Description des jeux de données classiques
      dataset_descriptions <- list(
        "HolzingerSwineford1939" = list(
          title = "Holzinger & Swineford (1939)",
          desc = "Données classiques de 301 étudiants sur 9 tests cognitifs. Souvent utilisées pour illustrer l'analyse factorielle en psychométrie. Source : package 'lavaan'.",
          variables = "Tests : x1-x3 (visuel), x4-x6 (textuel), x7-x9 (vitesse)",
          observations = "301 étudiants d'écoles du Midwest américain",
          source = "Package lavaan (Holzinger & Swineford, 1939)"
        ),
        "bfi_dataset" = list(
          title = "BFI Personality (25 items)",
          desc = "Inventaire des Cinq Facteurs (Big Five Inventory) mesurant les 5 grands traits de personnalité (Ouverture, Conscience, Extraversion, Agréabilité, Névrosisme).",
          variables = "25 items d'auto-évaluation de traits de personnalité",
          observations = "2800 participants",
          source = "Package psych (Revelle, W.)"
        ),
        "neo_n" = list(
          title = "NEO-PI-R (Neuroticism)",
          desc = "Données d'évaluation du facteur Neuroticism du NEO-PI-R, l'un des Big Five de la personnalité.",
          variables = "Items mesurant les 6 facettes du Neuroticism: Anxiété, Colère-Hostilité, Dépression, Conscience de soi, Impulsivité, Vulnérabilité",
          observations = "Environ 1000 participants",
          source = "Package EFA.dimensions (Costa & McCrae, 1992)"
        ),
        "neo_e" = list(
          title = "NEO-PI-R (Extraversion)",
          desc = "Données d'évaluation du facteur Extraversion du NEO-PI-R, l'un des Big Five de la personnalité.",
          variables = "Items mesurant les 6 facettes de l'Extraversion: Chaleur, Grégarité, Assertivité, Activité, Recherche d'excitation, Émotions positives",
          observations = "Environ 1000 participants",
          source = "Package EFA.dimensions (Costa & McCrae, 1992)"
        ),
        "neo_o" = list(
          title = "NEO-PI-R (Ouverture)",
          desc = "Données d'évaluation du facteur Ouverture du NEO-PI-R, l'un des Big Five de la personnalité.",
          variables = "Items mesurant les 6 facettes de l'Ouverture: Fantaisie, Esthétique, Sentiments, Actions, Idées, Valeurs",
          observations = "Environ 1000 participants",
          source = "Package EFA.dimensions (Costa & McCrae, 1992)"
        ),
        "msceit" = list(
          title = "Intelligence Émotionnelle (MSCEIT)",
          desc = "Données du test d'intelligence émotionnelle MSCEIT (Mayer-Salovey-Caruso Emotional Intelligence Test).",
          variables = "Items mesurant la perception, l'utilisation, la compréhension et la gestion des émotions",
          observations = "Échantillon représentatif d'adultes",
          source = "Package EFA.dimensions (Mayer, Salovey & Caruso, 2002)"
        ),
        "bfi2" = list(
          title = "BFI-2 (15 facets)",
          desc = "Big Five Inventory-2, une mesure hiérarchique de la personnalité avec 15 facettes (3 par facteur).",
          variables = "60 items mesurant 15 facettes des 5 grands facteurs de personnalité",
          observations = "Échantillon représentatif d'adultes",
          source = "Package EFAfactors (Soto & John, 2017)"
        ),
        "spi" = list(
          title = "SAPA Personality Inventory",
          desc = "Inventaire de personnalité issu du projet SAPA (Synthetic Aperture Personality Assessment).",
          variables = "Items mesurant divers traits de personnalité",
          observations = "Grand échantillon en ligne",
          source = "Package EFAfactors (Condon & Revelle, 2015)"
        ),
        "holzinger" = list(
          title = "Holzinger (1939) - Tests cognitifs",
          desc = "Données classiques de tests cognitifs, version disponible dans EFAtools.",
          variables = "9 tests cognitifs visuels, textuels et de vitesse",
          observations = "301 étudiants",
          source = "Package EFAtools (Holzinger & Swineford, 1939)"
        ),
        "WIT" = list(
          title = "WIT data",
          desc = "Données d'un test d'intelligence avec 14 variables.",
          variables = "14 mesures d'aptitudes cognitives diverses",
          observations = "Échantillon d'adultes",
          source = "Package EFAtools"
        ),
        "box_data" = list(
          title = "Thurstone (1947) - Box dataset",
          desc = "Données classiques de Thurstone utilisées dans les analyses factorielles.",
          variables = "Variables mesurant différentes aptitudes",
          observations = "Échantillon standardisé",
          source = "Package EFAutilities (Thurstone, 1947)"
        ),
        "wais_data" = list(
          title = "Jöreskog (2003) - WAIS",
          desc = "Données du test d'intelligence de Wechsler (WAIS) utilisées par Jöreskog.",
          variables = "Sous-tests du WAIS mesurant différentes composantes de l'intelligence",
          observations = "Échantillon clinique standardisé",
          source = "Package EFAutilities (Jöreskog, 2003)"
        )
      )
      
      # Afficher les informations sur le jeu de données sélectionné
      output$dataset_info <- renderUI({
        req(input$classic_dataset)
        
        if (input$classic_dataset == "") {
          return(NULL)
        }
        
        info <- dataset_descriptions[[input$classic_dataset]]
        
        tagList(
          h5(info$title, style = "margin-top: 0;"),
          p(info$desc),
          p(tags$b("Variables :"), info$variables),
          p(tags$b("Observations :"), info$observations),
          p(tags$b("Source :"), info$source)
        )
      })
      
      # Chargement d'un jeu de données classique
      observeEvent(input$load_classic, {
        req(input$classic_dataset, input$classic_dataset != "")
        
        # Afficher une notification pendant le chargement
        id <- showNotification("Chargement du jeu de données...", type = "message", duration = NULL)
        
        tryCatch({
          # Charger le jeu de données sélectionné
          data <- switch(input$classic_dataset,
                         "HolzingerSwineford1939" = {
                           # Assurez-vous que lavaan est chargé
                           require(lavaan)
                           # Charger les données
                           dat <- lavaan::HolzingerSwineford1939
                           # Sélectionner uniquement les variables de test (x1 à x9)
                           dat[, c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")]
                         },
                         "bfi_dataset" = {
                           require(psych)
                           data(bfi)
                           # Sélectionner uniquement les items, pas les données démographiques
                           bfi[, 1:25]
                         },
                         "neo_n" = {
                           require(EFA.dimensions)
                           data(neo_n)
                           neo_n
                         },
                         "neo_e" = {
                           require(EFA.dimensions)
                           data(neo_e)
                           neo_e
                         },
                         "neo_o" = {
                           require(EFA.dimensions)
                           data(neo_o)
                           neo_o
                         },
                         "msceit" = {
                           require(EFA.dimensions)
                           data(msceit)
                           msceit
                         },
                         "bfi2" = {
                           require(EFAfactors)
                           data(bfi2)
                           bfi2
                         },
                         "spi" = {
                           require(EFAfactors)
                           data(spi)
                           spi
                         },
                         "holzinger" = {
                           require(EFAtools)
                           data(holzinger)
                           holzinger$data
                         },
                         "WIT" = {
                           require(EFAtools)
                           data(WIT)
                           WIT
                         },
                         "box_data" = {
                           require(EFAutilities)
                           data(box)
                           box
                         },
                         "wais_data" = {
                           require(EFAutilities)
                           data(wais)
                           wais
                         },
                         {
                           # Si aucune correspondance
                           stop("Jeu de données non disponible")
                         }
          )
          
          # Stocker les données
          raw_data(data)
          filtered_data(data)
          
          # Mettre à jour les sélecteurs
          updateSelectInput(session, "filter_var", 
                            choices = c("Aucun" = "", names(data)))
          
          numeric_cols <- names(data)[sapply(data, is.numeric)]
          updateSelectInput(session, "efa_variables", 
                            choices = numeric_cols,
                            selected = numeric_cols)
          
          # Mettre à jour les sélecteurs pour la visualisation
          updateSelectInput(session, "viz_variables", 
                            choices = numeric_cols,
                            selected = numeric_cols[1:min(5, length(numeric_cols))])
          
          # Afficher l'onglet d'analyse
          showTab(inputId = "main_tabs", target = "Analyse factorielle")
          
          # Replier "Importation des données" et ouvrir "Paramètres de l'analyse"
          accordion_panel_set(id = "sidebar_accordion", values = "Paramètres de l'analyse")
          
          removeNotification(id)
          showNotification(paste("Données classiques chargées :", input$classic_dataset), 
                           type = "message")
          
        }, error = function(e) {
          # Gérer les erreurs
          removeNotification(id)
          showNotification(
            paste("Erreur lors du chargement du jeu de données:", e$message),
            type = "error",
            duration = 10
          )
        })
      })

} 

shinyApp(ui, server)