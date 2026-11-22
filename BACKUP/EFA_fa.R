
library(bslib)
library(psych)
library(DT)
library(ggplot2)
library(reshape2)
library(corrplot)
library(plotly)
library(readxl)
library(haven)
library(car)
library(nortest)
library(GGally)
library(EGAnet)
library(openxlsx)
library(flextable)

# Interface utilisateur
ui <- page_fluid(
  theme = bs_theme(bootswatch = "flatly"),
  
  # Titre de l'application
  h1("Analyse Factorielle Exploratoire"),
  
  # Panneau d'onglets principal
  navset_card_tab(
    # Nouvel onglet : Introduction
    nav_panel(
      title = "Introduction",
      card(
        layout_columns(
          col_widths = c(12),
          card(
            card_header("Présentation de l'analyse factorielle exploratoire"),
            div(
              style = "padding: 15px;",
              h4("Qu'est-ce que l'analyse factorielle exploratoire ?"),
              p("L'analyse factorielle exploratoire (AFE) est une technique statistique utilisée pour identifier les structures latentes (facteurs) sous-jacentes à un ensemble de variables observées. Contrairement à l'analyse factorielle confirmatoire, qui teste une structure préétablie, l'AFE vise à découvrir cette structure."),
              
              h4("Objectifs de l'analyse factorielle"),
              tags$ul(
                tags$li("Réduire un grand nombre de variables à un ensemble plus petit de facteurs"),
                tags$li("Identifier les dimensions latentes ou construits sous-jacents aux variables observées"),
                tags$li("Examiner les relations entre variables et générer des hypothèses sur ces relations"),
                tags$li("Évaluer la validité de construit d'un instrument de mesure")
              ),
              
              h4("Équations de base"),
              p("L'analyse factorielle est basée sur le modèle suivant :"),
              div(
                style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; font-family: monospace;",
                HTML("X<sub>i</sub> = &lambda;<sub>i1</sub>F<sub>1</sub> + &lambda;<sub>i2</sub>F<sub>2</sub> + ... + &lambda;<sub>im</sub>F<sub>m</sub> + &epsilon;<sub>i</sub>")
              ),
              p("Où :"),
              tags$ul(
                tags$li(HTML("X<sub>i</sub> est la i<sup>ème</sup> variable observée")),
                tags$li(HTML("F<sub>j</sub> est le j<sup>ème</sup> facteur commun")),
                tags$li(HTML("&lambda;<sub>ij</sub> est la charge factorielle de la i<sup>ème</sup> variable sur le j<sup>ème</sup> facteur")),
                tags$li(HTML("&epsilon;<sub>i</sub> est l'erreur unique (ou spécificité) associée à la i<sup>ème</sup> variable"))
              ),
              
              p("En notation matricielle :"),
              div(
                style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; font-family: monospace;",
                HTML("X = &Lambda;F + &epsilon;")
              ),
              
              h4("Étapes de l'analyse factorielle exploratoire"),
              div(
                style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px;",
                tags$ol(
                  tags$li(strong("Préparation des données"), 
                          p("Vérifier la taille de l'échantillon, examiner les distributions, gérer les valeurs manquantes.")),
                  tags$li(strong("Évaluation de l'adéquation des données"), 
                          p("Utiliser des tests comme celui de Kaiser-Meyer-Olkin (KMO) et le test de sphéricité de Bartlett.")),
                  tags$li(strong("Extraction des facteurs"), 
                          p("Appliquer une méthode d'extraction comme l'analyse en composantes principales, le maximum de vraisemblance, etc.")),
                  tags$li(strong("Détermination du nombre de facteurs"), 
                          p("Utiliser des critères comme la règle de Kaiser, le test du coude, l'analyse parallèle, ou le MAP test.")),
                  tags$li(strong("Rotation des facteurs"), 
                          p("Appliquer une rotation orthogonale (ex. Varimax) ou oblique (ex. Oblimin) pour faciliter l'interprétation.")),
                  tags$li(strong("Interprétation des facteurs"), 
                          p("Examiner les charges factorielles et donner un sens aux facteurs identifiés.")),
                  tags$li(strong("Validation de la solution factorielle"), 
                          p("Évaluer la qualité de l'ajustement, la fiabilité et la validité de la solution."))
                )
              ),
              
              h4("Différence entre l'Analyse en Composantes Principales (ACP) et l'Analyse Factorielle (AF)"),
              div(
                style = "display: flex; margin-top: 15px;",
                div(
                  style = "flex: 1; background-color: #e9f7fe; padding: 15px; border-radius: 5px; margin-right: 10px;",
                  h5("Analyse en Composantes Principales"),
                  tags$ul(
                    tags$li("Vise à expliquer la variance totale des variables"),
                    tags$li("Ne fait pas de distinction entre variance commune et variance unique"),
                    tags$li("Utile pour la réduction de dimensionnalité"),
                    tags$li("Suppose que toute la variance est commune")
                  )
                ),
                div(
                  style = "flex: 1; background-color: #f0f0e0; padding: 15px; border-radius: 5px;",
                  h5("Analyse Factorielle"),
                  tags$ul(
                    tags$li("Vise à expliquer les corrélations entre variables"),
                    tags$li("Distingue entre variance commune et variance unique"),
                    tags$li("Utile pour identifier des facteurs latents"),
                    tags$li("Suppose l'existence de facteurs communs sous-jacents")
                  )
                )
              )
            )
          ),
          
          card(
            card_header("Utilisation de cette application"),
            div(
              style = "padding: 15px;",
              p("Cette application vous guide à travers les étapes de l'analyse factorielle exploratoire :"),
              tags$ol(
                tags$li(strong("Onglet Données"), 
                        p("Chargez vos données brutes ou une matrice de corrélations.")),
                tags$li(strong("Onglet Paramètres"), 
                        p("Configurez votre analyse et explorez différentes méthodes pour déterminer le nombre optimal de facteurs.")),
                tags$li(strong("Onglet Résultats"), 
                        p("Examinez en détail les résultats de votre analyse factorielle."))
              ),
              
              h4("Types de fichiers supportés"),
              tags$ul(
                tags$li(strong("CSV, TXT"), " - Fichiers texte délimités"),
                tags$li(strong("XLSX, XLS"), " - Fichiers Excel"),
                tags$li(strong("SPSS (.sav)"), " - Fichiers SPSS"),
                tags$li(strong("STATA (.dta)"), " - Fichiers Stata")
              ),
              
              h4("Conseils d'utilisation"),
              tags$ul(
                tags$li("Commencez par explorer différentes méthodes pour déterminer le nombre optimal de facteurs (test parallèle, MAP, EGA)"),
                tags$li("Expérimentez avec différentes méthodes d'extraction et de rotation pour trouver la solution la plus interprétable"),
                tags$li("Examinez les charges factorielles et les communalités pour évaluer la qualité de votre solution"),
                tags$li("Vérifiez les indices d'ajustement pour confirmer l'adéquation de votre modèle")
              )
            )
          ),
          # Informations sur les méthodes (déplacées depuis l'onglet Paramètres)
          card(
            card_header("Informations sur les méthodes"),
            tabsetPanel(
              tabPanel("Extraction", 
                       div(
                         style = "padding: 15px;",
                         h4("Méthodes d'extraction factorielle"),
                         p("Le choix de la méthode d'extraction affecte comment les facteurs sont extraits des données:"),
                         tags$ul(
                           tags$li(strong("Minimum Residual (minres)"), ": Méthode robuste qui minimise les résidus. Recommandée comme option par défaut."),
                           tags$li(strong("Maximum Likelihood (ml)"), ": Estime les charges qui maximisent la probabilité d'observer la matrice de corrélation. Nécessite la normalité multivariée."),
                           tags$li(strong("Principal Axis Factoring (pa)"), ": Méthode itérative basée sur les communalités. Alternative robuste quand les données ne sont pas normales."),
                           tags$li(strong("Unweighted Least Squares (uls)"), ": Minimise la somme des carrés des différences entre les corrélations observées et reproduites."),
                           tags$li(strong("Generalized Least Squares (gls)"), ": Similaire à ULS mais pondère les corrélations inversement à leur unicité."),
                           tags$li(strong("Minimum Rank (minrank)"), ": Minimise le rang de la matrice de résidus."),
                           tags$li(strong("Alpha factoring (alpha)"), ": Maximise la fiabilité alpha des facteurs.")
                         )
                       )
              ),
              tabPanel("Rotation", 
                       div(
                         style = "padding: 15px;",
                         h4("Méthodes de rotation factorielle"),
                         p("Les rotations aident à simplifier l'interprétation des facteurs:"),
                         h5("Rotations orthogonales (facteurs non corrélés)"),
                         tags$ul(
                           tags$li(strong("Varimax"), ": Maximise la variance des charges au carré pour chaque facteur. La plus utilisée des rotations orthogonales."),
                           tags$li(strong("Quartimax"), ": Simplifie la structure des variables en maximisant la variance des charges au carré pour chaque variable."),
                           tags$li(strong("Equamax"), ": Compromis entre Varimax et Quartimax.")
                         ),
                         h5("Rotations obliques (facteurs pouvant être corrélés)"),
                         tags$ul(
                           tags$li(strong("Promax"), ": Rotation oblique rapide qui augmente l'interprétabilité des facteurs."),
                           tags$li(strong("Oblimin"), ": Produit une structure simple avec des facteurs corrélés."),
                           tags$li(strong("Simplimax"), ": Vise à maximiser le nombre de charges nulles.")
                         )
                       )
              ),
              tabPanel("Sélection", 
                       div(
                         style = "padding: 15px;",
                         h4("Critères pour déterminer le nombre de facteurs"),
                         tags$ul(
                           tags$li(strong("Règle de Kaiser"), ": Conserver les facteurs avec valeurs propres > 1."),
                           tags$li(strong("Test du coude (Scree test)"), ": Examiner le graphique des valeurs propres et identifier le 'coude'."),
                           tags$li(strong("Analyse parallèle"), ": Compare les valeurs propres observées avec celles de données aléatoires."),
                           tags$li(strong("Test MAP de Velicer"), ": Examine les corrélations partielles après extraction des facteurs. Le nombre optimal de facteurs correspond au minimum du critère MAP."),
                           tags$li(strong("Analyse Graphique Exploratoire (EGA)"), ": Utilise des techniques de détection de communautés dans les réseaux pour identifier le nombre de dimensions."),
                           tags$li(strong("Variance expliquée"), ": Conserver suffisamment de facteurs pour expliquer un pourcentage acceptable de la variance (souvent 60-80%)."),
                           tags$li(strong("Interprétabilité"), ": Les facteurs doivent avoir un sens théorique.")
                         ),
                         hr(),
                         h5("Analyse parallèle"),
                         p("Le", strong("test parallèle"), "est une méthode robuste pour déterminer le nombre optimal de facteurs à extraire. Il compare les valeurs propres de vos données avec celles obtenues à partir de données aléatoires ayant les mêmes caractéristiques. Les facteurs à retenir sont ceux dont les valeurs propres sont supérieures à celles des données aléatoires."),
                         hr(),
                         h5("Test MAP de Velicer"),
                         p("Le", strong("test MAP (Minimum Average Partial)"), "de Velicer est une méthode statistique basée sur l'analyse des corrélations partielles. Après extraction de chaque facteur, on calcule la moyenne des corrélations partielles au carré. Le nombre optimal de facteurs correspond au point où cette moyenne est minimale, indiquant que l'extraction de facteurs supplémentaires n'améliore plus l'explication des relations entre variables."),
                         hr(),
                         h5("Analyse Graphique Exploratoire (EGA)"),
                         p("L'", strong("Analyse Graphique Exploratoire"), "est une méthode récente qui utilise la théorie des réseaux pour déterminer le nombre de dimensions dans les données. Elle représente les variables comme des nœuds dans un réseau, et les corrélations comme des connexions entre ces nœuds. Des algorithmes de détection de communautés identifient ensuite des groupes de variables fortement interconnectées, correspondant à des facteurs latents."),
                         p("L'indice", strong("TEFI (Total Entropy Fit Index)"), "évalue la qualité de l'ajustement du modèle EGA, notamment en ce qui concerne la structure de communautés (ou de dimensions/facteurs) détectée dans le graphe des corrélations ou des associations entre variables. TEFI n'a pas de seuil universel: il est interprété comparativement entre plusieurs modèles. Une faible valeur indique un meilleur ajustement.")
                       )
              )
            )
          ),
          card(
            card_header("Ressources et références"),
            div(
              style = "padding: 15px;",
              h4("Ouvrages de référence"),
              tags$ul(
                tags$li("Brown, T. A. (2015). Confirmatory factor analysis for applied research (2nd ed.). Guilford Press."),
                tags$li("Fabrigar, L. R., & Wegener, D. T. (2011). Exploratory factor analysis. Oxford University Press."),
                tags$li("Field, A. (2018). Discovering Statistics Using IBM SPSS Statistics (5th ed.). SAGE Publications."),
                tags$li("Hair, J. F., Black, W. C., Babin, B. J., & Anderson, R. E. (2019). Multivariate data analysis (8th ed.). Cengage Learning."),
                tags$li("Tabachnick, B. G., & Fidell, L. S. (2019). Using multivariate statistics (7th ed.). Pearson.")
              ),
              
              h4("Articles scientifiques importants"),
              tags$ul(
                tags$li("Costello, A. B., & Osborne, J. W. (2005). Best practices in exploratory factor analysis: Four recommendations for getting the most from your analysis. Practical Assessment, Research & Evaluation, 10(7), 1-9."),
                tags$li("Fabrigar, L. R., Wegener, D. T., MacCallum, R. C., & Strahan, E. J. (1999). Evaluating the use of exploratory factor analysis in psychological research. Psychological Methods, 4(3), 272-299."),
                tags$li("Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30(2), 179-185."),
                tags$li("Velicer, W. F. (1976). Determining the number of components from the matrix of partial correlations. Psychometrika, 41(3), 321-327.")
              ),
              
              h4("Ressources en ligne"),
              tags$ul(
                tags$li(tags$a(href = "https://stats.idre.ucla.edu/spss/seminars/introduction-to-factor-analysis/", "Introduction à l'analyse factorielle (UCLA)")),
                tags$li(tags$a(href = "https://personality-project.org/r/psych/", "Package psych (William Revelle)")),
                tags$li(tags$a(href = "https://www.r-bloggers.com/2019/05/principal-component-analysis-in-r/", "ACP avec R (R-bloggers)")),
                tags$li(tags$a(href = "https://www.youtube.com/playlist?list=PLkk92zzyru5OAtc_ItUubaSSq6S_TGmCP", "Tutoriels vidéo sur l'analyse factorielle")),
                tags$li(tags$a(href = "https://egatutorial.com/", "Tutoriel sur l'analyse graphique exploratoire (EGA)")),
                tags$li(tags$a(href = "https://cran.r-project.org/web/packages/EGAnet/vignettes/EGA.html", "Documentation du package EGAnet")),
                tags$li(tags$a(href = "https://www.rdocumentation.org/packages/corrplot/versions/0.92/topics/corrplot", "Documentation du package corrplot")),
                tags$li(tags$a(href = "https://cran.r-project.org/web/packages/lavaan/vignettes/tutorial.pdf", "Tutoriel lavaan pour les modèles de mesure"))
              )
            )
          ),
          
          card(
            card_header("À propos de cette application"),
            div(
              style = "padding: 15px;",
              p("Cette application a été développée pour faciliter l'exploration et l'analyse factorielle de données multivariées."),
              p("Auteur : Daniel Coulombe, Ph.D"),
              p("Version : 1.0.0"),
              p("Dernière mise à jour : ", Sys.Date()),
              hr(),
              p("Cette application utilise les packages R suivants :"),
              tags$ul(
                tags$li(strong("psych"), " - Pour l'analyse factorielle et les tests associés"),
                tags$li(strong("corrplot"), " - Pour la visualisation des matrices de corrélation"),
                tags$li(strong("EGAnet"), " - Pour l'analyse graphique exploratoire"),
                tags$li(strong("bslib"), " - Pour l'interface utilisateur moderne"),
                tags$li(strong("plotly"), " - Pour les graphiques interactifs")
              ),
              p("Pour citer cette application :"),
              div(
                style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px; font-family: monospace; font-size: 0.9em;",
                "Daniel Coulombe. (", format(Sys.Date(), "%Y"), "). Analyse Factorielle Exploratoire [Application Shiny]. URL: [votre-url]"
              )
            )
          )
        )
      )
    ),
    
    # Premier onglet : Chargement des données
    nav_panel(
      title = "Données",
      layout_sidebar(
        sidebar = sidebar(
          width = 350,
          card(
            card_header("Chargement des données"),
            radioButtons("data_type", "Type de données :", 
                         choices = c("Données brutes" = "raw", 
                                     "Matrice de corrélations" = "cormat")),
            conditionalPanel(
              condition = "input.data_type == 'raw'",
              fileInput("rawfile", "Charger les données :",
                        accept = c(
                          "text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt",
                          ".xlsx", ".xls",           # Excel
                          ".sav",                    # SPSS
                          ".dta"                     # Stata
                        )),
              conditionalPanel(
                condition = "input.rawfile && (input.rawfile.name.endsWith('.csv') || input.rawfile.name.endsWith('.txt'))",
                checkboxInput("header", "En-têtes", TRUE),
                radioButtons("sep", "Séparateur :",
                             choices = c(Virgule = ",", "Point-virgule" = ";", Tabulation = "\t"),
                             selected = ","),
                radioButtons("dec", "Décimal :",
                             choices = c(Point = ".", Virgule = ","),
                             selected = ".")
              ),
              conditionalPanel(
                condition = "output.data_loaded",
                uiOutput("var_selection_ui")
              )
            ),
            conditionalPanel(
              condition = "input.data_type == 'cormat'",
              fileInput("corfile", "Charger la matrice (CSV, TXT) :",
                        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")),
              numericInput("n_obs", "Nombre d'observations :", 100, min = 10),
              # Ajouter l'UI pour la sélection des variables quand une matrice est chargée
              conditionalPanel(
                condition = "output.cormat_loaded",
                uiOutput("cormat_var_selection_ui")
              )
            )
          )
        ),
        card(
          card_header("Aperçu des données"),
          conditionalPanel(
            condition = "output.data_loaded",
            uiOutput("data_tabs")
          ),
          conditionalPanel(
            condition = "!output.data_loaded",
            p("Aucune donnée chargée.")
          )
        ),

          conditionalPanel(
            condition = "!output.data_loaded",
            div(
              style = "padding: 30px; text-align: center; color: #888;",
              h4("Aucune donnée chargée"),
              p("Veuillez d'abord charger des données dans l'onglet 'Données' pour voir les tests d'adéquation.")
            )
          )
        
      )
    ),
    
    # Deuxième onglet : Paramètres de l'analyse factorielle
    nav_panel(
      title = "Paramètres",
      layout_sidebar(
        sidebar = sidebar(
          width = 300,
          card(
            card_header("Paramètres de l'analyse"),
            # Choix du type d'analyse
            radioButtons("analysis_type", "Type d'analyse :", 
                         choices = c("Analyse Factorielle (AF)" = "fa", 
                                     "Composantes Principales (ACP)" = "pca")),
            # Options du modèle
            numericInput("nfactors", "Nombre de facteurs :", 3, min = 1, max = 20),
            selectInput("fm", "Méthode d'extraction :", 
                        choices = c("Minimum Residual" = "minres",
                                    "Maximum Likelihood" = "ml",
                                    "Principal Axis Factoring" = "pa",
                                    "Unweighted Least Squares" = "uls",
                                    "Generalized Least Squares" = "gls",
                                    "Minimum Rank" = "minrank",
                                    "Alpha factoring" = "alpha")),
            selectInput("rotate", "Méthode de rotation :", 
                        choices = c("Aucune" = "none",
                                    "Varimax" = "varimax",
                                    "Quartimax" = "quartimax", 
                                    "Equamax" = "equamax",
                                    "Promax" = "promax",
                                    "Oblimin" = "oblimin",
                                    "Simplimax" = "simplimax"),
                        selected = "oblimin"),
            checkboxInput("kaiser", "Normalisation de Kaiser", TRUE),
            conditionalPanel(
              condition = "input.rotate == 'promax'",
              numericInput("promax_k", "Exposant Promax (k) :", 4, min = 2, max = 10)
            ),
            selectInput("cor", "Matrice analysée :", 
                        choices = c("Corrélation de Pearson" = "cor",
                                    "Corrélation de Spearman" = "spearman",
                                    "Covariance" = "cov"),
                        selected = "cor"),
            checkboxInput("use_smc", "Utiliser SMC comme estimations initiales", TRUE),
            numericInput("max_iter", "Nombre max d'itérations :", 100, min = 10, max = 1000),
            actionButton("run_analysis", "Lancer l'analyse", class = "btn-primary")
          )
        ),
        
        # Zone principale avec contenu informatif
  #      layout_columns(
  #        col_widths = c(4, 8),
          
          # Colonne 2: Aides visuelles pour déterminer le nombre de facteurs
          card(
            card_header("Aide au choix du nombre de facteurs"),
            conditionalPanel(
              condition = "output.data_loaded && input.data_type == 'raw'",
              tabsetPanel(
                id = "factor_choice_tabs",
                tabPanel("Scree Plot", 
                         layout_columns(
                           col_widths = c(7, 5), # Répartition : graphique à gauche (7) et eigenvalues à droite (5)
                           
                           # Card pour le graphique Scree Plot (côté gauche)
                           card(
                             card_header("Graphique des valeurs propres (Scree Plot)"),
                             div(
                               class = "p-3",
                               plotOutput("initial_scree_plot", height = "400px")
                             ),
                             div(
                               style = "padding: 10px; font-style: italic; font-size: 90%;",
                               "Le graphique des valeurs propres (scree plot) basé sur les données brutes. La ligne horizontale rouge représente la règle de Kaiser (valeurs propres > 1) et la ligne verticale verte indique le nombre de facteurs retenus."
                             )
                           ),
                           
                           # Card pour la liste des eigenvalues (côté droit)
                           card(
                             card_header("Valeurs propres (Eigenvalues)"),
                             div(
                               class = "p-3",
                               DTOutput("eigenvalues_table_raw")
                             )
                           )
                         )
                ),
                tabPanel("Critère de Kaiser", 
                         tableOutput("kaiser_criterion_table"),
                         div(
                           style = "padding: 10px; font-style: italic; font-size: 90%;",
                           "La règle de Kaiser suggère de conserver les facteurs dont les valeurs propres sont supérieures à 1. Ce tableau montre les valeurs propres des facteurs potentiels et indique combien devraient être retenus selon ce critère."
                         )
                ),
                tabPanel("Test MAP", 
                         plotOutput("map_test_plot", height = "300px"),
                         verbatimTextOutput("map_test_result"),
                         div(
                           style = "padding: 10px; font-style: italic; font-size: 90%;",
                           "Le test MAP (Minimum Average Partial) de Velicer est une méthode statistique pour déterminer le nombre optimal de facteurs. Il examine les corrélations partielles après avoir retiré successivement l'influence de chaque facteur."
                         )
                ),
                tabPanel("Test parallèle", 
                         div(
                           style = "padding: 10px; display: flex; justify-content: space-between; align-items: center;",
                           div(
                             style = "flex-grow: 1;",
                             p("Le test parallèle compare les valeurs propres observées avec celles obtenues à partir de données aléatoires.")
                           ),
                           div(
                             style = "width: 280px;",
                             numericInput("n_datasets_auto", "Jeux de données simulés:", 100, min = 20, max = 1000),
                             numericInput("percentile_auto", "Percentile:", 95, min = 50, max = 100)
                           )
                         ),
                         plotOutput("parallel_plot_auto", height = "300px"),
                         verbatimTextOutput("parallel_summary_auto"),
#                         actionButton("run_parallel_auto", "Actualiser le test parallèle", icon = icon("refresh"))
                ),

tabPanel("EGA",
         layout_sidebar(
           sidebar = sidebar(
             width = 300,
             # Mettre tous les contrôles dans la sidebar
             h4("Paramètres EGA"),
             radioButtons("ega_model", "Modèle d'estimation:", 
                          choices = c("GLASSO" = "glasso", 
                                     "TMFG" = "TMFG"),
                          selected = "glasso",
                          inline = TRUE),
             radioButtons("ega_algorithm", "Algorithme de détection:", 
                          choices = c("Walktrap" = "walktrap", 
                                     "Louvain" = "louvain",
                                     "Leiden" = "leiden"),
                          selected = "louvain",
                          inline = TRUE),

             h4("Options graphiques"),
             selectInput("ega_layout", "Layout du graphique:", 
                         choices = c("Spring" = "spring", 
                                   "Cercle" = "circle",
                                   "MDS" = "mds",
                                   "Fruchterman-Reingold" = "fruchtermanreingold",
                                   "Kamada-Kawai" = "kamadakawai"),
                         selected = "kamadakawai"),
             sliderInput("ega_node_size", "Taille des nœuds:", 15, min = 5, max = 25),
             numericInput("ega_label_size", "Taille des étiquettes:", 4, min = 1, max = 10),
             sliderInput("ega_edge_alpha", "Transparence:", min = 0, max = 1, value = 0.8, step = 0.1)
           ),
           # Zone principale - Tabset avec sous-onglets
           tabsetPanel(
             id = "ega_inner_tabs",
             
             # Sous-onglet du graphique
             tabPanel("Graphique",
                      plotOutput("ega_plot", height = "600px")
             ),
             
             # Sous-onglet des résultats
             tabPanel("Résultats",
                      uiOutput("ega_summary")
             ),
             
             # Sous-onglet pour les Network Loads
             tabPanel("Network Loads",
                      uiOutput("netloads_ui")
             ),
             
             # Nouveau sous-onglet pour la matrice de réseau
             tabPanel("Matrice de réseau",
                      div(
                        style = "padding: 10px;",
                        p("Cette matrice représente les relations (arêtes) entre les variables (nœuds) dans le réseau EGA."),
                        p("Les valeurs plus élevées indiquent des relations plus fortes entre les variables.")
                      ),
                      DTOutput("ega_network_matrix")
             )
           )
         )
)


              )
              ),


conditionalPanel(
  condition = "output.data_loaded && input.data_type == 'cormat'",
  tabsetPanel(
    tabPanel("Scree Plot", 
             layout_columns(
               col_widths = c(7, 5), # Répartition : graphique à gauche (7) et eigenvalues à droite (5)
               
               # Card pour le graphique Scree Plot (côté gauche)
               card(
                 card_header("Graphique des valeurs propres (Scree Plot)"),
                 div(
                   class = "p-3",
                   plotOutput("initial_scree_plot_cormat", height = "400px")
                 ),
                 div(
                   style = "padding: 10px; font-style: italic; font-size: 90%;",
                   "Le graphique des valeurs propres (scree plot) basé sur la matrice de corrélation chargée. La ligne horizontale rouge représente la règle de Kaiser (valeurs propres > 1) et la ligne verticale verte indique le nombre de facteurs retenus."
                 )
               ),
               
               # Card pour la liste des eigenvalues (côté droit)
               card(
                 card_header("Valeurs propres (Eigenvalues)"),
                 div(
                   class = "p-3",
                   DTOutput("eigenvalues_table_cormat")
                 )
               )
             )
    ),
    tabPanel("Critère de Kaiser", 
             tableOutput("kaiser_criterion_table_cormat"),
             div(
               style = "padding: 10px; font-style: italic; font-size: 90%;",
               "La règle de Kaiser suggère de conserver les facteurs dont les valeurs propres sont supérieures à 1. Ce tableau montre les valeurs propres des facteurs potentiels et indique combien devraient être retenus selon ce critère."
             )
    ),
    tabPanel("Test MAP", 
             plotOutput("map_test_plot_cormat", height = "300px"),
             verbatimTextOutput("map_test_result_cormat"),
             div(
               style = "padding: 10px; font-style: italic; font-size: 90%;",
               "Le test MAP (Minimum Average Partial) de Velicer est une méthode statistique pour déterminer le nombre optimal de facteurs. Il examine les corrélations partielles après avoir retiré successivement l'influence de chaque facteur."
             )
    ),
    tabPanel("Test parallèle", 
             div(
               style = "padding: 10px; display: flex; justify-content: space-between; align-items: center;",
               div(
                 style = "flex-grow: 1;",
                 p("Le test parallèle compare les valeurs propres observées avec celles obtenues à partir de données aléatoires.")
               ),
               div(
                 style = "width: 280px;",
                 numericInput("n_datasets_auto_cormat", "Jeux de données simulés:", 100, min = 20, max = 1000),
                 numericInput("percentile_auto_cormat", "Percentile:", 95, min = 50, max = 100)
               )
             ),
             plotOutput("parallel_plot_auto_cormat", height = "300px"),
             verbatimTextOutput("parallel_summary_auto_cormat")
    ),
    tabPanel("EGA",
             layout_sidebar(
               sidebar = sidebar(
                 width = 300,
                 # Mettre tous les contrôles dans la sidebar
                 h4("Paramètres EGA"),
                 radioButtons("ega_model_cormat", "Modèle d'estimation:", 
                              choices = c("GLASSO" = "glasso", 
                                          "TMFG" = "TMFG"),
                              selected = "glasso",
                              inline = TRUE),
                 radioButtons("ega_algorithm_cormat", "Algorithme de détection:", 
                              choices = c("Walktrap" = "walktrap", 
                                          "Louvain" = "louvain",
                                          "Leiden" = "leiden"),
                              selected = "louvain",
                              inline = TRUE),
                 
                 h4("Options graphiques"),
                 selectInput("ega_layout_cormat", "Layout du graphique:", 
                             choices = c("Spring" = "spring", 
                                         "Cercle" = "circle",
                                         "MDS" = "mds",
                                         "Fruchterman-Reingold" = "fruchtermanreingold",
                                         "Kamada-Kawai" = "kamadakawai"),
                             selected = "kamadakawai"),
                 sliderInput("ega_node_size_cormat", "Taille des nœuds:", 15, min = 5, max = 25),
                 numericInput("ega_label_size_cormat", "Taille des étiquettes:", 4, min = 1, max = 10),
                 sliderInput("ega_edge_alpha_cormat", "Transparence:", min = 0, max = 1, value = 0.8, step = 0.1)
               ),
               # Zone principale - Tabset avec sous-onglets
               tabsetPanel(
                 id = "ega_inner_tabs_cormat",
                 
                 # Sous-onglet du graphique
                 tabPanel("Graphique",
                          plotOutput("ega_plot_cormat", height = "600px")
                 ),
                 
                 # Sous-onglet des résultats
                 tabPanel("Résultats",
                          uiOutput("ega_summary_cormat")
                 ),
                 
                 # Sous-onglet pour les Network Loads
                 tabPanel("Network Loads",
                          uiOutput("netloads_ui_cormat")
                 ),
                 
                 # Nouveau sous-onglet pour la matrice de réseau
                 tabPanel("Matrice de réseau",
                          div(
                            style = "padding: 10px;",
                            p("Cette matrice représente les relations (arêtes) entre les variables (nœuds) dans le réseau EGA."),
                            p("Les valeurs plus élevées indiquent des relations plus fortes entre les variables.")
                          ),
                          DTOutput("ega_network_matrix_cormat")
                 )
               )
             )
    )
  )
),
            conditionalPanel(
              condition = "!output.data_loaded",
              div(
                style = "padding: 30px; text-align: center; color: #888;",
                h4("Aucune donnée chargée"),
                p("Veuillez d'abord charger des données dans l'onglet 'Données' pour voir les outils d'aide au choix du nombre de facteurs.")
              )
            )
          )
        )
        
      
    ),
    
# Troisième onglet : Résultats
nav_panel(
  title = "Résultats",
  conditionalPanel(
    condition = "output.analysis_done",
    uiOutput("results_content")  # Utiliser un uiOutput dynamique ici
  ),
  conditionalPanel(
    condition = "!output.analysis_done",
    p("Analyse factorielle non exécutée. Veuillez aller à l'onglet Paramètres et lancer l'analyse.")
  )
)
)
  )



# Logique serveur
server <- function(input, output, session) {
  
  # Réactifs pour les données
  data_raw <- reactive({
    req(input$data_type == "raw")
    req(input$rawfile)
    
    # Récupérer l'extension du fichier
    file_ext <- tools::file_ext(input$rawfile$name)
    
    # Charger les données en fonction de l'extension
    if (file_ext %in% c("csv", "txt")) {
      # Fichiers CSV/TXT
      tryCatch({
        data <- read.table(input$rawfile$datapath, 
                           header = input$header, 
                           sep = input$sep, 
                           dec = input$dec,
                           stringsAsFactors = FALSE)
      }, error = function(e) {
        showNotification(paste("Erreur lors du chargement du fichier CSV:", e$message), type = "error")
        return(NULL)
      })
    } else if (file_ext %in% c("xlsx", "xls")) {
      # Fichiers Excel
      tryCatch({
        data <- readxl::read_excel(input$rawfile$datapath)
        data <- as.data.frame(data)  # Convertir en data.frame standard
      }, error = function(e) {
        showNotification(paste("Erreur lors du chargement du fichier Excel:", e$message), type = "error")
        return(NULL)
      })
    } else if (file_ext == "sav") {
      # Fichiers SPSS
      tryCatch({
        data <- haven::read_spss(input$rawfile$datapath)
        data <- as.data.frame(data)  # Convertir en data.frame standard
        
        # Conserver les valeurs numériques pour les variables ordinales
        # On conserve les valeurs brutes des variables labellisées au lieu de les convertir en facteurs
        for (col in names(data)) {
          if (inherits(data[[col]], "haven_labelled")) {
            # Conserver les valeurs numériques telles quelles
            # On extrait simplement la valeur numérique de chaque élément labellisé
            data[[col]] <- as.numeric(data[[col]])
          }
        }
      }, error = function(e) {
        showNotification(paste("Erreur lors du chargement du fichier SPSS:", e$message), type = "error")
        return(NULL)
      })
    } else if (file_ext == "dta") {
      # Fichiers Stata
      tryCatch({
        data <- haven::read_dta(input$rawfile$datapath)
        data <- as.data.frame(data)  # Convertir en data.frame standard
        
        # Conserver les valeurs numériques pour les variables ordinales
        for (col in names(data)) {
          if (inherits(data[[col]], "haven_labelled")) {
            data[[col]] <- as.numeric(data[[col]])
          }
        }
      }, error = function(e) {
        showNotification(paste("Erreur lors du chargement du fichier Stata:", e$message), type = "error")
        return(NULL)
      })
    }
    
    return(data)
  })
  
  # Fonction pour obtenir la matrice de corrélation appropriée selon le type de données
  correlation_matrix_for_analysis <- reactive({
    if (input$data_type == "raw") {
      req(analysis_data())
      
      # Calculer la matrice de corrélation selon le type spécifié
      if (input$cor == "cor") {
        return(cor(analysis_data(), use = "pairwise.complete.obs"))
      } else if (input$cor == "spearman") {
        return(cor(analysis_data(), use = "pairwise.complete.obs", method = "spearman"))
      } else if (input$cor == "cov") {
        return(cov(analysis_data(), use = "pairwise.complete.obs"))
      }
    } else {
      req(cormat_data())
      req(selected_cormat_vars())
      
      return(cormat_data()[selected_cormat_vars(), selected_cormat_vars(), drop = FALSE])
    }
  })
  
  # Fonction pour obtenir le nombre d'observations approprié
  n_observations <- reactive({
    if (input$data_type == "raw") {
      req(analysis_data())
      return(nrow(analysis_data()))
    } else {
      req(input$n_obs)
      return(input$n_obs)
    }
  })
  
  # Remplacer le code actuel de output$data_tabs par:
  output$data_tabs <- renderUI({
    if (input$data_type == "raw") {
      # Pour les données brutes, afficher tous les onglets
      tabsetPanel(
        tabPanel("Données", DTOutput("data_preview")),
        tabPanel("Statistiques", 
                 layout_columns(
                   col_widths = c(7, 5),
                   card(uiOutput("summary")),
                   card(
                     card_header("Tests de normalité"),
                     selectInput("normality_test", "Type de test:",
                                 choices = c("Shapiro-Wilk" = "shapiro",
                                             "Kolmogorov-Smirnov" = "ks",
                                             "Anderson-Darling" = "ad"),
                                 selected = "shapiro"),
                     DTOutput("normality_results")
                   )
                 )
        ),
        tabPanel("Corrélations", 
                 layout_columns(
                   col_widths = c(2, 10),  # Réduire la largeur de la colonne des contrôles
                   card(
                     card_header("Options d'affichage"),
                     div(
                       style = "padding: 10px;",
                       selectInput("corr_method", "Représentation:", 
                                   choices = c("Cercle" = "circle", 
                                               "Carré" = "square",
                                               "Ellipse" = "ellipse", 
                                               "Nombre" = "number",
                                               "Pie" = "pie",
                                               "Couleur" = "color"),
                                   selected = "ellipse"),
                       selectInput("corr_order", "Ordre:", 
                                   choices = c("Original" = "original", 
                                               "AOE" = "AOE",
                                               "FPC" = "FPC", 
                                               "Hiérarchique" = "hclust",
                                               "Alphabet" = "alphabet"),
                                   selected = "hclust"),
                       selectInput("corr_color_palette", "Palette de couleurs:", 
                                   choices = c("Divergente (RdBu)" = "RdBu",
                                               "Divergente (RdYlBu)" = "RdYlBu",
                                               "Divergente (PiYG)" = "PiYG",
                                               "Séquent. (Blues)" = "Blues",
                                               "Séquent. (Reds)" = "Reds"),
                                   selected = "RdYlBu"),
                       checkboxInput("corr_show_numbers", "Afficher les coefficients", TRUE),
                       numericInput("corr_number_digits", "Décimales:", 2, min = 0, max = 4, width = "75%"),
                       hr(),
                       downloadButton("download_cormat", "Télécharger la matrice", class = "btn-sm btn-info w-100")
                     )
                   ),
                   card(
                     card_header("Matrice de corrélation"),
                     plotOutput("corplot", height = "800px")
                   )
                 )
        ),

        # Nouvel onglet pour les tests d'adéquation
        tabPanel("Adéquation des données",
                 layout_columns(
                   col_widths = c(6, 3),
                   card(
                     card_header("Tests d'adéquation principaux"),
                     div(
                       class = "p-3",
                       popover(
                         actionButton("info_adequacy", "À propos des tests d'adéquation", 
                                      icon = icon("info-circle"), 
                                      class = "btn-sm btn-outline-info"),
                         div(
                           p("L'indice KMO mesure l'adéquation de l'échantillonnage. Valeurs conventionnelles:"),
                           tags$ul(
                             tags$li("0.9 - 1.0 : Excellente"),
                             tags$li("0.8 - 0.9 : Bonne"),
                             tags$li("0.7 - 0.8 : Moyenne"),
                             tags$li("0.6 - 0.7 : Médiocre"),
                             tags$li("0.5 - 0.6 : Faible"),
                             tags$li("< 0.5 : Inacceptable")
                           ),
                           hr(),
                           p("Le test de Bartlett vérifie si la matrice de corrélation est différente d'une matrice identité."),
                           p("Une p-valeur < 0.05 indique que les corrélations entre variables sont significativement différentes de zéro.")
                         ),
                         title = "Tests d'adéquation pour l'analyse factorielle",
                         placement = "right",
                         options = list(html = TRUE, width = "400px")
                       ),
                       DTOutput("adequacy_summary_table"),
                       hr(),
                       DTOutput("bartlett_result")
                     )
                   ),
                   # Card pour les MSA par variable (côté droit)
                   card(
                     card_header("Mesures MSA par variable"),
                     div(
                       class = "p-3",
                       # Ajouter un bouton d'information pour les MSA
                       popover(
                         actionButton("info_msa", "À propos des mesures MSA", 
                                      icon = icon("info-circle"), 
                                      class = "btn-sm btn-outline-info"),
                         div(
                           p("Les mesures MSA (Measure of Sampling Adequacy) indiquent l'adéquation de chaque variable pour l'analyse factorielle:"),
                           tags$ul(
                             tags$li("0.9 - 1.0 : Excellente"),
                             tags$li("0.8 - 0.9 : Bonne"),
                             tags$li("0.7 - 0.8 : Moyenne"),
                             tags$li("0.6 - 0.7 : Médiocre"),
                             tags$li("0.5 - 0.6 : Faible"),
                             tags$li("< 0.5 : Inacceptable - Considérer l'exclusion de cette variable")
                           ),
                           p("Une valeur MSA < 0.5 suggère que la variable peut ne pas convenir à l'analyse factorielle.")
                         ),
                         title = "Mesures MSA par variable",
                         placement = "right",
                         options = list(html = TRUE, width = "400px")
                       ),
                       br(),
                       # Tableau MSA (sans graphique)
                       DTOutput("msa_table")
                     )
                   )
                 )
        )
      )
    } else {
      # Pour les matrices de corrélations, afficher uniquement l'onglet Corrélations
      tabsetPanel(
        tabPanel("Corrélations", 
                 layout_columns(
                   col_widths = c(2, 10),  # Réduire la largeur de la colonne des contrôles
                   card(
                     card_header("Options d'affichage"),
                     div(
                       style = "padding: 10px;",
                       selectInput("corr_method", "Représentation:", 
                                   choices = c("Cercle" = "circle", 
                                               "Carré" = "square",
                                               "Ellipse" = "ellipse", 
                                               "Nombre" = "number",
                                               "Pie" = "pie",
                                               "Couleur" = "color"),
                                   selected = "ellipse"),
                       selectInput("corr_order", "Ordre:", 
                                   choices = c("Original" = "original", 
                                               "AOE" = "AOE",
                                               "FPC" = "FPC", 
                                               "Hiérarchique" = "hclust",
                                               "Alphabet" = "alphabet"),
                                   selected = "hclust"),
                       selectInput("corr_color_palette", "Palette de couleurs:", 
                                   choices = c("Divergente (RdBu)" = "RdBu",
                                               "Divergente (RdYlBu)" = "RdYlBu",
                                               "Divergente (PiYG)" = "PiYG",
                                               "Séquent. (Blues)" = "Blues",
                                               "Séquent. (Reds)" = "Reds"),
                                   selected = "RdYlBu"),
                       checkboxInput("corr_show_numbers", "Afficher les coefficients", TRUE),
                       numericInput("corr_number_digits", "Décimales:", 2, min = 0, max = 4, width = "75%")
                     )
                   ),
                   card(
                     card_header("Matrice de corrélation"),
                     plotOutput("corplot", height = "800px")
                   )
                 )
        )
      )
    }
  })
  
  # Onglet Adéquation des données - Modification pour éviter la redondance du test de Bartlett
  output$adequacy_summary_table <- renderDT({
    req(adequacy_tests_data())
    
    # Extraire les résultats
    results <- adequacy_tests_data()
    
    # Interpréter l'indice KMO
    kmo_interp <- function(kmo) {
      if (kmo >= 0.9) return("Excellente")
      if (kmo >= 0.8) return("Bonne")
      if (kmo >= 0.7) return("Moyenne")
      if (kmo >= 0.6) return("Médiocre")
      if (kmo >= 0.5) return("Faible")
      return("Inacceptable")
    }
    
    # Créer le dataframe pour le tableau - SUPPRESSION DU TEST DE BARTLETT
    adequacy_df <- data.frame(
      Test = c("Kaiser-Meyer-Olkin (KMO)"),
      Valeur = c(sprintf("%.3f", results$kmo)),
      Interprétation = c(paste0("Adéquation ", tolower(kmo_interp(results$kmo)))),
      Statut = c(ifelse(results$kmo >= 0.5, "Acceptable", "Non acceptable"))
    )
    
    # Créer le datatable avec formatage de couleur
    datatable(adequacy_df,
              options = list(dom = 't', ordering = FALSE, paging = FALSE),
              rownames = FALSE) %>%
      formatStyle(
        'Statut',
        backgroundColor = styleEqual(
          c("Acceptable", "Non acceptable"),
          c('#d4edda', '#f8d7da')
        ),
        fontWeight = 'bold'
      ) %>%
      formatStyle(
        'Test',
        fontWeight = 'bold'
      )
  })
  
  # Données pour l'analyse
  analysis_data <- reactive({
    if (input$data_type == "raw") {
      # Si données brutes, prendre le sous-ensemble sélectionné (code existant)
      req(data_raw())
      req(selected_vars())
      
      # Vérifier qu'au moins 3 variables sont sélectionnées
      req(length(selected_vars()) >= 3)
      
      # Extraire seulement les colonnes numériques sélectionnées
      numeric_data <- data_raw()[, selected_vars(), drop = FALSE]
      
      # Vérifier que toutes les variables sont numériques
      non_numeric <- sapply(numeric_data, function(x) !is.numeric(x))
      if (any(non_numeric)) {
        stop("Toutes les variables sélectionnées doivent être numériques pour l'analyse factorielle")
      }
      
      return(numeric_data)
    } else {
      # Si matrice de corrélations, appliquer la sélection des variables
      req(cormat_data())
      req(selected_cormat_vars())
      
      # Vérifier qu'au moins 3 variables sont sélectionnées
      req(length(selected_cormat_vars()) >= 3)
      
      # Extraire le sous-ensemble de la matrice
      cormat_subset <- cormat_data()[selected_cormat_vars(), selected_cormat_vars(), drop = FALSE]
      
      return(cormat_subset)
    }
  })
  
  # Interface pour sélectionner les variables
  output$var_selection_ui <- renderUI({
    req(data_raw())
    
    # Identifier les colonnes numériques
    numeric_cols <- sapply(data_raw(), is.numeric)
    
    if (sum(numeric_cols) < 3) {
      return(p("Pas assez de variables numériques dans les données (minimum 3 nécessaires)"))
    }
    
    # Liste des variables numériques
    numeric_vars <- names(data_raw())[numeric_cols]
    
    tagList(
      div(
        style = "margin-bottom: 10px;",
        strong("Sélectionner les variables :"),
        div(
          style = "display: flex; justify-content: space-between; margin-top: 5px; margin-bottom: 10px;",
          actionButton("select_all_vars", "Tout sélectionner", class = "btn-sm btn-info"),
          actionButton("deselect_all_vars", "Tout désélectionner", class = "btn-sm btn-outline-secondary")
        )
      ),
      # Au lieu d'utiliser deux colonnes de checkboxes, utiliser un seul groupe
      # mais les afficher sur deux colonnes avec CSS
      div(
        style = "column-count: 2; column-gap: 10px;",
        checkboxGroupInput("selected_vars", 
                           label = NULL,
                           choices = numeric_vars,
                           selected = numeric_vars)
      )
    )
  })
  
  selected_vars <- reactive({
    input$selected_vars
  })
  
  selected_cormat_vars <- reactive({
    input$selected_cormat_vars
  })
  
  # Mettre à jour les observateurs
  observeEvent(input$select_all_vars, {
    req(data_raw())
    numeric_cols <- sapply(data_raw(), is.numeric)
    numeric_vars <- names(data_raw())[numeric_cols]
    updateCheckboxGroupInput(session, "selected_vars", selected = numeric_vars)
  })
  
  observeEvent(input$deselect_all_vars, {
    updateCheckboxGroupInput(session, "selected_vars", selected = character(0))
  })
  
  observeEvent(input$select_all_cormat_vars, {
    req(cormat_data())
    vars <- colnames(cormat_data())
    updateCheckboxGroupInput(session, "selected_cormat_vars", selected = vars)
  })
  
  observeEvent(input$deselect_all_cormat_vars, {
    updateCheckboxGroupInput(session, "selected_cormat_vars", selected = character(0))
  })
  
  # Indicateurs pour les sorties conditionnelles
  output$data_loaded <- reactive({
    if (input$data_type == "raw") {
      return(!is.null(data_raw()))
    } else {
      return(!is.null(cormat_data()))
    }
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  # Aperçu des données
  output$data_preview <- renderDT({
    # Vérifier si des données sont disponibles
    if (input$data_type == "raw") {
      req(data_raw())
      data <- data_raw()
    } else {
      req(cormat_data()) # Remplacer data_cormat() par cormat_data()
      data <- cormat_data() # Remplacer data_cormat() par cormat_data()
    }
    
    # Formater les nombres à 4 décimales pour l'affichage
    datatable(data, options = list(
      pageLength = 10,
      autoWidth = TRUE
    )) %>% 
      formatRound(columns = which(sapply(data, is.numeric)), digits = 4)
  })
  
  # Statistiques descriptives
  output$summary <- renderUI({
    if (input$data_type == "raw") {
      req(analysis_data())
      
      # Calculer les statistiques descriptives avec psych::describe pour les variables sélectionnées
      stats_df <- psych::describe(analysis_data())
      
      # Convertir en data.frame et ajuster pour flextable
      stats_df <- as.data.frame(stats_df)
      stats_df$variable <- rownames(stats_df)
      rownames(stats_df) <- NULL
      
      # Réorganiser les colonnes pour avoir la variable en premier
      stats_df <- stats_df[, c("variable", "n", "mean", "sd", "median", "min", "max", "skew", "kurtosis")]
      
      # Renommer les colonnes pour plus de clarté
      names(stats_df) <- c("Variable", "N", "Moyenne", "Écart-type", "Médiane", "Minimum", "Maximum", "Asymétrie", "Aplatissement")
      
      # Créer le flextable
      ft <- flextable(stats_df) %>%
        colformat_double(j = 3:9, digits = 3) %>%
        theme_booktabs() %>%
        autofit() %>%
        add_header_lines("Statistiques descriptives des variables sélectionnées")
      
      # Appliquer des styles
      ft <- ft %>%
        bg(bg = "#f5f5f5", part = "header") %>%
        bold(part = "header")
      
      # Retourner le flextable dans un élément HTML
      flextable::htmltools_value(ft)
      
    } else {
      req(cormat_data())
      req(selected_cormat_vars())
      
      # Afficher des informations sur les variables sélectionnées dans la matrice de corrélation
      tagList(
        h4("Matrice de corrélations chargée"),
        p(paste("Variables sélectionnées:", length(selected_cormat_vars()))),
        p(paste("Nombre d'observations spécifié:", input$n_obs))
      )
    }
  })
  
  
  # Tests de normalité
  output$normality_results <- renderDT({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Préparer le tableau des résultats
    results <- data.frame(
      Variable = character(),
      'Test Statistic' = numeric(),
      'p-value' = numeric(),
      'Normality' = character(),
      stringsAsFactors = FALSE
    )
    
    # Réaliser les tests de normalité pour chaque variable
    for (col in names(analysis_data())) {
      x <- analysis_data()[[col]]
      x <- x[!is.na(x)] # Supprimer les valeurs manquantes
      
      if (length(x) < 3) next # Ignorer les variables avec trop peu de données
      
      if (input$normality_test == "shapiro") {
        # Test de Shapiro-Wilk (limité à 5000 observations)
        if (length(x) > 5000) x <- sample(x, 5000)
        test <- shapiro.test(x)
        stat <- test$statistic
        p_val <- test$p.value
      } else if (input$normality_test == "ks") {
        # Test de Kolmogorov-Smirnov
        test <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x))
        stat <- test$statistic
        p_val <- test$p.value
      } else if (input$normality_test == "ad") {
        # Test d'Anderson-Darling
        if (requireNamespace("nortest", quietly = TRUE)) {
          test <- nortest::ad.test(x)
          stat <- test$statistic
          p_val <- test$p.value
        } else {
          # Fallback au Shapiro-Wilk si nortest n'est pas disponible
          if (length(x) > 5000) x <- sample(x, 5000)
          test <- shapiro.test(x)
          stat <- test$statistic
          p_val <- test$p.value
        }
      }
      
      # Interprétation du test
      normality <- ifelse(p_val > 0.05, "Distribution normale", "Non normal")
      
      # Ajouter les résultats au tableau
      results <- rbind(results, data.frame(
        Variable = col,
        'Test Statistic' = round(as.numeric(stat), 4),
        'p-value' = round(p_val, 4),
        'Normality' = normality,
        stringsAsFactors = FALSE
      ))
    }
    
    # Formatter et afficher le tableau
    datatable(results, options = list(
      pageLength = 10,
      autoWidth = TRUE
    )) %>%
      formatStyle('Normality',
                  backgroundColor = styleEqual(
                    c("Distribution normale", "Non normal"),
                    c('#d4edda', '#f8d7da')
                  ))
  })
  
  cormat_data <- reactive({
    req(input$data_type == "cormat")
    req(input$corfile)
    
    # Charger la matrice de corrélations
    mat <- read.table(input$corfile$datapath, 
                      header = TRUE, 
                      sep = ",", 
                      row.names = 1)
    
    # Vérifier que c'est une matrice carrée
    if (nrow(mat) != ncol(mat)) {
      stop("La matrice chargée n'est pas carrée (nombre de lignes ≠ nombre de colonnes)")
    }
    
    # Vérifier la symétrie
    if (!all(abs(mat - t(mat)) < 1e-10)) {
      stop("La matrice chargée n'est pas symétrique")
    }
    
    return(mat)
  })
  
  # Mettre à jour les choix de variables pour la matrice de diagrammes de dispersion
  observe({
    req(input$data_type == "raw")
    req(data_raw())
    
    data <- data_raw()
    num_vars <- names(data)[sapply(data, is.numeric)]
    
    updateSelectInput(session, "scatter_vars", 
                      choices = num_vars,
                      selected = if(length(num_vars) > 6) num_vars[1:6] else num_vars)
  })
  
  # Générer la matrice de diagrammes de dispersion
  output$scatter_matrix <- renderPlot({
    req(input$data_type == "raw")
    req(data_raw())
    req(input$scatter_vars)
    
    # Limiter à 6 variables maximum pour des raisons de performance et de lisibilité
    selected_vars <- input$scatter_vars
    if(length(selected_vars) > 6) {
      selected_vars <- selected_vars[1:6]
      showNotification("Limité à 6 variables pour la matrice de dispersion.", type = "warning")
    }
    
    data <- data_raw()[, selected_vars, drop = FALSE]
    
    # Configuration pour le panneau diagonal
    diag_func <- input$scatter_diag_type
    
    # Configuration pour le panneau inférieur
    lower_func <- switch(input$scatter_lower_type,
                         "points" = wrap("points", alpha = input$scatter_alpha),
                         "smooth" = wrap("smooth_loess", 
                                         alpha = input$scatter_alpha, 
                                         color = "blue",           # Couleur de la ligne de lissage
                                         fill = "lightblue",       # Couleur de la bande de confiance
                                         size = 1)
                         )
    
    # Configuration pour le panneau supérieur - coefficients uniquement
    upper_func <- if(input$scatter_show_corr) {
      function(data, mapping, ...) {
        x <- eval_data_col(data, mapping$x)
        y <- eval_data_col(data, mapping$y)
        corr <- cor(x, y, use = "pairwise.complete.obs")
        
        # Créer un graphique vide avec seulement le coefficient de corrélation
        ggplot(data = data, mapping = mapping) + 
          annotate("text", x = mean(range(x, na.rm = TRUE)), y = mean(range(y, na.rm = TRUE)), 
                   label = sprintf(paste0("%.", input$scatter_corr_digits, "f"), corr), 
                   size = input$scatter_corr_size) +
          theme_void()
      }
    } else {
      "blank"
    }
    
    # Créer la matrice de diagrammes de dispersion avec GGally::ggpairs
    ggpairs(
      data,
      upper = list(continuous = upper_func),
      lower = list(continuous = lower_func),
      diag = list(continuous = diag_func),
      axisLabels = "show"
    ) + 
      theme_bw() + 
      theme(
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 9),
        panel.grid.minor = element_blank()
      )
  })
  
  output$cormat_loaded <- reactive({
    return(!is.null(cormat_data()))
  })
  outputOptions(output, "cormat_loaded", suspendWhenHidden = FALSE)
  
  # Ajouter l'UI pour sélectionner les variables de la matrice de corrélation :
  output$cormat_var_selection_ui <- renderUI({
    req(cormat_data())
    
    # Liste des variables dans la matrice de corrélation
    vars <- colnames(cormat_data())
    
    tagList(
      div(
        style = "margin-bottom: 10px;",
        strong("Sélectionner les variables :"),
        div(
          style = "display: flex; justify-content: space-between; margin-top: 5px; margin-bottom: 10px;",
          actionButton("select_all_cormat_vars", "Tout sélectionner", class = "btn-sm btn-info"),
          actionButton("deselect_all_cormat_vars", "Tout désélectionner", class = "btn-sm btn-outline-secondary")
        )
      ),
      # Utiliser un seul groupe de checkboxes affiché sur deux colonnes
      div(
        style = "column-count: 2; column-gap: 10px;",
        checkboxGroupInput("selected_cormat_vars", 
                           label = NULL,
                           choices = vars,
                           selected = vars)
      )
    )
  })
  
  
  # Ajouter les observateurs pour les boutons de sélection/désélection
  observeEvent(input$select_all_cormat_vars, {
    req(cormat_data())
    vars <- colnames(cormat_data())
    
    # Séparer les variables en deux groupes pour les deux colonnes
    vars_left <- vars[seq(1, length(vars), by = 2)]
    vars_right <- vars[seq(2, length(vars), by = 2)]
    
    updateCheckboxGroupInput(session, "selected_cormat_vars_left", selected = vars_left)
    updateCheckboxGroupInput(session, "selected_cormat_vars_right", selected = vars_right)
  })
  
  observeEvent(input$deselect_all_cormat_vars, {
    updateCheckboxGroupInput(session, "selected_cormat_vars_left", selected = character(0))
    updateCheckboxGroupInput(session, "selected_cormat_vars_right", selected = character(0))
  })
  
  
  # Matrice de corrélation
  correlation_matrix <- reactive({
    if (input$data_type == "raw") {
      req(analysis_data())
      # Calculer la matrice de corrélation selon le type spécifié
      if (input$cor == "cor") {
        cor(analysis_data(), use = "pairwise.complete.obs")
      } else if (input$cor == "spearman") {
        cor(analysis_data(), use = "pairwise.complete.obs", method = "spearman")
      } else if (input$cor == "cov") {
        cov(analysis_data(), use = "pairwise.complete.obs")
      }
    } else {
      req(cormat_data())
      req(selected_cormat_vars())
      
      # Utiliser le sous-ensemble sélectionné de la matrice
      cormat_data()[selected_cormat_vars(), selected_cormat_vars(), drop = FALSE]
    }
  })
  
  # Plot de la matrice de corrélation
  output$corplot <- renderPlot({
    # Utiliser la fonction réactive correlation_matrix() déjà définie
    # qui contient la matrice de corrélation des variables sélectionnées
    req(correlation_matrix())
    corr_to_plot <- correlation_matrix()
    
    # S'assurer que corr_to_plot est bien une matrice et non une liste
    if (is.data.frame(corr_to_plot)) {
      corr_to_plot <- as.matrix(corr_to_plot)
    }
    
    # Créer la visualisation
    corrplot(corr_to_plot, 
             method = input$corr_method,
             order = input$corr_order,
             type = "lower",
             diag = FALSE,
             tl.col = "black",
             tl.srt = 45,
             col = colorRampPalette(RColorBrewer::brewer.pal(9, input$corr_color_palette))(200),
             addCoef.col = ifelse(input$corr_show_numbers, "black", "transparent"),
             number.digits = input$corr_number_digits,
             number.cex = 0.7)
  })
  
  # Téléchargement de la matrice de corrélations
  output$download_cormat <- downloadHandler(
    filename = function() {
      paste("correlation_matrix_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(correlation_matrix(), file)
    }
  )
  
  # Scree plot initial pour les données brutes (avant l'analyse factorielle)
  output$initial_scree_plot <- renderPlot({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Calculer la matrice de corrélation
    if (input$cor == "cor") {
      cor_mat <- cor(analysis_data(), use = "pairwise.complete.obs", method = "pearson")
    } else if (input$cor == "spearman") {
      cor_mat <- cor(analysis_data(), use = "pairwise.complete.obs", method = "spearman")
    } else {
      cor_mat <- cov(analysis_data(), use = "pairwise.complete.obs")
    }
    
    # Calcul des valeurs propres
    eigenvalues <- eigen(cor_mat)$values
    
    # Déterminer le nombre de facteurs selon la règle de Kaiser ou la sélection de l'utilisateur
    n_factors <- sum(eigenvalues >= 1)  # Nombre selon le critère de Kaiser
    
    # Si aucun facteur n'est au-dessus de 1, prendre au moins 1 facteur
    if (n_factors == 0) n_factors <- 1
    
    # Utiliser le nombre de facteurs défini par l'utilisateur s'il est disponible
    if (!is.null(input$nfactors)) {
      n_factors <- input$nfactors
    }
    
    # Créer un dataframe pour le graphique
    scree_data <- data.frame(
      Factor = 1:length(eigenvalues),
      Eigenvalue = eigenvalues
    )
    
    # Règle de Kaiser (valeurs propres > 1)
    kaiser_line <- data.frame(x = c(0, length(eigenvalues)), y = c(1, 1))
    
    # Créer le graphique
    ggplot(scree_data, aes(x = Factor, y = Eigenvalue)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_line(data = kaiser_line, aes(x = x, y = y), linetype = "dashed", color = "red") +
      geom_vline(xintercept = n_factors + 0.5, linetype = "dashed", color = "green") +
      scale_x_continuous(breaks = 1:length(eigenvalues)) +
      labs(title = "Graphique des valeurs propres (Scree Plot)",
           subtitle = "Ligne rouge = Règle de Kaiser (VP > 1) | Ligne verte = Facteurs retenus",
           x = "Facteur",
           y = "Valeur propre") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))
  })

  
  # Scree plot initial pour la matrice de corrélation chargée
  output$initial_scree_plot_cormat <- renderPlot({
    req(input$data_type == "cormat")
    req(cormat_data())
    
    # Récupérer la matrice de corrélation
    cor_mat <- cormat_data()
    
    # Calcul des valeurs propres
    eigenvalues <- eigen(cor_mat)$values
    
    # Déterminer le nombre de facteurs selon la règle de Kaiser ou la sélection de l'utilisateur
    n_factors <- sum(eigenvalues >= 1)  # Nombre selon le critère de Kaiser
    
    # Si aucun facteur n'est au-dessus de 1, prendre au moins 1 facteur
    if (n_factors == 0) n_factors <- 1
    
    # Utiliser le nombre de facteurs défini par l'utilisateur s'il est disponible
    if (!is.null(input$nfactors)) {
      n_factors <- input$nfactors
    }
    
    # Créer un dataframe pour le graphique
    scree_data <- data.frame(
      Factor = 1:length(eigenvalues),
      Eigenvalue = eigenvalues
    )
    
    # Règle de Kaiser (valeurs propres > 1)
    kaiser_line <- data.frame(x = c(0, length(eigenvalues)), y = c(1, 1))
    
    # Créer le graphique
    ggplot(scree_data, aes(x = Factor, y = Eigenvalue)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_line(data = kaiser_line, aes(x = x, y = y), linetype = "dashed", color = "red") +
      geom_vline(xintercept = n_factors + 0.5, linetype = "dashed", color = "green") +
      scale_x_continuous(breaks = 1:length(eigenvalues)) +
      labs(title = "Graphique des valeurs propres (Scree Plot)",
           subtitle = "Ligne rouge = Règle de Kaiser (VP > 1) | Ligne verte = Facteurs retenus",
           x = "Facteur",
           y = "Valeur propre") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))
  })
  
  # Tableau du critère de Kaiser pour les données brutes
  output$kaiser_criterion_table <- renderTable({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Calculer la matrice de corrélation
    if (input$cor == "cor") {
      cor_mat <- cor(analysis_data(), use = "pairwise.complete.obs", method = "pearson")
    } else if (input$cor == "spearman") {
      cor_mat <- cor(analysis_data(), use = "pairwise.complete.obs", method = "spearman")
    } else {
      cor_mat <- cov(analysis_data(), use = "pairwise.complete.obs")
    }
    
    # Calcul des valeurs propres
    eigenvalues <- eigen(cor_mat)$values
    
    # Déterminer combien de facteurs ont des valeurs propres > 1
    factors_to_retain <- sum(eigenvalues > 1)
    
    # Créer un dataframe pour le tableau
    kaiser_table <- data.frame(
      Facteur = 1:min(10, length(eigenvalues)),
      `Valeur propre` = eigenvalues[1:min(10, length(eigenvalues))],
      `À conserver` = ifelse(eigenvalues[1:min(10, length(eigenvalues))] > 1, "Oui", "Non")
    )
    
    # Ajouter une note de résumé au tableau
    summary_row <- data.frame(
      Facteur = "Résumé",
      `Valeur propre` = NA,
      `À conserver` = paste("Conserver", factors_to_retain, "facteurs")
    )
    
    rbind(kaiser_table, summary_row)
  }, rownames = FALSE, digits = 3)
  
  # Tableau du critère de Kaiser pour les données de matrice de corrélation
  output$kaiser_criterion_table_cormat <- renderTable({
    req(input$data_type == "cormat")
    req(cormat_data())
    
    # Récupérer la matrice de corrélation
    cor_mat <- cormat_data()
    
    # Calcul des valeurs propres
    eigenvalues <- eigen(cor_mat)$values
    
    # Déterminer combien de facteurs ont des valeurs propres > 1
    factors_to_retain <- sum(eigenvalues > 1)
    
    # Créer un dataframe pour le tableau
    kaiser_table <- data.frame(
      Facteur = 1:min(10, length(eigenvalues)),
      `Valeur propre` = eigenvalues[1:min(10, length(eigenvalues))],
      `À conserver` = ifelse(eigenvalues[1:min(10, length(eigenvalues))] > 1, "Oui", "Non")
    )
    
    # Ajouter une note de résumé au tableau
    summary_row <- data.frame(
      Facteur = "Résumé",
      `Valeur propre` = NA,
      `À conserver` = paste("Conserver", factors_to_retain, "facteurs")
    )
    
    rbind(kaiser_table, summary_row)
  }, rownames = FALSE, digits = 3)
  
  # Test MAP (Minimum Average Partial) de Velicer
  output$map_test_plot <- renderPlot({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Exécuter le test MAP de Velicer
    tryCatch({
      map_result <- psych::VSS(analysis_data(), n = min(ncol(analysis_data())-1, 10), plot = FALSE)
      
      # Vérifier si le test a réussi
      if (!is.null(map_result$map)) {
        # Créer un dataframe pour le graphique
        map_data <- data.frame(
          Factors = 0:(length(map_result$map)-1),
          MAP_Criterion = c(1, map_result$map)
        )
        
        # Déterminer le nombre optimal de facteurs (le point le plus bas)
        opt_factors <- which.min(map_data$MAP_Criterion) - 1
        
        # Créer le graphique
        ggplot(map_data, aes(x = Factors, y = MAP_Criterion)) +
          geom_line(linewidth = 1.2) +
          geom_point(size = 3) +
          geom_vline(xintercept = opt_factors, linetype = "dashed", color = "red") +
          scale_x_continuous(breaks = 0:max(map_data$Factors)) +
          labs(title = "Test MAP de Velicer",
               subtitle = paste("Le nombre optimal de facteurs est:", opt_factors),
               x = "Nombre de facteurs",
               y = "Critère MAP") +
          theme_minimal() +
          theme(plot.title = element_text(face = "bold"),
                axis.title = element_text(face = "bold"))
      } else {
        # Si le test n'a pas de résultats mais n'a pas généré d'erreur
        plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
        text(1, 1, "Impossible de calculer le test MAP.\nEssayez avec un jeu de données plus large ou\navec moins de corrélations nulles ou négatives.")
      }
    }, error = function(e) {
      # En cas d'erreur, afficher un message
      plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
      text(1, 1, "Impossible de calculer le test MAP.\nEssayez avec un jeu de données plus large ou\navec moins de corrélations nulles ou négatives.")
    })
  })
  
  # Test MAP (Minimum Average Partial) de Velicer
  map_result_reactive <- reactive({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Exécuter le test MAP de Velicer
    tryCatch({
      result <- psych::VSS(analysis_data(), n = min(ncol(analysis_data())-1, 10), plot = FALSE)
      return(result)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  output$map_test_plot <- renderPlot({
    req(map_result_reactive())
    result <- map_result_reactive()
    
    # Vérifier si le test a réussi
    if (!is.null(result$map)) {
      # Créer un dataframe pour le graphique - CORRECTION : s'assurer que les deux vecteurs ont la même longueur
      map_values <- result$map
      map_data <- data.frame(
        Factors = 0:(length(map_values)),  # Inclure le facteur 0 pour l'indice de départ
        MAP_Criterion = c(1, map_values)   # Ajouter 1 comme première valeur
      )
      
      # Déterminer le nombre optimal de facteurs (le point le plus bas)
      opt_factors <- which.min(map_data$MAP_Criterion) - 1
      
      # Créer le graphique
      ggplot(map_data, aes(x = Factors, y = MAP_Criterion)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        geom_vline(xintercept = opt_factors, linetype = "dashed", color = "red") +
        scale_x_continuous(breaks = 0:max(map_data$Factors)) +
        labs(title = "Test MAP de Velicer",
             subtitle = paste("Le nombre optimal de facteurs est:", opt_factors),
             x = "Nombre de facteurs",
             y = "Critère MAP") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              axis.title = element_text(face = "bold"))
    } else {
      # Afficher un message d'erreur
      plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
      text(1, 1, "Impossible de calculer le test MAP.\nEssayez avec un jeu de données plus large ou\navec moins de corrélations nulles ou négatives.")
    }
  })
  
  # Résultats du test MAP
  output$map_test_result <- renderPrint({
    req(map_result_reactive())
    result <- map_result_reactive()
    
    # Vérifier si le test a réussi
    if (!is.null(result$map)) {
      # Ajouter 1 comme première valeur pour le critère MAP complet
      map_values <- c(1, result$map)
      
      # Déterminer le nombre optimal de facteurs (le point le plus bas)
      opt_factors <- which.min(map_values) - 1
      
      cat("Test MAP de Velicer\n")
      cat("------------------\n")
      cat("Le nombre optimal de facteurs selon le critère MAP est:", opt_factors, "\n")
      cat("\n")
      cat("Valeurs du critère MAP:\n")
      print(round(map_values, 4))
    } else {
      cat("Impossible de calculer le test MAP.\n")
      cat("Essayez avec un jeu de données plus large ou avec moins de corrélations nulles ou négatives.")
    }
  })
  # Définir un réactif pour savoir si nous devons forcer le recalcul
  force_parallel_calc <- reactiveVal(FALSE)
  
  # Observer pour l'onglet de test parallèle
  observeEvent(input$factor_choice_tabs, {
    if (input$factor_choice_tabs == "Test parallèle" && input$data_type == "raw" && !is.null(analysis_data())) {
      # Forcer le calcul du test parallèle quand on arrive sur cet onglet
      force_parallel_calc(TRUE)
    }
  })
  
  # Observer pour le bouton d'actualisation du test parallèle
  observeEvent(input$run_parallel_auto, {
    if (input$data_type == "raw" && !is.null(analysis_data())) {
      # Forcer la réexécution du test parallèle
      force_parallel_calc(TRUE)
    }
  })
  
  # Fonction réactive pour le test parallèle automatique
  parallel_results_auto <- eventReactive(
    c(input$n_datasets_auto, input$percentile_auto, input$run_parallel_auto, force_parallel_calc()), 
    {
      # Ne pas exécuter si ce n'est pas l'onglet actif, sauf si forcé
      if (input$factor_choice_tabs != "Test parallèle" && !force_parallel_calc()) {
        return(NULL)
      }
      
      # Réinitialiser le drapeau de forçage
      force_parallel_calc(FALSE)
      
      req(input$data_type == "raw")
      req(analysis_data())
      
      withProgress(message = 'Exécution du test parallèle...', {
        tryCatch({
          # Exécuter l'analyse parallèle
          result <- fa.parallel(analysis_data(), 
                                fm = input$fm, 
                                n.iter = input$n_datasets_auto,
                                cor = input$cor,
                                use = "pairwise.complete.obs",
                                plot = FALSE)  # Ne pas afficher le graphique ici
          
          return(result)
        }, error = function(e) {
          showNotification(paste("Erreur dans le test parallèle:", e$message), type = "error")
          return(NULL)
        })
      })
    }
  )
  
  # Graphique du test parallèle intégré à l'interface
  output$parallel_plot_auto <- renderPlot({
    req(parallel_results_auto())
    
    result <- parallel_results_auto()
    
    # Extraire les valeurs propres observées
    observed <- result$fa.values
    
    # Calculer le percentile spécifié des valeurs simulées
    percentile_value <- input$percentile_auto / 100
    simulated_percentiles <- apply(result$values, 2, function(x) quantile(x, probs = percentile_value))
    
    # S'assurer que nous avons assez de valeurs simulées
    n_factors <- length(observed)
    if (length(simulated_percentiles) < n_factors) {
      # Compléter avec la dernière valeur si nécessaire
      simulated_percentiles <- c(simulated_percentiles, 
                                 rep(simulated_percentiles[length(simulated_percentiles)], 
                                     n_factors - length(simulated_percentiles)))
    } else {
      # Limiter au nombre de facteurs observés
      simulated_percentiles <- simulated_percentiles[1:n_factors]
    }
    
    # Créer un dataframe pour le graphique
    df <- data.frame(
      Factor = c(1:n_factors, 1:n_factors),
      Eigenvalue = c(observed, simulated_percentiles),
      Type = c(rep("Données observées", n_factors), 
               rep(paste0("Données simulées (P", input$percentile_auto, ")"), n_factors))
    )
    
    # Déterminer le nombre de facteurs à retenir (où valeurs observées > percentile simulé)
    greater_than_random <- observed > simulated_percentiles
    factors_to_retain <- sum(greater_than_random)
    
    # Nous devons créer le mapping des couleurs avant de l'utiliser
    color_values <- c("blue", "red")
    names(color_values) <- c("Données observées", paste0("Données simulées (P", input$percentile_auto, ")"))
    
    # Créer le graphique avec ggplot2
    p <- ggplot(df, aes(x = Factor, y = Eigenvalue, group = Type, color = Type)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_color_manual(values = color_values) +
      labs(title = "Analyse parallèle",
           subtitle = paste("Comparaison des valeurs propres observées avec le", 
                            input$percentile_auto, "e percentile de", 
                            input$n_datasets_auto, "jeux de données simulés"),
           x = "Facteur", y = "Valeur propre") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_blank())
    
    # Ajouter la ligne verticale et l'annotation pour le nombre de facteurs à retenir
    if (factors_to_retain > 0) {
      p <- p + 
        geom_vline(xintercept = factors_to_retain + 0.5, linetype = "dashed", color = "darkgreen") +
        annotate("text", x = factors_to_retain + 1, y = max(df$Eigenvalue) * 0.9, 
                 label = paste("Retenir", factors_to_retain, "facteurs"), 
                 color = "darkgreen", fontface = "bold")
    }
    
    return(p)
  })
  
  # Résultats textuels du test parallèle
  output$parallel_summary_auto <- renderPrint({
    req(parallel_results_auto())
    
    result <- parallel_results_auto()
    
    # Extraire les valeurs propres observées
    observed <- result$fa.values
    
    # Calculer le percentile spécifié des valeurs simulées
    percentile_value <- input$percentile_auto / 100
    simulated_percentiles <- apply(result$values, 2, function(x) quantile(x, probs = percentile_value))
    
    # S'assurer que nous avons assez de valeurs simulées
    n_factors <- length(observed)
    if (length(simulated_percentiles) < n_factors) {
      # Compléter avec la dernière valeur si nécessaire
      simulated_percentiles <- c(simulated_percentiles, 
                                 rep(simulated_percentiles[length(simulated_percentiles)], 
                                     n_factors - length(simulated_percentiles)))
    } else {
      # Limiter au nombre de facteurs observés
      simulated_percentiles <- simulated_percentiles[1:n_factors]
    }
    
    # Déterminer le nombre de facteurs à retenir (où valeurs observées > percentile simulé)
    greater_than_random <- observed > simulated_percentiles
    factors_to_retain <- sum(greater_than_random)
    
    cat("Résultats de l'analyse parallèle\n")
    cat("------------------------------\n")
    cat("Nombre de facteurs à retenir selon l'analyse parallèle:", factors_to_retain, "\n")
    cat("(facteurs dont la valeur propre observée > percentile", input$percentile_auto, "des valeurs simulées)\n\n")
    
    # Notons aussi la valeur renvoyée par fa.parallel() qui peut être différente
    if (result$nfact != factors_to_retain) {
      cat("Note: La fonction fa.parallel() suggère", result$nfact, "facteurs, mais l'examen des valeurs\n")
      cat("montre que", factors_to_retain, "facteurs ont des valeurs propres supérieures au percentile", 
          input$percentile_auto, "des données simulées.\n\n")
    }
    
    cat("Comparaison des valeurs propres:\n")
    
    # S'assurer que les dimensions correspondent
    n_to_compare <- min(length(observed), length(simulated_percentiles), 10)
    
    # Créer d'abord les noms de colonnes
    col_names <- c("Facteur", "Valeur propre observée", 
                   paste0("Valeur simulée (P", input$percentile_auto, ")"), 
                   "Différence", "À retenir")
    
    # Créer le dataframe avec des noms de colonnes standard
    comparison <- data.frame(
      "Facteur" = 1:n_to_compare,
      "Valeur_propre_observee" = observed[1:n_to_compare],
      "Valeur_simulee" = simulated_percentiles[1:n_to_compare],
      "Difference" = observed[1:n_to_compare] - simulated_percentiles[1:n_to_compare],
      "A_retenir" = ifelse(observed[1:n_to_compare] > simulated_percentiles[1:n_to_compare], "Oui", "Non")
    )
    
    # Renommer les colonnes après la création
    names(comparison) <- col_names
    
    print(comparison)
  })

  # Réactif pour l'analyse EGA - version améliorée pour gérer GLASSO sur matrices de corrélation
  ega_results <- reactive({
    req(correlation_matrix_for_analysis())
    
    # Surveiller les changements dans les paramètres EGA selon le type de données
    if (input$data_type == "raw") {
      model <- input$ega_model
      algorithm <- input$ega_algorithm
    } else {
      model <- input$ega_model_cormat
      algorithm <- input$ega_algorithm_cormat
    }
    
    withProgress(message = 'Exécution de l\'analyse EGA...', {
      tryCatch({
        # Récupérer la matrice de corrélation
        cor_mat <- correlation_matrix_for_analysis()
          ega_result <- EGA(
            data = cor_mat,
            n = n_observations(),
            model = model,
            algorithm = algorithm,
            plot.EGA = FALSE
          )
          
          return(ega_result)
      }, error = function(e) {
        showNotification(paste("Erreur dans l'analyse EGA:", e$message), type = "error")
        return(NULL)
      })
    })
  })
  
  # Graphique EGA avec options contrôlables - unified for both data types
  output$ega_plot <- renderPlot({
    req(ega_results())
    
    # Déterminer les paramètres d'affichage en fonction du type de données
    if (input$data_type == "raw") {
      node_size <- input$ega_node_size
      label_size <- input$ega_label_size
      edge_alpha <- input$ega_edge_alpha
      layout <- input$ega_layout
    } else {
      node_size <- input$ega_node_size_cormat
      label_size <- input$ega_label_size_cormat
      edge_alpha <- input$ega_edge_alpha_cormat
      layout <- input$ega_layout_cormat
    }
    
    # Afficher le graphique EGA avec les options choisies
    plot(ega_results(), 
         title = "Analyse Graphique Exploratoire (EGA)",
         vsize = node_size,
         label.size = label_size,
         edge.alpha = edge_alpha,
         layout = layout
    )
  })
  
  # Rediriger output$ega_plot_cormat vers output$ega_plot
  output$ega_plot_cormat <- renderPlot({
    req(ega_results())
    output$ega_plot()
  })
  
  # Résumé EGA
  output$ega_summary <- renderUI({
    req(ega_results())
    
    result <- ega_results()
    
    # Calculer les statistiques pour chaque dimension
    dim_table <- table(result$wc)
    
    # Créer un dataframe pour les variables par dimension
    var_dims <- data.frame(
      Variable = names(result$wc),
      Dimension = result$wc
    )
    
    # Trier par dimension, puis par nom de variable
    var_dims <- var_dims[order(var_dims$Dimension, var_dims$Variable), ]
    
    # Créer une liste de variables par dimension pour l'affichage en accordion
    dim_vars_list <- split(var_dims$Variable, var_dims$Dimension)
    
    # Couleurs pour les dimensions (palette RColorBrewer Set1)
    dim_colors <- RColorBrewer::brewer.pal(min(9, result$n.dim), "Set1")
    if (result$n.dim > 9) {
      dim_colors <- colorRampPalette(dim_colors)(result$n.dim)
    }
    
    # Créer l'interface utilisateur
    tagList(
      div(
        class = "mb-4",
        h3("Résultats de l'Analyse Graphique Exploratoire (EGA)"),
        
        # Informations principales dans une card
        card(
          card_header("Informations générales"),
          div(
            class = "p-3",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h4("Nombre de dimensions détectées"),
                div(
                  class = "display-4 text-center mb-3",
                  result$n.dim
                ),
                
                hr(),
                
                # Ajouter l'indice TEFI
                h5("Indice d'ajustement TEFI"),
                div(
                  class = "text-center",
                  div(
                    style = paste0(
                      "display: inline-block; padding: 8px 16px; border-radius: 5px; ",
                      "background-color: ", 
                      ifelse(result$TEFI < 0.3, "#d4edda", 
                             ifelse(result$TEFI < 0.5, "#fff3cd", "#f8d7da")),
                      "; font-weight: bold; font-size: 18px;"
                    ),
                    round(result$TEFI, 4)
                  ),
                  div(
                    class = "small text-muted mt-2",
                    "Indice de qualité d'ajustement (valeurs plus basses indiquent un meilleur ajustement)"
                  )
                )
              ),
              
              # Graphique de répartition des variables par dimension
              div(
                h4("Répartition des variables par dimension"),
                plotOutput("ega_dimension_distribution", height = "200px")
              )
            )
          )
        ),
        
        # Tableau détaillé des variables par dimension
        card(
          card_header("Distribution des variables par dimension"),
          div(
            class = "p-3",
            tabsetPanel(
              # Onglet 1: Tableau des variables
              tabPanel("Tableau des variables",
                       div(
                         class = "mt-3",
                         DTOutput("ega_dimension_table")
                       )
              ),
              
              # Onglet 2: Visualisation par dimension (accordion)
              tabPanel("Visualisation par dimension",
                       div(
                         class = "mt-3",
                         accordion(
                           id = "dim_accordion",
                           multiple = TRUE,
                           # Créer un panneau d'accordéon pour chaque dimension
                           lapply(1:result$n.dim, function(dim) {
                             accordion_panel(
                               title = paste0("Dimension ", dim, " (", length(dim_vars_list[[as.character(dim)]]), " variables)"),
                               div(
                                 style = paste0("border-left: 5px solid ", dim_colors[dim], "; padding-left: 15px;"),
                                 # Afficher les variables de cette dimension avec des puces
                                 lapply(dim_vars_list[[as.character(dim)]], function(var) {
                                   div(
                                     class = "mb-1",
                                     icon("circle"), " ", var
                                   )
                                 })
                               )
                             )
                           })
                         )
                       )
              ),
              
              # Onglet 3: Statistiques des dimensions
              tabPanel("Statistiques des dimensions",
                       div(
                         class = "mt-3",
                         h5("Nombre de variables par dimension"),
                         layout_columns(
                           col_widths = c(5, 7),
                           datatable(
                             data.frame(
                               Dimension = 1:length(dim_table),
                               `Nombre de variables` = as.numeric(dim_table),
                               `Pourcentage` = paste0(round(100 * as.numeric(dim_table) / sum(dim_table), 1), "%")
                             ),
                             options = list(dom = 't', ordering = FALSE),
                             rownames = FALSE
                           ),
                           plotOutput("ega_dim_barplot", height = "200px")
                         )
                       )
              )
            )
          )
        ),
        
        # Section supplémentaire pour l'interprétation
        card(
          card_header("Interprétation des résultats"),
          div(
            class = "p-3",
            p("L'analyse EGA a identifié ", tags$strong(result$n.dim), " dimensions distinctes dans les données:"),
            
            div(
              class = "mb-3 p-3 bg-light rounded",
              h5("Signification des dimensions"),
              p("Les dimensions identifiées par EGA représentent des clusters cohérents de variables, reflétant potentiellement des construits latents dans les données. Elles peuvent être comparées aux facteurs identifiés par l'analyse factorielle traditionnelle."),
              
              h5("Comment interpréter l'indice TEFI"),
              p("L'indice TEFI (Total Entropy Fit Index) mesure la qualité d'ajustement du modèle EGA:"),
              tags$ul(
                tags$li("Valeurs plus basses indiquent un meilleur ajustement"),
                tags$li("Généralement interprété en comparant différents modèles"),
                tags$li("Pas de seuil universel, mais les modèles avec TEFI < 0.3 sont souvent considérés comme bien ajustés")
              )
            )
          )
        )
      )
    )
  })
  
  # Ajouter ces nouveaux outputs pour les visualisations
  output$ega_dimension_distribution <- renderPlot({
    req(ega_results())
    
    # Créer un dataframe pour le graphique
    dim_counts <- as.data.frame(table(ega_results()$wc))
    names(dim_counts) <- c("Dimension", "Count")
    dim_counts$Dimension <- as.numeric(as.character(dim_counts$Dimension))
    
    # Couleurs pour les dimensions (palette RColorBrewer Set1)
    dim_colors <- RColorBrewer::brewer.pal(min(9, ega_results()$n.dim), "Set1")
    if (ega_results()$n.dim > 9) {
      dim_colors <- colorRampPalette(dim_colors)(ega_results()$n.dim)
    }
    
    # Créer un graphique en secteurs
    ggplot(dim_counts, aes(x = "", y = Count, fill = factor(Dimension))) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = dim_colors) +
      labs(fill = "Dimension") +
      theme_void() +
      theme(legend.position = "right")
  })
  
  output$ega_dimension_table <- renderDT({
    req(ega_results())
    
    # Créer un dataframe pour le tableau
    var_dims <- data.frame(
      Variable = names(ega_results()$wc),
      Dimension = ega_results()$wc
    )
    
    # Trier par dimension, puis par nom de variable
    var_dims <- var_dims[order(var_dims$Dimension, var_dims$Variable), ]
    
    # Créer un tableau interactif
    datatable(
      var_dims,
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50),
        order = list(list(1, 'asc')),
        dom = 'lfrtip'
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'Dimension',
        backgroundColor = styleEqual(
          1:ega_results()$n.dim,
          RColorBrewer::brewer.pal(min(9, ega_results()$n.dim), "Pastel1") %>%
            {if(ega_results()$n.dim > 9) colorRampPalette(.)(ega_results()$n.dim) else .}
        )
      )
  })
  
  output$ega_dim_barplot <- renderPlot({
    req(ega_results())
    
    # Créer un dataframe pour le graphique
    dim_counts <- as.data.frame(table(ega_results()$wc))
    names(dim_counts) <- c("Dimension", "Count")
    dim_counts$Dimension <- as.numeric(as.character(dim_counts$Dimension))
    
    # Couleurs pour les dimensions (palette RColorBrewer Set1)
    dim_colors <- RColorBrewer::brewer.pal(min(9, ega_results()$n.dim), "Set1")
    if (ega_results()$n.dim > 9) {
      dim_colors <- colorRampPalette(dim_colors)(ega_results()$n.dim)
    }
    
    # Créer un graphique en barres
    ggplot(dim_counts, aes(x = factor(Dimension), y = Count, fill = factor(Dimension))) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = dim_colors) +
      labs(x = "Dimension", y = "Nombre de variables") +
      theme_minimal() +
      theme(legend.position = "none",
            panel.grid.minor = element_blank(),
            axis.title = element_text(face = "bold"))
  })
  
  output$ega_plot <- renderPlot({
    req(ega_results())
    
    # Afficher le graphique EGA avec les options choisies par l'utilisateur
    plot(ega_results(), 
         title = "Analyse Graphique Exploratoire (EGA)",
         vsize = input$ega_node_size,   # taille des nœuds contrôlable
         label.size = input$ega_label_size,  # taille des étiquettes contrôlable
         edge.alpha = input$ega_edge_alpha,  # transparence des arêtes contrôlable 
         layout = input$ega_layout  # layout contrôlable
    )
  })
  
  # Analyse factorielle ou ACP
  fa_results <- eventReactive(input$run_analysis, {
    req(analysis_data())
    
    # Paramètres communs pour les deux types d'analyse
    params <- list(
      nfactors = input$nfactors,
      rotate = input$rotate
    )
    
    # Paramètres spécifiques selon le type d'analyse
    if (input$analysis_type == "fa") {
      # Paramètres spécifiques à l'analyse factorielle
      params$fm <- input$fm
      params$SMC <- input$use_smc
      params$max.iter <- input$max_iter  # Seulement pour FA
    }
    
    # Paramètres spécifiques selon le type de données
    if (input$data_type == "raw") {
      params$r <- analysis_data()
      params$cor <- input$cor
      params$use <- "pairwise.complete.obs"
    } else {
      params$r <- correlation_matrix()
      params$n.obs <- input$n_obs
    }
    
    # Ajouter les paramètres spécifiques à la rotation Promax
    if (input$rotate == "promax" && !is.null(input$promax_k)) {
      params$promax.k <- input$promax_k
    }
    
    # Paramètre pour la normalisation de Kaiser
    if (input$kaiser && input$rotate %in% c("varimax", "quartimax", "equamax", "promax", "oblimin")) {
      params$normalize <- TRUE
    } else {
      params$normalize <- FALSE
    }
    
    withProgress(message = 'Exécution de l\'analyse...', {
      tryCatch({
        # Choisir la fonction à utiliser selon le type d'analyse
        func_to_call <- if(input$analysis_type == "fa") psych::fa else psych::principal
        do.call(func_to_call, params)
      }, error = function(e) {
        showNotification(paste("Erreur dans l'analyse:", e$message), type = "error")
        return(NULL)
      })
    })
  })
  
  output$analysis_done <- reactive({
    return(!is.null(fa_results()))
  })
  
  # Output UI dynamique pour les résultats
  output$results_content <- renderUI({
    req(fa_results())
    
    # Créer le contenu de base du tabset panel
    main_tabs <- tabsetPanel(
      id = "results_tabset",
      tabPanel("Indices d'ajustement", 
               uiOutput("fit_indices")),
      tabPanel("Variance expliquée", 
               layout_columns(
                 col_widths = c(6, 6),
                 card(
                   card_header("Graphique des valeurs propres"),
                   plotOutput("scree_plot", height = "300px")
                 ),
                 card(
                   card_header("Détails sur la variance expliquée"),
                   verbatimTextOutput("var_explained"),
                   div(
                     style = "padding: 10px; font-style: italic; font-size: 90%; color: #666;",
                     p("Interprétation de la variance expliquée:"),
                     tags$ul(
                       tags$li("La colonne SS loadings indique la somme des carrés des charges factorielles pour chaque facteur"),
                       tags$li("Proportion Var représente le pourcentage de variance expliquée par chaque facteur"),
                       tags$li("Cumulative Var indique la variance expliquée cumulée")
                     ),
                     p("NOTE: Il n'existe pas de seuil universellement accepté, car cela dépend du domaine d'application, du type de données, et du nombre de variables. En sciences sociales, on considère généralement qu'une solution factorielle expliquant 60% ou plus de la variance totale est satisfaisante."),
                     tags$ul(
                       tags$li("Plus le nombre de variables est élevé, plus il est difficile d'atteindre une grande proportion de variance expliquée."),
                       tags$li("Une bonne interprétabilité des facteurs peut compenser une variance expliquée un peu plus faible."),
                       tags$li("Si la solution factorielle est stable, interprétable et confirmée par d'autres indices (MAP, VSS, RMSEA, etc.), alors la variance expliquée n'a pas besoin d'être élevée.")
                     )
                    )
                 )
               )
              ),
      tabPanel("Charges factorielles", 
               layout_sidebar(
                 sidebar = sidebar(
                   width = 300,
                   card(
                     card_header("Options d'affichage"),
                     h5("Options des charges factorielles"),
                     numericInput("loading_cutoff", "Seuil des charges factorielles :", 0.3, min = 0, max = 1, step = 0.05),
                     selectInput("sort_loadings", "Trier les charges :", 
                                 choices = c("Par taille" = "size", 
                                             "Par facteur" = "factor", 
                                             "Par variable" = "variable"),
                                 selected = "factor"),
                     checkboxInput("show_complex", "Afficher les variables complexes", TRUE),
                     checkboxInput("show_heywood", "Mettre en évidence cas Heywood", TRUE),
                     checkboxInput("highlight_crossloadings", "Identifier charges croisées", TRUE),
                     numericInput("crossloading_diff", "Seuil différence charges croisées:", 0.2, min = 0.1, max = 0.5, step = 0.05),
                     checkboxInput("highlight_factors", "Mettre en couleur par facteur", TRUE),
                     
                     br()
                   )
                 ),
                 layout_columns(
                   col_widths = c(6, 6),
                   card(
                     card_header("Tableau des charges"),
                     DTOutput("loadings_table")
                   ),
                   card(
                     card_header("Représentation graphique"),
                     plotOutput("loadings_plot", height = "500px"),
                     div(
                       style = "text-align: right; margin-top: 10px;",
                       actionButton("show_lavaan_syntax", "Voir syntaxe lavaan", 
                                    icon = icon("code"), class = "btn-sm btn-info")
                     )
                   )
                 )
               )),
      tabPanel("Structure factorielle", 
               card(
                 card_header("Structure factorielle"),
                 uiOutput("factor_structure")
               )),
      tabPanel("Corrélations entre facteurs", 
               card(
                 card_header("Corrélations entre facteurs"),
                 conditionalPanel(
                   condition = "input.rotate != 'none' && input.rotate != 'varimax' && input.rotate != 'quartimax' && input.rotate != 'equamax'",
                   plotOutput("factor_cors_plot", height = "400px"),
                   DTOutput("factor_cors_table")
                 ),
                 conditionalPanel(
                   condition = "input.rotate == 'none' || input.rotate == 'varimax' || input.rotate == 'quartimax' || input.rotate == 'equamax'",
                   p("Les rotations orthogonales produisent des facteurs non corrélés.")
                 )
               )),
      tabPanel("Résidus", 
               card(
                 card_header("Résidus"),
                 plotOutput("residuals_plot", height = "500px"),
                 DTOutput("residuals_table")
               )),
      tabPanel("Communalités", 
               card(
                 card_header("Communalités"),
                 layout_columns(
                   col_widths = c(6, 6),
                   card(
                     card_header("Graphique des communalités"),
                     plotOutput("communalities_plot", height = "400px")
                   ),
                   card(
                     card_header("Tableau des communalités"),
                     DTOutput("communalities_table")
                   )
                 )
               )),
      tabPanel("Résumé complet", 
               card(
                 card_header("Résumé complet"),
                 verbatimTextOutput("full_summary")
               ))
    )
    
    # Emballer dans une carte
    card(
      card_header("Résultats de l'analyse factorielle"),
      main_tabs
    )
  })
  
  
  outputOptions(output, "analysis_done", suspendWhenHidden = FALSE)

  output$fit_indices <- renderUI({
    req(fa_results())
    
    if (input$analysis_type == "fa") {
      # Pour l'analyse factorielle
      
      # Définir les seuils pour chaque indice
      thresholds <- list(
        pval = c(0.01, 0.05),
        rms = c(0.10, 0.08),
        GFI = c(0.90, 0.95),
        AGFI = c(0.90, 0.95),
        rmsr_ratio = c(1.0, 0.8),
        TLI = c(0.90, 0.95),
        CFI = c(0.90, 0.95),
        RMSEA = c(0.08, 0.05)
      )
      
      # Créer une fonction pour déterminer le statut d'un indice (et la couleur associée)
      get_status <- function(value, thresholds, is_decreasing = FALSE) {
        if (is.na(value)) return(list(status = "Non disponible", color = "#e9ecef"))
        
        if (is_decreasing) {
          # Pour indices où valeurs plus basses sont meilleures
          if (value <= thresholds[2]) {
            return(list(status = "Bon", color = "#d4edda"))
          } else if (value <= thresholds[1]) {
            return(list(status = "Acceptable", color = "#fff3cd"))
          } else {
            return(list(status = "Médiocre", color = "#f8d7da"))
          }
        } else {
          # Pour indices où valeurs plus élevées sont meilleures
          if (value >= thresholds[2]) {
            return(list(status = "Bon", color = "#d4edda"))
          } else if (value >= thresholds[1]) {
            return(list(status = "Acceptable", color = "#fff3cd"))
          } else {
            return(list(status = "Médiocre", color = "#f8d7da"))
          }
        }
      }
      
      # Extraire les valeurs des indices d'ajustement
      indices <- list(
        pval = fa_results()$PVAL,
        rms = fa_results()$rms,
        GFI = if(is.null(fa_results()$GFI)) NA else fa_results()$GFI,
        AGFI = if(is.null(fa_results()$AGFI)) NA else fa_results()$AGFI,
        rmsr_ratio = fa_results()$rms/fa_results()$crms,
        TLI = fa_results()$TLI,
        CFI = if(is.null(fa_results()$CFI)) NA else fa_results()$CFI,
        RMSEA = fa_results()$RMSEA[1]
      )
      
      # Interprétation globale
      if (indices$pval > 0.05 && indices$TLI > 0.95 && indices$RMSEA < 0.05) {
        global_interpretation <- "Le modèle présente un bon ajustement selon plusieurs indices."
        global_color <- "#d4edda"
      } else if (indices$pval > 0.01 && indices$TLI > 0.90 && indices$RMSEA < 0.08) {
        global_interpretation <- "Le modèle présente un ajustement acceptable."
        global_color <- "#fff3cd"
      } else {
        global_interpretation <- "Certains indices suggèrent un ajustement insuffisant du modèle. Considérez de modifier le nombre de facteurs ou la méthode d'extraction."
        global_color <- "#f8d7da"
      }
      
      # Créer la table des indices d'ajustement
      div(
        h3("Indices d'ajustement du modèle", class = "mt-3 mb-4"),
        
        # Interprétation globale en haut
        div(class = "alert", style = paste0("background-color:", global_color, "; border-color:", global_color),
            h4("Interprétation globale:", class = "alert-heading"),
            p(global_interpretation)
        ),
        
        # Indices d'ajustement absolus
        card(
          card_header(h4("Indices d'ajustement absolus", class = "m-0")),
          div(class = "table-responsive",
              tags$table(class = "table table-hover",
                         tags$thead(
                           tags$tr(
                             tags$th("Indice"),
                             tags$th("Valeur"),
                             tags$th("Interprétation"),
                             tags$th("Statut")
                           )
                         ),
                         tags$tbody(
                           # Chi-carré
                           tags$tr(
                             tags$td(
                               "Chi-carré",
                               tags$div(class = "small text-muted", "Test d'adéquation du modèle")
                             ),
                             tags$td(paste0(
                               round(fa_results()$STATISTIC, 2), 
                               " (ddl = ", fa_results()$dof, 
                               ", p = ", format.pval(indices$pval, digits = 3), ")"
                             )),
                             tags$td("Une p-value > 0.05 suggère un bon ajustement. Cependant, ce test est sensible à la taille de l'échantillon."),
                             tags$td(style = paste0("background-color:", get_status(indices$pval, thresholds$pval, FALSE)$color),
                                     get_status(indices$pval, thresholds$pval, FALSE)$status)
                           ),
                           
                           # SRMR
                           tags$tr(
                             tags$td(
                               "SRMR",
                               tags$div(class = "small text-muted", "Standardized Root Mean Square Residual")
                             ),
                             tags$td(round(indices$rms, 3)),
                             tags$td("Racine carrée de la moyenne des résidus standardisés au carré. Valeurs < 0.08 suggèrent un bon ajustement."),
                             tags$td(style = paste0("background-color:", get_status(indices$rms, thresholds$rms, TRUE)$color),
                                     get_status(indices$rms, thresholds$rms, TRUE)$status)
                           ),
                           
                           # GFI
                           if (!is.null(fa_results()$GFI)) {
                             tags$tr(
                               tags$td(
                                 "GFI",
                                 tags$div(class = "small text-muted", "Goodness of Fit Index")
                               ),
                               tags$td(round(indices$GFI, 3)),
                               tags$td("Mesure la variance/covariance relative expliquée par le modèle. Valeurs > 0.95 indiquent un bon ajustement."),
                               tags$td(style = paste0("background-color:", get_status(indices$GFI, thresholds$GFI, FALSE)$color),
                                       get_status(indices$GFI, thresholds$GFI, FALSE)$status)
                             )
                           },
                           
                           # AGFI
                           if (!is.null(fa_results()$AGFI)) {
                             tags$tr(
                               tags$td(
                                 "AGFI",
                                 tags$div(class = "small text-muted", "Adjusted Goodness of Fit Index")
                               ),
                               tags$td(round(indices$AGFI, 3)),
                               tags$td("GFI ajusté pour les degrés de liberté du modèle. Valeurs > 0.95 indiquent un bon ajustement."),
                               tags$td(style = paste0("background-color:", get_status(indices$AGFI, thresholds$AGFI, FALSE)$color),
                                       get_status(indices$AGFI, thresholds$AGFI, FALSE)$status)
                             )
                           },
                           
                           # RMSR / RMSE ratio
                           tags$tr(
                             tags$td(
                               "RMSR / RMSE ratio",
                               tags$div(class = "small text-muted", "Rapport entre RMSR observé et RMSE attendu")
                             ),
                             tags$td(round(indices$rmsr_ratio, 3)),
                             tags$td("Rapport entre RMSR observé et RMSE attendu. Valeurs < 1.0 indiquent un bon ajustement, < 0.8 excellent."),
                             tags$td(style = paste0("background-color:", get_status(indices$rmsr_ratio, thresholds$rmsr_ratio, TRUE)$color),
                                     get_status(indices$rmsr_ratio, thresholds$rmsr_ratio, TRUE)$status)
                           )
                         )
              )
          )
        ),
        
        # Indices d'ajustement incrémentaux
        card(
          class = "mt-4",
          card_header(h4("Indices d'ajustement incrémentaux", class = "m-0")),
          div(class = "table-responsive",
              tags$table(class = "table table-hover",
                         tags$thead(
                           tags$tr(
                             tags$th("Indice"),
                             tags$th("Valeur"),
                             tags$th("Interprétation"),
                             tags$th("Statut")
                           )
                         ),
                         tags$tbody(
                           # TLI
                           tags$tr(
                             tags$td(
                               "TLI",
                               tags$div(class = "small text-muted", "Tucker Lewis Index")
                             ),
                             tags$td(round(indices$TLI, 3)),
                             tags$td("Compare l'ajustement du modèle à celui d'un modèle nul. Les valeurs > 0.95 indiquent un bon ajustement."),
                             tags$td(style = paste0("background-color:", get_status(indices$TLI, thresholds$TLI, FALSE)$color),
                                     get_status(indices$TLI, thresholds$TLI, FALSE)$status)
                           ),
                           
                           # CFI
                           if (!is.null(fa_results()$CFI)) {
                             tags$tr(
                               tags$td(
                                 "CFI",
                                 tags$div(class = "small text-muted", "Comparative Fit Index")
                               ),
                               tags$td(round(indices$CFI, 3)),
                               tags$td("Compare l'ajustement du modèle à celui d'un modèle d'indépendance. Les valeurs > 0.95 indiquent un bon ajustement."),
                               tags$td(style = paste0("background-color:", get_status(indices$CFI, thresholds$CFI, FALSE)$color),
                                       get_status(indices$CFI, thresholds$CFI, FALSE)$status)
                             )
                           }
                         )
              )
          )
        ),
        
        # Indices d'ajustement de parcimonie
        card(
          class = "mt-4",
          card_header(h4("Indices d'ajustement de parcimonie", class = "m-0")),
          div(class = "table-responsive",
              tags$table(class = "table table-hover",
                         tags$thead(
                           tags$tr(
                             tags$th("Indice"),
                             tags$th("Valeur"),
                             tags$th("Interprétation"),
                             tags$th("Statut")
                           )
                         ),
                         tags$tbody(
                           # RMSEA
                           tags$tr(
                             tags$td(
                               "RMSEA",
                               tags$div(class = "small text-muted", "Root Mean Square Error of Approximation")
                             ),
                             tags$td(
                               HTML(paste0(
                                 round(indices$RMSEA, 3),
                                 "<br><span class='small text-muted'>IC 90% [", 
                                 round(fa_results()$RMSEA[2], 3), "-", 
                                 round(fa_results()$RMSEA[3], 3), "]</span>"
                               ))
                             ),
                             tags$td("Mesure l'erreur d'approximation du modèle par degré de liberté. Valeurs < 0.05 indiquent un bon ajustement."),
                             tags$td(style = paste0("background-color:", get_status(indices$RMSEA, thresholds$RMSEA, TRUE)$color),
                                     get_status(indices$RMSEA, thresholds$RMSEA, TRUE)$status)
                           ),
                           
                           # BIC
                           tags$tr(
                             tags$td(
                               "BIC",
                               tags$div(class = "small text-muted", "Bayesian Information Criterion")
                             ),
                             tags$td(round(fa_results()$BIC, 2)),
                             tags$td("Critère d'information bayésien prenant en compte la complexité du modèle. Utile pour comparer différents modèles (valeur plus basse = meilleur modèle)."),
                             tags$td("Comparatif")
                           )
                         )
              )
          )
        ),
        
        # Informations supplémentaires et références
        card(
          class = "mt-4",
          card_header(h4("Pour aller plus loin", class = "m-0")),
          bslib::accordion(
            multiple = TRUE,
            bslib::accordion_panel(
              "Quels indices privilégier?",
              p("Il n'existe pas de consensus absolu sur les indices d'ajustement les plus importants, mais plusieurs auteurs recommandent de consulter:"),
              tags$ul(
                tags$li(strong("La p-value du Chi-carré"), " - bien que sensible à la taille de l'échantillon"),
                tags$li(strong("RMSEA avec son intervalle de confiance"), " - un indice de parcimonie robuste"),
                tags$li(strong("TLI (ou NNFI)"), " - moins affecté par la taille de l'échantillon"),
                tags$li(strong("SRMR"), " - un indice basé sur les résidus")
              ),
              p("Il est recommandé d'examiner plusieurs indices de différentes familles plutôt que de se fier à un seul.")
            ),
            bslib::accordion_panel(
              "Références bibliographiques",
              tags$ul(
                tags$li("Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis: Conventional criteria versus new alternatives. Structural Equation Modeling, 6(1), 1-55."),
                tags$li("Brown, T. A. (2015). Confirmatory factor analysis for applied research (2nd ed.). Guilford Press."),
                tags$li("Kline, R. B. (2016). Principles and practice of structural equation modeling (4th ed.). Guilford Press.")
              )
            )
          )
        )
      )
    } else {
      # Pour les analyses en composantes principales (ACP)
      # Les indices d'ajustement classiques ne s'appliquent pas de la même façon
      div(
        h3("Indices d'ajustement du modèle ACP", class = "mt-3 mb-4"),
        
        div(class = "alert alert-info",
            h4("Note sur l'ajustement en ACP:", class = "alert-heading"),
            p("Contrairement à l'analyse factorielle, l'ACP est une méthode de réduction de dimensionnalité plutôt qu'un modèle confirmatoire. L'évaluation de son ajustement se concentre principalement sur la proportion de variance expliquée.")
        ),
        
        # Tableau des indices pertinents pour ACP
        card(
          card_header(h4("Indices d'ajustement pour ACP", class = "m-0")),
          div(class = "table-responsive",
              tags$table(class = "table table-hover",
                         tags$thead(
                           tags$tr(
                             tags$th("Indice"),
                             tags$th("Valeur"),
                             tags$th("Interprétation")
                           )
                         ),
                         tags$tbody(
                           # Variance expliquée
                           tags$tr(
                             tags$td(
                               "Variance expliquée",
                               tags$div(class = "small text-muted", "Proportion de variance expliquée")
                             ),
                             tags$td(HTML(paste0(
                               "<strong>", round(fa_results()$Vaccounted[2, input$nfactors]*100, 1), "%</strong>",
                               "<br><span class='small text-muted'>(par les ", input$nfactors, " composantes)</span>"
                             ))),
                             tags$td(
                               paste0(
                                 "Proportion de la variance des variables expliquée par les composantes. ",
                                 "Valeurs autour de 70-80% indiquent généralement une bonne solution."
                               )
                             )
                           ),
                           
                           # RMSR si disponible
                           if (!is.null(fa_results()$residual)) {
                             rmsr_value <- sqrt(mean(fa_results()$residual^2))
                             tags$tr(
                               tags$td(
                                 "RMSR",
                                 tags$div(class = "small text-muted", "Root Mean Square Residual")
                               ),
                               tags$td(round(rmsr_value, 3)),
                               tags$td("Racine carrée de la moyenne des résidus au carré. Valeurs plus basses indiquent un meilleur ajustement.")
                             )
                           }
                         )
              )
          )
        )
      )
    }
  })


  # Variance expliquée et graphique des valeurs propres
  output$var_explained <- renderPrint({
    req(fa_results())
    
    cat("VARIANCE EXPLIQUÉE PAR LES FACTEURS\n")
    cat("----------------------------------\n")
    print(fa_results()$Vaccounted)
    
    cat("\nValeurs propres des facteurs:\n")
    print(fa_results()$values[1:input$nfactors])
  })
  
  output$scree_plot <- renderPlot({
    req(fa_results())
    
    # Extraire les valeurs propres
    eigenvalues <- fa_results()$values
    
    # Créer un dataframe pour le graphique
    scree_data <- data.frame(
      Factor = 1:length(eigenvalues),
      Eigenvalue = eigenvalues
    )
    
    # Règle de Kaiser (valeurs propres > 1)
    kaiser_line <- data.frame(x = c(0, length(eigenvalues)), y = c(1, 1))
    
    # Créer le graphique
    ggplot(scree_data, aes(x = Factor, y = Eigenvalue)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_line(data = kaiser_line, aes(x = x, y = y), linetype = "dashed", color = "red") +
      geom_vline(xintercept = input$nfactors + 0.5, linetype = "dashed", color = "blue") +
      scale_x_continuous(breaks = 1:length(eigenvalues)) +
      labs(title = "Graphique des valeurs propres (Scree Plot)",
           subtitle = "Ligne rouge = Règle de Kaiser (VP > 1), Ligne bleue = Nombre choisi de facteurs",
           x = "Facteur",
           y = "Valeur propre") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))
  })
  
  # Charges factorielles
  loadings_formatted <- reactive({
    req(fa_results())
    req(input$loading_cutoff)
    
    # Extraire la matrice de charges
    if (input$rotate %in% c("promax", "oblimin", "simplimax")) {
      # Pour les rotations obliques, utiliser la matrice pattern
      loadings_mat <- fa_results()$loadings
    } else {
      # Pour les rotations orthogonales, utiliser la matrice des charges
      loadings_mat <- fa_results()$loadings
    }
    
    # Convertir en dataframe
    loadings_df <- as.data.frame(unclass(loadings_mat))
    loadings_df$Variable <- rownames(loadings_mat)
    
    # Calculer la complexité et l'unicité
    loadings_df$Complexity <- fa_results()$complexity
    loadings_df$Uniqueness <- fa_results()$uniquenesses
    loadings_df$Communality <- 1 - loadings_df$Uniqueness
    
    # Tri selon l'option sélectionnée
    if (input$sort_loadings == "size") {
      # Trier par la plus grande charge absolue
      max_loading <- apply(loadings_mat, 1, function(x) max(abs(x)))
      loadings_df <- loadings_df[order(max_loading, decreasing = TRUE), ]
    } else if (input$sort_loadings == "factor") {
      # Trier par facteur principal
      primary_factor <- apply(loadings_mat, 1, function(x) which.max(abs(x)))
      secondary_val <- apply(loadings_mat, 1, function(x) -max(abs(x)))
      loadings_df <- loadings_df[order(primary_factor, secondary_val), ]
    }
    
    # Réorganiser les colonnes
    col_order <- c("Variable", colnames(loadings_df)[1:ncol(loadings_mat)], 
                   "Communality", "Uniqueness", "Complexity")
    loadings_df <- loadings_df[, col_order]
    
    # Formater les charges (masquer celles inférieures au seuil)
    loadings_clean <- loadings_df
    for (i in 2:(ncol(loadings_mat) + 1)) {
      loadings_clean[abs(loadings_clean[, i]) < input$loading_cutoff, i] <- NA
    }
    
    return(list(full = loadings_df, clean = loadings_clean))
  })
  
  # Modification du rendu de output$loadings_table
  output$loadings_table <- renderDT({
    req(loadings_formatted())
    
    # Obtenir les données formatées
    loadings_data <- loadings_formatted()$clean
    
    # Préparer un dataframe de base pour styliser
    loadings_styled <- loadings_data
    
    # Identifier les facteurs pour chaque variable
    factor_columns <- 2:(input$nfactors + 1)
    var_factors <- apply(loadings_styled[, factor_columns], 1, function(row) {
      if(all(is.na(row))) return(NA)
      which.max(abs(row))
    })
    
    # Créer un dataTable
    dt <- datatable(
      loadings_styled, 
      options = list(
        scrollX = TRUE, 
        pageLength = 20,
        columnDefs = list(
          list(targets = factor_columns, className = "dt-center")
        )
      ),
      rownames = FALSE
    ) %>%
      formatRound(columns = 2:(ncol(loadings_styled)), digits = 3)
    
    # Mettre en évidence cas Heywood (charges > 1.0)
    if(input$show_heywood) {
      for(i in factor_columns) {
        dt <- dt %>% formatStyle(
          columns = i,
          backgroundColor = styleInterval(c(1.0), c("white", "#ffcccc")),
          fontWeight = styleInterval(c(1.0), c("normal", "bold"))
        )
      }
    }
    
    # Identifier les variables complexes
    if(input$show_complex) {
      # Calculer le nombre de charges significatives par variable
      significant_loadings <- apply(loadings_data[, factor_columns], 1, function(row) {
        sum(!is.na(row))
      })
      
      # Mettre en évidence les variables qui chargent sur plusieurs facteurs
      dt <- dt %>% formatStyle(
        columns = "Variable",
        backgroundColor = styleEqual(
          loadings_data$Variable[significant_loadings > 1], 
          rep("#fff3cd", sum(significant_loadings > 1))
        ),
        fontWeight = styleEqual(
          loadings_data$Variable[significant_loadings > 1], 
          rep("bold", sum(significant_loadings > 1))
        )
      )
    }
    
    # Identifier les charges croisées
    if(input$highlight_crossloadings) {
      # Pour chaque variable, créer un vecteur de styles pour les cellules correspondantes
      for(row_idx in 1:nrow(loadings_data)) {
        row_values <- unlist(loadings_data[row_idx, factor_columns])
        row_values <- row_values[!is.na(row_values)]
        
        if(length(row_values) > 1) {
          # Trier les valeurs absolues
          sorted_abs <- sort(abs(row_values), decreasing = TRUE)
          
          # Vérifier si la différence entre les deux plus grandes charges est inférieure au seuil
          if(sorted_abs[1] - sorted_abs[2] < input$crossloading_diff) {
            # Trouver les colonnes des charges croisées
            cross_cols <- factor_columns[!is.na(loadings_data[row_idx, factor_columns])]
            
            # Appliquer le style à chaque colonne croisée
            for(col in cross_cols) {
              dt <- dt %>% formatStyle(
                columns = col,
                target = "cell",
                selector = paste0('tr:eq(', row_idx - 1, ')'),
                backgroundColor = "#e6f7ff"
              )
            }
          }
        }
      }
    }
    
    # Mettre en couleur par facteur
    if(input$highlight_factors) {
      # Générer une palette de couleurs pour les facteurs
      factor_colors <- RColorBrewer::brewer.pal(min(input$nfactors, 9), "Set1")
      if(input$nfactors > 9) {
        # Étendre la palette si nécessaire
        factor_colors <- colorRampPalette(factor_colors)(input$nfactors)
      }
      
      # Pour chaque variable, colorer la cellule correspondant à sa charge principale
      for(row_idx in 1:nrow(loadings_data)) {
        # Trouver le facteur dominant pour cette variable
        dominant_factor <- var_factors[row_idx]
        
        if(!is.na(dominant_factor)) {
          # Appliquer un style à la cellule correspondante
          col_idx <- factor_columns[dominant_factor]
          
          # Créer une couleur avec transparence (20%)
          bg_color <- paste0(substr(factor_colors[dominant_factor], 1, 7), "33")
          
          dt <- dt %>% formatStyle(
            columns = col_idx,
            target = "cell",
            selector = paste0('tr:eq(', row_idx - 1, ')'),
            backgroundColor = bg_color
          )
        }
      }
    }
    
    return(dt)
  })
  
  output$loadings_plot <- renderPlot({
    req(fa_results())
    req(input$loading_cutoff)
    
    # Extraire les charges pour le graphique
    if (input$nfactors < 2) {
      return(NULL) # Ne pas afficher le graphique si un seul facteur
    }
    
    # Version compatible avec les versions récentes et anciennes de psych
    if (packageVersion("psych") >= "2.0.6") {
      # Pour psych >= 2.0.6, utiliser l'argument loading.method
      fa.diagram(fa_results(), 
                 cut = input$loading_cutoff, 
                 simple = TRUE,
                 digits = 2,
                 main = "Diagramme Factoriel",
                 rsize = 0.15,
                 marg = c(0.5, 1, 0.5, 1),
                 adj = 4,
                 cex = 1.2,
                 l.cex = 0.9,
                 loading.method = "original")
    } else {
      # Pour les anciennes versions de psych
      fa.diagram(fa_results(), 
                 cut = input$loading_cutoff, 
                 simple = TRUE,
                 digits = 2,
                 main = "Diagramme Factoriel",
                 rsize = 0.15,
                 marg = c(0.5, 1, 0.5, 1),
                 adj = 4,
                 cex = 1.2,
                 l.cex = 0.9)
    }
  })
  
  observeEvent(input$show_lavaan_syntax, {
    req(fa_results())
    req(loadings_formatted())
    
    # N'afficher que pour les analyses factorielles, pas pour les ACP
    if(input$analysis_type != "fa") {
      showModal(modalDialog(
        title = "Syntaxe lavaan non disponible",
        "La syntaxe lavaan est uniquement disponible pour les analyses factorielles, pas pour les ACP.",
        easyClose = TRUE,
        footer = modalButton("Fermer")
      ))
      return()
    }
    
    # Générer la syntaxe lavaan avec coefficients
    syntax_with_coef <- capture.output({
      cat("# Modèle d'analyse factorielle en syntaxe lavaan (avec coefficients)\n")
      cat("# Cette version peut être utilisée pour la simulation de données\n\n")
      
      # Pour chaque facteur, créer la définition
      loadings <- loadings_formatted()$clean
      
      for (i in 1:input$nfactors) {
        factor_name <- paste0("F", i)
        
        # Récupérer les variables qui chargent sur ce facteur
        factor_col <- i + 1  # +1 car la première colonne est "Variable"
        variables <- loadings$Variable[!is.na(loadings[, factor_col])]
        
        if (length(variables) > 0) {
          cat(factor_name, " =~ ")
          
          # Ajouter les variables avec leurs charges
          for (j in 1:length(variables)) {
            var_name <- variables[j]
            loading_val <- loadings[loadings$Variable == var_name, factor_col]
            
            # Formater la charge pour plus de clarté
            loading_str <- sprintf("%.3f", loading_val)
            
            # Ajouter la variable à la définition du facteur
            cat(loading_str, "*", var_name)
            
            # Ajouter le séparateur approprié
            if (j < length(variables)) {
              cat(" + ")
            } else {
              cat("\n")
            }
          }
        } else {
          cat(factor_name, " =~ 0\n")
        }
      }
      
      # Si rotation oblique, ajouter les corrélations entre facteurs
      if (input$rotate %in% c("promax", "oblimin", "simplimax") && input$nfactors > 1) {
        cat("\n# Corrélations entre facteurs\n")
        
        phi_matrix <- fa_results()$Phi
        for (i in 1:(input$nfactors-1)) {
          for (j in (i+1):input$nfactors) {
            corr_value <- phi_matrix[i, j]
            if (!is.na(corr_value)) {
              cat(sprintf("F%d ~~ %.3f*F%d\n", i, corr_value, j))
            }
          }
        }
      }
    })
    
    # Générer la syntaxe lavaan sans coefficients
    syntax_without_coef <- capture.output({
      cat("# Modèle d'analyse factorielle en syntaxe lavaan (sans coefficients)\n")
      cat("# Cette version peut être utilisée pour l'analyse de données\n\n")
      
      # Pour chaque facteur, créer la définition
      loadings <- loadings_formatted()$clean
      
      for (i in 1:input$nfactors) {
        factor_name <- paste0("F", i)
        
        # Récupérer les variables qui chargent sur ce facteur
        factor_col <- i + 1  # +1 car la première colonne est "Variable"
        variables <- loadings$Variable[!is.na(loadings[, factor_col])]
        
        if (length(variables) > 0) {
          cat(factor_name, " =~ ")
          
          # Ajouter les variables sans leurs charges
          for (j in 1:length(variables)) {
            var_name <- variables[j]
            
            # Ajouter la variable à la définition du facteur
            cat(var_name)
            
            # Ajouter le séparateur approprié
            if (j < length(variables)) {
              cat(" + ")
            } else {
              cat("\n")
            }
          }
        } else {
          cat(factor_name, " =~ 0\n")
        }
      }
      
      # Si rotation oblique, ajouter les corrélations entre facteurs sans coefficients
      if (input$rotate %in% c("promax", "oblimin", "simplimax") && input$nfactors > 1) {
        cat("\n# Corrélations entre facteurs\n")
        
        phi_matrix <- fa_results()$Phi
        for (i in 1:(input$nfactors-1)) {
          for (j in (i+1):input$nfactors) {
            corr_value <- phi_matrix[i, j]
            if (!is.na(corr_value) && abs(corr_value) > 0.1) { # Seuil minimal pour inclure une corrélation
              cat(sprintf("F%d ~~ F%d\n", i, j))
            }
          }
        }
      }
    })
    
    # Créer et afficher la boîte de dialogue modale avec onglets
    showModal(modalDialog(
      title = "Syntaxe lavaan du modèle factoriel",
      tabsetPanel(
        tabPanel("Avec coefficients", 
                 pre(paste(syntax_with_coef, collapse = "\n"))
        ),
        tabPanel("Sans coefficients", 
                 pre(paste(syntax_without_coef, collapse = "\n"))
        )
      ),
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        downloadButton("download_lavaan_with_coef", "Télécharger (avec coef)", class = "btn-info"),
        downloadButton("download_lavaan_without_coef", "Télécharger (sans coef)", class = "btn-outline-info"),
        modalButton("Fermer")
      )
    ))
  })
  
  # Gestion des téléchargements pour les deux types de syntaxe
  output$download_lavaan_with_coef <- downloadHandler(
    filename = function() {
      paste("lavaan_model_with_coef_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Régénérer la syntaxe avec coefficients
      syntax <- capture.output({
        cat("# Modèle d'analyse factorielle en syntaxe lavaan (avec coefficients)\n")
        cat("# Cette version peut être utilisée pour la simulation de données\n\n")
        
        # Pour chaque facteur, créer la définition
        loadings <- loadings_formatted()$clean
        
        for (i in 1:input$nfactors) {
          factor_name <- paste0("F", i)
          
          # Récupérer les variables qui chargent sur ce facteur
          factor_col <- i + 1  # +1 car la première colonne est "Variable"
          variables <- loadings$Variable[!is.na(loadings[, factor_col])]
          
          if (length(variables) > 0) {
            cat(factor_name, " =~ ")
            
            # Ajouter les variables avec leurs charges
            for (j in 1:length(variables)) {
              var_name <- variables[j]
              loading_val <- loadings[loadings$Variable == var_name, factor_col]
              
              # Formater la charge pour plus de clarté
              loading_str <- sprintf("%.3f", loading_val)
              
              # Ajouter la variable à la définition du facteur
              cat(loading_str, "*", var_name)
              
              # Ajouter le séparateur approprié
              if (j < length(variables)) {
                cat(" + ")
              } else {
                cat("\n")
              }
            }
          } else {
            cat(factor_name, " =~ 0\n")
          }
        }
        
        # Si rotation oblique, ajouter les corrélations entre facteurs
        if (input$rotate %in% c("promax", "oblimin", "simplimax") && input$nfactors > 1) {
          cat("\n# Corrélations entre facteurs\n")
          
          phi_matrix <- fa_results()$Phi
          for (i in 1:(input$nfactors-1)) {
            for (j in (i+1):input$nfactors) {
              corr_value <- phi_matrix[i, j]
              if (!is.na(corr_value)) {
                cat(sprintf("F%d ~~ %.3f*F%d\n", i, corr_value, j))
              }
            }
          }
        }
      })
      
      writeLines(syntax, file)
    }
  )
  
  output$download_lavaan_without_coef <- downloadHandler(
    filename = function() {
      paste("lavaan_model_without_coef_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Régénérer la syntaxe sans coefficients
      syntax <- capture.output({
        cat("# Modèle d'analyse factorielle en syntaxe lavaan (sans coefficients)\n")
        cat("# Cette version peut être utilisée pour l'analyse de données\n\n")
        
        # Pour chaque facteur, créer la définition
        loadings <- loadings_formatted()$clean
        
        for (i in 1:input$nfactors) {
          factor_name <- paste0("F", i)
          
          # Récupérer les variables qui chargent sur ce facteur
          factor_col <- i + 1  # +1 car la première colonne est "Variable"
          variables <- loadings$Variable[!is.na(loadings[, factor_col])]
          
          if (length(variables) > 0) {
            cat(factor_name, " =~ ")
            
            # Ajouter les variables sans leurs charges
            for (j in 1:length(variables)) {
              var_name <- variables[j]
              
              # Ajouter la variable à la définition du facteur
              cat(var_name)
              
              # Ajouter le séparateur approprié
              if (j < length(variables)) {
                cat(" + ")
              } else {
                cat("\n")
              }
            }
          } else {
            cat(factor_name, " =~ 0\n")
          }
        }
        
        # Si rotation oblique, ajouter les corrélations entre facteurs sans coefficients
        if (input$rotate %in% c("promax", "oblimin", "simplimax") && input$nfactors > 1) {
          cat("\n# Corrélations entre facteurs\n")
          
          phi_matrix <- fa_results()$Phi
          for (i in 1:(input$nfactors-1)) {
            for (j in (i+1):input$nfactors) {
              corr_value <- phi_matrix[i, j]
              if (!is.na(corr_value) && abs(corr_value) > 0.1) { # Seuil minimal pour inclure une corrélation
                cat(sprintf("F%d ~~ F%d\n", i, j))
              }
            }
          }
        }
      })
      
      writeLines(syntax, file)
    }
  )
  
  # Ajout du handler pour le téléchargement de la syntaxe lavaan
  output$download_lavaan <- downloadHandler(
    filename = function() {
      paste("lavaan_model_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Récupérer les chargements formatés
      loadings <- loadings_formatted()$clean
      
      # Ouvrir une connexion de fichier
      con <- file(file, "w")
      
      # Écrire le modèle lavaan
      cat("# Modèle d'analyse factorielle en syntaxe lavaan\n\n", file = con)
      
      # Pour chaque facteur, créer la définition
      for (i in 1:input$nfactors) {
        factor_name <- paste0("F", i)
        
        # Récupérer les variables qui chargent sur ce facteur
        factor_col <- i + 1  # +1 car la première colonne est "Variable"
        variables <- loadings$Variable[!is.na(loadings[, factor_col])]
        
        if (length(variables) > 0) {
          cat(factor_name, " =~ ", file = con)
          
          # Ajouter les variables avec leurs charges
          for (j in 1:length(variables)) {
            var_name <- variables[j]
            loading_val <- loadings[loadings$Variable == var_name, factor_col]
            
            # Formater la charge pour plus de clarté
            loading_str <- sprintf("%.3f", loading_val)
            
            # Ajouter la variable à la définition du facteur
            cat(loading_str, "*", var_name, file = con)
            
            # Ajouter le séparateur approprié
            if (j < length(variables)) {
              cat(" + ", file = con)
            } else {
              cat("\n", file = con)
            }
          }
        } else {
          cat(factor_name, " =~ 0\n", file = con)
        }
      }
      
      # Si rotation oblique, ajouter les corrélations entre facteurs
      if (input$rotate %in% c("promax", "oblimin", "simplimax") && input$nfactors > 1) {
        cat("\n# Corrélations entre facteurs\n", file = con)
        
        phi_matrix <- fa_results()$Phi
        for (i in 1:(input$nfactors-1)) {
          for (j in (i+1):input$nfactors) {
            corr_value <- phi_matrix[i, j]
            if (!is.na(corr_value)) {
              cat(sprintf("F%d ~~ %.3f*F%d\n", i, corr_value, j), file = con)
            }
          }
        }
      }
      
      # Fermer la connexion
      close(con)
    }
  )
  
  # Syntaxe lavaan pour le modèle d'analyse factorielle
  output$lavaan_syntax2 <- renderPrint({
    req(fa_results())
    req(loadings_formatted())
    
    # Extraire les charges factorielles nettoyées
    loadings <- loadings_formatted()$clean
    
    # Générer la syntaxe lavaan
    cat("# Modèle d'analyse factorielle en syntaxe lavaan\n\n")
    
    # Pour chaque facteur, créer la définition
    for (i in 1:input$nfactors) {
      factor_name <- paste0("F", i)
      
      # Récupérer les variables qui chargent sur ce facteur
      factor_col <- i + 1  # +1 car la première colonne est "Variable"
      variables <- loadings$Variable[!is.na(loadings[, factor_col])]
      
      if (length(variables) > 0) {
        cat(factor_name, " =~ ")
        
        # Ajouter les variables avec leurs charges
        for (j in 1:length(variables)) {
          var_name <- variables[j]
          loading_val <- loadings[loadings$Variable == var_name, factor_col]
          
          # Formater la charge pour plus de clarté
          loading_str <- sprintf("%.3f", loading_val)
          
          # Ajouter la variable à la définition du facteur
          cat(loading_str, "*", var_name)
          
          # Ajouter le séparateur approprié
          if (j < length(variables)) {
            cat(" + ")
          } else {
            cat("\n")
          }
        }
      } else {
        cat(factor_name, " =~ 0\n")
      }
    }
    
    # Si rotation oblique, ajouter les corrélations entre facteurs
    if (input$rotate %in% c("promax", "oblimin", "simplimax") && input$nfactors > 1) {
      cat("\n# Corrélations entre facteurs\n")
      
      phi_matrix <- fa_results()$Phi
      for (i in 1:(input$nfactors-1)) {
        for (j in (i+1):input$nfactors) {
          corr_value <- phi_matrix[i, j]
          if (!is.na(corr_value)) {
            cat(sprintf("F%d ~~ %.3f*F%d\n", i, corr_value, j))
          }
        }
      }
    }
  })
  
  # Structure factorielle
  output$factor_structure <- renderUI({
    req(fa_results())
    
    # Déterminer si nous avons une rotation oblique ou orthogonale
    if (input$rotate %in% c("promax", "oblimin", "simplimax")) {
      # Pour les rotations obliques (afficher la matrice pattern et la matrice de structure côte à côte)
      layout_columns(
        col_widths = c(6, 6),
        
        # Card pour la matrice pattern
        card(
          card_header("Matrice pattern"),
          div(
            style = "padding: 10px;",
            p("La matrice pattern contient les coefficients de régression standardisés des variables sur les facteurs."),
            flextable_from_loadings(fa_results()$loadings, input$loading_cutoff, "pattern")
          )
        ),
        
        # Card pour la matrice de structure
        card(
          card_header("Matrice de structure"),
          div(
            style = "padding: 10px;",
            p("La matrice de structure contient les corrélations entre variables et facteurs."),
            flextable_from_loadings(fa_results()$Structure, input$loading_cutoff, "structure")
          )
        )
      )
    } else {
      # Pour les rotations orthogonales (afficher uniquement la matrice des charges dans une card)
      card(
        card_header("Matrice des charges factorielles"),
        div(
          style = "padding: 10px;",
          p("Dans une rotation orthogonale, les charges factorielles peuvent être interprétées à la fois comme des coefficients de régression standardisés et comme des corrélations entre variables et facteurs."),
          flextable_from_loadings(fa_results()$loadings, input$loading_cutoff, "loadings")
        )
      )
    }
  })
  
  # Fonction helper pour créer des flextables à partir des matrices de charges factorielles
  flextable_from_loadings <- function(loadings_mat, cutoff = 0.3, type = "loadings") {
    # Vérifier si la matrice existe et est valide
    if (is.null(loadings_mat) || !is.matrix(loadings_mat)) {
      # Retourner un message d'erreur formaté si la matrice est invalide
      return(htmltools::div(
        style = "color: red; padding: 15px; border: 1px solid #f8d7da; background-color: #f8d7da; border-radius: 5px;",
        htmltools::h4("Erreur de chargement des données"),
        htmltools::p("Impossible d'afficher la matrice des charges factorielles. Vérifiez vos paramètres d'analyse.")
      ))
    }
    
    # Convertir la matrice de charges en dataframe
    loadings_df <- as.data.frame(unclass(loadings_mat))
    
    # Ajouter une colonne pour les noms de variables
    loadings_df$Variable <- rownames(loadings_mat)
    
    # Réorganiser les colonnes (Variable en premier)
    loadings_df <- loadings_df[, c("Variable", setdiff(names(loadings_df), "Variable"))]
    
    # Créer le flextable
    ft <- flextable(loadings_df) %>%
      theme_booktabs() %>%
      autofit()
    
    # Formater les valeurs numériques avec 3 décimales
    ft <- ft %>% colformat_double(j = 2:ncol(loadings_df), digits = 3)
    
    # Identifier les indices des lignes où les valeurs sont >= cutoff (en valeur absolue)
    for (i in 2:ncol(loadings_df)) {
      col_name <- names(loadings_df)[i]
      
      # Identifier les indices des lignes où les valeurs sont >= cutoff (en valeur absolue)
      sig_indices <- which(abs(loadings_df[[col_name]]) >= cutoff)
      nonsig_indices <- which(abs(loadings_df[[col_name]]) < cutoff)
      
      # Appliquer le style pour les valeurs significatives
      if (length(sig_indices) > 0) {
        ft <- ft %>%
          bold(i = sig_indices, j = i) %>%
          bg(i = sig_indices, j = i, bg = "#e6f3ff")
      }
      
      # Appliquer le style pour les valeurs non significatives
      if (length(nonsig_indices) > 0) {
        ft <- ft %>%
          color(i = nonsig_indices, j = i, color = "#888888")
      }
    }
    
    # Styles d'en-tête différents selon le type de matrice
    header_bg <- switch(type,
                        "pattern" = "#daeeff",
                        "structure" = "#ffedd8",
                        "loadings" = "#e2f4e5")
    
    # Appliquer des styles d'en-tête
    ft <- ft %>%
      bg(part = "header", bg = header_bg) %>%
      bold(part = "header") %>%
      align(part = "header", align = "center") %>%
      align(j = 1, align = "left") %>%
      align(j = 2:ncol(loadings_df), align = "right") %>%
      add_footer_lines(paste("Les charges factorielles supérieures à", cutoff, "sont affichées en bleu."))
    
    # Appliquer des bordures correctement, sans utiliser l'argument i pour hline_bottom
    ft <- ft %>%
      border_remove() %>%
      hline_top(part = "header", border = officer::fp_border(color = "gray30")) %>%
      hline_bottom(part = "header", border = officer::fp_border(color = "gray30")) %>%
      hline_bottom(part = "body", border = officer::fp_border(color = "gray30"))
    
    # Retourner le flextable comme HTML pour l'affichage
    flextable::htmltools_value(ft)
  }
  
  # Corrélations entre facteurs
  output$factor_cors_plot <- renderPlot({
    req(fa_results())
    
    # Vérifier si nous avons des corrélations entre facteurs
    if (!input$rotate %in% c("promax", "oblimin", "simplimax") || 
        is.null(fa_results()$Phi) || input$nfactors < 2) {
      return(NULL)
    }
    
    # Créer une matrice de corrélation
    cor_mat <- fa_results()$Phi
    
    # Afficher avec corrplot
    corrplot(cor_mat, method = "color", 
             type = "full", 
             tl.col = "black", 
             addCoef.col = "black",
             title = "Corrélations entre facteurs",
             mar = c(0, 0, 1, 0))
  })
  
  output$factor_cors_table <- renderDT({
    req(fa_results())
    
    # Vérifier si nous avons des corrélations entre facteurs
    if (!input$rotate %in% c("promax", "oblimin", "simplimax") || 
        is.null(fa_results()$Phi) || input$nfactors < 2) {
      return(NULL)
    }
    
    # Convertir la matrice en dataframe
    cor_df <- as.data.frame(fa_results()$Phi)
    rownames(cor_df) <- colnames(cor_df) <- paste0("Facteur ", 1:input$nfactors)
    
    # Afficher le tableau
    datatable(cor_df, options = list(pageLength = 10)) %>%
      formatRound(columns = 1:ncol(cor_df), digits = 3)
  })
  
  # Résidus
  output$residuals_plot <- renderPlot({
    req(fa_results())
    
    # Extraire les résidus
    residuals <- fa_results()$residual
    
    # Convertir en dataframe pour ggplot
    resid_df <- melt(residuals)
    colnames(resid_df) <- c("Var1", "Var2", "Residual")
    
    # Créer le graphique
    ggplot(resid_df, aes(x = Var1, y = Var2, fill = Residual)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, limits = max(abs(residuals)) * c(-1, 1)) +
      labs(title = "Matrice des résidus",
           subtitle = "Rouge = résidu positif, Bleu = résidu négatif",
           x = "", y = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold"))
  })
  
  output$residuals_table <- renderDT({
    req(fa_results())
    
    # Créer un dataframe des résidus
    resid_df <- as.data.frame(fa_results()$residual)
    
    # Afficher le tableau
    datatable(resid_df, options = list(scrollX = TRUE, pageLength = 10)) %>%
      formatRound(columns = 1:ncol(resid_df), digits = 3)
  })
  
  # Communalités
  output$communalities_plot <- renderPlot({
    req(fa_results())
    
    # Créer un dataframe pour le graphique
    comm_df <- data.frame(
      Variable = names(fa_results()$communality),
      Communality = fa_results()$communality,
      Uniqueness = fa_results()$uniquenesses
    )
    
    # Trier par communalité
    comm_df <- comm_df[order(comm_df$Communality, decreasing = TRUE), ]
    
    # Convertir en format long pour ggplot
    comm_long <- reshape2::melt(comm_df, id.vars = "Variable")
    
    # Créer le graphique
    ggplot(comm_long, aes(x = reorder(Variable, value), y = value, fill = variable)) +
      geom_bar(stat = "identity", position = "stack") +
      coord_flip() +
      scale_fill_manual(values = c("Communality" = "steelblue", "Uniqueness" = "lightgrey")) +
      labs(title = "Communalités et unicités des variables",
           x = "", y = "Proportion de variance", fill = "") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")
  })
  
  output$communalities_table <- renderDT({
    req(fa_results())
    
    # Créer un dataframe
    comm_df <- data.frame(
      Variable = names(fa_results()$communality),
      Communality = fa_results()$communality,
      Uniqueness = fa_results()$uniquenesses
    )
    
    # Trier par communalité
    comm_df <- comm_df[order(comm_df$Communality, decreasing = TRUE), ]
    
    # Afficher le tableau
    datatable(comm_df, options = list(pageLength = 20), rownames = FALSE) %>%
      formatRound(columns = c("Communality", "Uniqueness"), digits = 3)
  })
  
  # Résumé complet
  output$full_summary <- renderPrint({
    req(fa_results())
    print(fa_results())
  })
  
  # Fonction pour calculer les network loads - version unifiée
  netloads_result <- reactive({
    req(ega_results())
    
    # Vérifier si l'onglet Network Loads est actif selon le type de données
    tab_active <- if (input$data_type == "raw") {
      input$ega_inner_tabs == "Network Loads"
    } else {
      input$ega_inner_tabs_cormat == "Network Loads"
    }
    
    if (!tab_active) {
      return(NULL)
    }
    
    withProgress(message = 'Calcul automatique des network loads...', {
      tryCatch({
        loads <- net.loads(
          A = ega_results(),
          wc = ega_results()$wc
        )
        
        return(loads)
      }, error = function(e) {
        showNotification(paste("Erreur lors du calcul des network loads:", e$message), type = "error")
        return(NULL)
      })
    })
  })
  
  output$ega_network_matrix <- renderDT({
    req(ega_results())
    
    # Obtenir la matrice du réseau
    network_matrix <- ega_results()$network
    
    # Convertir en format triangulaire inférieur (sans diagonale)
    network_lower <- network_matrix
    network_lower[upper.tri(network_lower, diag = TRUE)] <- NA
    
    # Convertir en dataframe pour l'affichage
    network_df <- as.data.frame(network_lower)
    
    # Créer un datatable avec la matrice triangulaire inférieure
    datatable(network_df, 
              options = list(
                scrollX = TRUE, 
                pageLength = 20,
                searching = FALSE,  # Désactiver la recherche
                lengthChange = FALSE  # Désactiver le changement de nombre d'éléments par page
              )) %>%
      formatRound(columns = 1:ncol(network_df), digits = 4) %>%
      formatStyle(
        columns = 1:ncol(network_df),
        backgroundColor = styleInterval(
          c(-0.5, -0.2, 0.2, 0.5), 
          c("#d73027", "#f46d43", "#ffffff", "#74add1", "#4575b4")
        )
      )
  })
  
  # UI pour afficher les network loads
  output$netloads_ui <- renderUI({
    req(ega_results())
    
    # Si les network loads sont calculés, les afficher
    if (!is.null(netloads_result())) {
      loads <- netloads_result()
      
      tagList(
        hr(),
        h4("Network Loads"),
        tabsetPanel(
          tabPanel("Loads standardisés", 
                   DTOutput("netloads_std_table")
          ),
          tabPanel("Loads non-standardisés", 
                   DTOutput("netloads_unstd_table")
          ),
          tabPanel("Loads pivotés", 
                   if(!is.null(loads$rotated)) {
                     tagList(
                       h5("Charges pivotées"),
                       DTOutput("netloads_rotated_table"),
                       h5("Corrélations entre dimensions"),
                       DTOutput("netloads_phi_table")
                     )
                   } else {
                     p("Aucune rotation n'a été appliquée. Choisissez une méthode de rotation dans les paramètres pour afficher les charges pivotées.")
                   }
          )
        )
      )
    } else {
      # Message d'attente au cas où le calcul prend du temps
      tagList(
        div(
          style = "text-align: center; padding: 30px;",
          h4("Calcul des Network Loads en cours..."),
          p("Les network loads sont calculés automatiquement. Veuillez patienter quelques instants.")
        )
      )
    }
  })
  
  # Rendu des tableaux
  output$netloads_std_table <- renderDT({
    req(netloads_result())
    datatable(netloads_result()$std, 
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatRound(columns = 1:ncol(netloads_result()$std), digits = 3)
  })
  
  output$netloads_unstd_table <- renderDT({
    req(netloads_result())
    datatable(netloads_result()$unstd, 
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatRound(columns = 1:ncol(netloads_result()$unstd), digits = 3)
  })
  
  output$netloads_rotated_table <- renderDT({
    req(netloads_result())
    req(netloads_result()$rotated)
    
    datatable(netloads_result()$rotated$loadings, 
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatRound(columns = 1:ncol(netloads_result()$rotated$loadings), digits = 3)
  })
  
  output$netloads_phi_table <- renderDT({
    req(netloads_result())
    req(netloads_result()$rotated)
    
    datatable(netloads_result()$rotated$Phi, 
              options = list(scrollX = TRUE, pageLength = 10)) %>%
      formatRound(columns = 1:ncol(netloads_result()$rotated$Phi), digits = 3)
  })
  
  # Handler pour télécharger les network loads
  output$download_netloads <- downloadHandler(
    filename = function() {
      paste("network_loads_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      req(netloads_result())
      
      # Créer un workbook Excel
      wb <- openxlsx::createWorkbook()
      
      # Ajouter les feuilles
      openxlsx::addWorksheet(wb, "Loads standardisés")
      openxlsx::writeData(wb, "Loads standardisés", netloads_result()$std)
      
      openxlsx::addWorksheet(wb, "Loads non-standardisés")
      openxlsx::writeData(wb, "Loads non-standardisés", netloads_result()$unstd)
      
      if(!is.null(netloads_result()$rotated)) {
        openxlsx::addWorksheet(wb, "Loads pivotés")
        openxlsx::writeData(wb, "Loads pivotés", netloads_result()$rotated$loadings)
        
        openxlsx::addWorksheet(wb, "Corrélations entre dimensions")
        openxlsx::writeData(wb, "Corrélations entre dimensions", netloads_result()$rotated$Phi)
      }
      
      # Sauvegarder le workbook
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # Fonction améliorée pour calculer le KMO et le test de Bartlett
  adequacy_tests_data <- reactive({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Obtenir les données d'analyse
    data <- analysis_data()
    n <- nrow(data)
    p <- ncol(data)
    
    # Calculer la matrice de corrélation
    cor_mat <- cor(data, use = "pairwise.complete.obs")
    
    # Utiliser psych::KMO pour calculer l'indice KMO
    kmo_result <- psych::KMO(r = cor_mat)
    
    # Si les résultats du test de Bartlett sont NULL, les calculer manuellement
    if (is.null(kmo_result$chi.sq) || is.null(kmo_result$p.value) || is.null(kmo_result$dof)) {
      # Calculer le déterminant de la matrice de corrélation
      det_cor <- det(cor_mat)
      
      # Calculer la statistique de Bartlett
      chi_square <- -((n - 1) - (2 * p + 5) / 6) * log(det_cor)
      
      # Degrés de liberté
      df <- p * (p - 1) / 2
      
      # Calculer la p-valeur
      p_value <- pchisq(chi_square, df, lower.tail = FALSE)
      
      # Ajouter les résultats de Bartlett à la liste des résultats
      bartlett_chi <- chi_square
      bartlett_df <- df
      bartlett_p <- p_value
    } else {
      # Utiliser les résultats de psych::KMO si disponibles
      bartlett_chi <- kmo_result$chi.sq
      bartlett_df <- kmo_result$dof
      bartlett_p <- kmo_result$p.value
    }
    
    # Retourner les résultats
    return(list(
      kmo = kmo_result$MSA,  # Indice KMO global
      bartlett = bartlett_chi,  # Chi-carré de Bartlett
      bartlett_df = bartlett_df,  # Degrés de liberté
      bartlett_p = bartlett_p,  # P-valeur
      n_obs = n,  # Nombre d'observations
      n_vars = p,  # Nombre de variables
      msa = kmo_result$MSAi  # Mesures MSA individuelles
    ))
  })
  
  # Affichage des tests d'adéquation
  output$adequacy_tests <- renderUI({
    req(adequacy_tests_data())
    
    # Extraire les résultats
    results <- adequacy_tests_data()
    
    # Interpréter l'indice KMO
    kmo_interp <- function(kmo) {
      if (kmo >= 0.9) return(list(text = "Excellente", color = "#28a745"))
      if (kmo >= 0.8) return(list(text = "Bonne", color = "#28a745"))
      if (kmo >= 0.7) return(list(text = "Moyenne", color = "#ffc107"))
      if (kmo >= 0.6) return(list(text = "Médiocre", color = "#fd7e14"))
      if (kmo >= 0.5) return(list(text = "Faible", color = "#dc3545"))
      return(list(text = "Inacceptable", color = "#dc3545"))
    }
    
    kmo_rating <- kmo_interp(results$kmo)
    
    div(
      h4("Tests d'adéquation pour l'analyse factorielle", class = "mb-4"),
      
      # Card pour le KMO
      card(
        card_header("Indice KMO (Kaiser-Meyer-Olkin)"),
        div(class = "p-3",
            layout_columns(
              col_widths = c(8, 4),
              div(
                p("L'indice KMO mesure l'adéquation de l'échantillonnage pour l'analyse factorielle. Les valeurs supérieures à 0.6 sont généralement considérées comme acceptables, et celles supérieures à 0.8 comme bonnes."),
                p("Interprétation conventionnelle :"),
                tags$ul(
                  tags$li("0.9 - 1.0 : Excellente"),
                  tags$li("0.8 - 0.9 : Bonne"),
                  tags$li("0.7 - 0.8 : Moyenne"),
                  tags$li("0.6 - 0.7 : Médiocre"),
                  tags$li("0.5 - 0.6 : Faible"),
                  tags$li("< 0.5 : Inacceptable")
                )
              ),
              div(
                style = "display: flex; flex-direction: column; align-items: center; justify-content: center;",
                div(
                  style = paste0("display: inline-block; padding: 25px; border-radius: 50%; background-color: ", 
                                 kmo_rating$color, "; color: white; font-size: 24px; font-weight: bold;"),
                  sprintf("%.3f", results$kmo)
                ),
                div(
                  style = "margin-top: 15px; text-align: center;",
                  tags$strong("Adéquation ", kmo_rating$text)
                )
              )
            )
        )
      ),
      
      # Card pour le test de Bartlett
      card(
        class = "mt-4",
        card_header("Test de sphéricité de Bartlett"),
        div(class = "p-3",
            layout_columns(
              col_widths = c(8, 4),
              div(
                p("Le test de sphéricité de Bartlett teste l'hypothèse que la matrice de corrélation est une matrice identité, indiquant que les variables ne sont pas corrélées."),
                p("Une p-valeur significative (< 0.05) rejette cette hypothèse et suggère que les données sont appropriées pour une analyse factorielle.")
              ),
              div(
                style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px;",
                div(style = "margin-bottom: 10px;",
                    tags$strong("Chi-carré : "), sprintf("%.2f", results$bartlett)
                ),
                div(style = "margin-bottom: 10px;",
                    tags$strong("Degrés de liberté : "), results$bartlett_df
                ),
                div(style = "margin-bottom: 10px;",
                    tags$strong("P-valeur : "), format.pval(results$bartlett_p, digits = 3)
                ),
                div(
                  style = paste0("margin-top: 15px; padding: 8px; border-radius: 5px; text-align: center; font-weight: bold; ",
                                 "background-color: ", ifelse(results$bartlett_p < 0.05, "#d4edda", "#f8d7da"), ";"),
                  ifelse(results$bartlett_p < 0.05, 
                         "Test significatif - Données appropriées", 
                         "Test non significatif - Données non appropriées")
                )
              )
            )
        )
      )
    )
  })
  
  # Affichage des mesures MSA individuelles
  output$msa_measures <- renderUI({
    req(adequacy_tests_data())
    
    div(
      h4("Mesures d'adéquation de l'échantillonnage (MSA) par variable"),
      p("Ces valeurs indiquent l'adéquation de chaque variable pour l'analyse factorielle. Une valeur MSA < 0.5 pour une variable suggère de l'exclure de l'analyse."),
      
      # Affichage des MSA sous forme de tableau uniquement (sans graphique)
      div(
        DTOutput("msa_table")
      )
    )
  })
  
  # Tableau MSA
  output$msa_table <- renderDT({
    req(adequacy_tests_data())
    
    # Extraire les mesures MSA
    msa_values <- adequacy_tests_data()$msa
    
    # Créer le dataframe
    msa_df <- data.frame(
      Variable = names(msa_values),
      MSA = msa_values
    )
    
    # Trier par MSA décroissant
#    msa_df <- msa_df[order(msa_df$MSA, decreasing = TRUE), ]
    
    # Créer des catégories d'adéquation
    msa_df$Adequation <- cut(msa_df$MSA, 
                             breaks = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                             labels = c("Inacceptable", "Faible", "Médiocre", "Moyenne", "Bonne", "Excellente"),
                             include.lowest = TRUE)
    
    # Créer un tableau interactif avec formatage - sans recherche
    datatable(msa_df, 
              options = list(pageLength = 15,
                             searching = FALSE,  # Désactiver la recherche
                             dom = 'tip'),      # Enlever la boîte de recherche 'f'
              rownames = FALSE) %>%
      formatRound(columns = "MSA", digits = 3) %>%
      formatStyle(
        'MSA',
        backgroundColor = styleInterval(
          c(0.5, 0.6, 0.7, 0.8, 0.9),
          c('#f8d7da', '#fff3cd', '#fffde7', '#e8f5e9', '#e3f2fd', '#e1f5fe')  # Couleurs plus douces
        )
      ) %>%
      formatStyle(
        'Adequation',
        backgroundColor = styleEqual(
          c("Inacceptable", "Faible", "Médiocre", "Moyenne", "Bonne", "Excellente"),
          c('#f8d7da', '#fff3cd', '#fffde7', '#e8f5e9', '#e3f2fd', '#e1f5fe')  # Couleurs plus douces
        )
      )
  })
  
  # Tableau des résultats du test de Bartlett
  output$bartlett_result <- renderDT({
    req(adequacy_tests_data())
    
    # Extraire les résultats du test de Bartlett
    results <- adequacy_tests_data()
    
    # Créer un dataframe pour le tableau - nous nous assurons que tous les vecteurs ont la même longueur
    bartlett_df <- data.frame(
      Statistique = c("Chi-carré", "Degrés de liberté", "P-valeur"),
      Valeur = c(
        sprintf("%.2f", results$bartlett),
        as.character(results$bartlett_df),
        format.pval(results$bartlett_p, digits = 4)
      )
    )
    
    # Ajouter la colonne d'interprétation séparément pour s'assurer que les dimensions correspondent
    bartlett_df$Interprétation <- c(
      "Mesure de déviation par rapport à l'indépendance",
      paste0(results$n_vars * (results$n_vars - 1) / 2, " comparaisons de paires de variables"),
      ifelse(results$bartlett_p < 0.05, 
             "Test significatif - Rejette l'hypothèse nulle d'indépendance", 
             "Test non significatif - Ne peut pas rejeter l'indépendance")
    )
    
    # Formater le tableau
    datatable(bartlett_df, 
              options = list(dom = 't', ordering = FALSE, paging = FALSE),
              rownames = FALSE,
              caption = htmltools::tags$caption(
                style = 'caption-side: top; text-align: center; font-size: 120%; font-weight: bold;',
                'Test de sphéricité de Bartlett'
              )) %>%
      formatStyle(
        'Interprétation',
        target = 'row',
        backgroundColor = styleEqual(
          c("Test significatif - Rejette l'hypothèse nulle d'indépendance", 
            "Test non significatif - Ne peut pas rejeter l'indépendance"),
          c('#d4edda', '#f8d7da')
        )
      )
  })

  # Création du tableau des eigenvalues pour la matrice de corrélation
  output$eigenvalues_table_cormat <- renderDT({
    req(input$data_type == "cormat")
    req(cormat_data())
    
    # Récupérer la matrice de corrélation
    cor_mat <- cormat_data()
    
    # Calcul des valeurs propres
    eigenvalues <- eigen(cor_mat)$values
    
    # Calculer le pourcentage de variance expliquée
    variance_pct <- (eigenvalues / sum(eigenvalues)) * 100
    
    # Calculer le pourcentage cumulatif
    cumul_variance_pct <- cumsum(variance_pct)
    
    # Créer le dataframe pour le tableau
    eigenvalues_df <- data.frame(
      Facteur = 1:length(eigenvalues),
      Valeur_propre = eigenvalues,
      Variance_pct = variance_pct,
      Variance_cumul_pct = cumul_variance_pct
    )
    
    # Créer le datatable
    datatable(eigenvalues_df,
              colnames = c("Facteur", "Valeur propre", "% Variance", "% Variance cumulée"),
              options = list(
                pageLength = 10,
                dom = 'tip',
                ordering = TRUE,
                searching = FALSE
              ),
              rownames = FALSE) %>%
      formatRound(columns = c("Valeur_propre"), digits = 3) %>%
      formatRound(columns = c("Variance_pct", "Variance_cumul_pct"), digits = 2) %>%
      formatStyle(
        'Valeur_propre',
        backgroundColor = styleInterval(
          c(1),
          c('#fff3cd', '#e8f5e9')
        ),
        fontWeight = styleInterval(c(1), c('normal', 'bold'))
      )
  })
  
  # Création du tableau des eigenvalues pour les données brutes
  output$eigenvalues_table_raw <- renderDT({
    req(input$data_type == "raw")
    req(analysis_data())
    
    # Calculer la matrice de corrélation
    if (input$cor == "cor") {
      cor_mat <- cor(analysis_data(), use = "pairwise.complete.obs", method = "pearson")
    } else if (input$cor == "spearman") {
      cor_mat <- cor(analysis_data(), use = "pairwise.complete.obs", method = "spearman")
    } else {
      cor_mat <- cov(analysis_data(), use = "pairwise.complete.obs")
    }
    
    # Calcul des valeurs propres
    eigenvalues <- eigen(cor_mat)$values
    
    # Calculer le pourcentage de variance expliquée
    variance_pct <- (eigenvalues / sum(eigenvalues)) * 100
    
    # Calculer le pourcentage cumulatif
    cumul_variance_pct <- cumsum(variance_pct)
    
    # Créer le dataframe pour le tableau
    eigenvalues_df <- data.frame(
      Facteur = 1:length(eigenvalues),
      Valeur_propre = eigenvalues,
      Variance_pct = variance_pct,
      Variance_cumul_pct = cumul_variance_pct
    )
    
    # Créer le datatable
    datatable(eigenvalues_df,
              colnames = c("Facteur", "Valeur propre", "% Variance", "% Variance cumulée"),
              options = list(
                pageLength = 10,
                dom = 'tip',
                ordering = TRUE,
                searching = FALSE
              ),
              rownames = FALSE) %>%
      formatRound(columns = c("Valeur_propre"), digits = 3) %>%
      formatRound(columns = c("Variance_pct", "Variance_cumul_pct"), digits = 2) %>%
      formatStyle(
        'Valeur_propre',
        backgroundColor = styleInterval(
          c(1),
          c('#fff3cd', '#e8f5e9')
        ),
        fontWeight = styleInterval(c(1), c('normal', 'bold'))
      )
  })
  
  # Test MAP pour la matrice de corrélation
  map_result_reactive_cormat <- reactive({
    req(input$data_type == "cormat")
    req(cormat_data())
    
    # Créer un dataframe temporaire pour le test MAP
    # On doit simuler des données ayant cette matrice de corrélation
    
    # D'abord, vérifier si la matrice est définie positive
    eigen_values <- eigen(cormat_data())$values
    if (any(eigen_values <= 0)) {
      showNotification("La matrice de corrélation n'est pas définie positive. Test MAP impossible.", type = "error")
      return(NULL)
    }
    
    # Tenter d'exécuter le test MAP de Velicer
    tryCatch({
      result <- psych::VSS(cormat_data(), n = min(ncol(cormat_data())-1, 10), 
                           plot = FALSE, cor = cormat_data(), n.obs = input$n_obs)
      return(result)
    }, error = function(e) {
      showNotification(paste("Erreur lors de l'exécution du test MAP:", e$message), type = "error")
      return(NULL)
    })
  })
  
  output$map_test_plot_cormat <- renderPlot({
    req(map_result_reactive_cormat())
    result <- map_result_reactive_cormat()
    
    # Vérifier si le test a réussi
    if (!is.null(result$map)) {
      # Créer un dataframe pour le graphique
      map_values <- result$map
      map_data <- data.frame(
        Factors = 0:(length(map_values)),  # Inclure le facteur 0 pour l'indice de départ
        MAP_Criterion = c(1, map_values)   # Ajouter 1 comme première valeur
      )
      
      # Déterminer le nombre optimal de facteurs (le point le plus bas)
      opt_factors <- which.min(map_data$MAP_Criterion) - 1
      
      # Créer le graphique
      ggplot(map_data, aes(x = Factors, y = MAP_Criterion)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        geom_vline(xintercept = opt_factors, linetype = "dashed", color = "red") +
        scale_x_continuous(breaks = 0:max(map_data$Factors)) +
        labs(title = "Test MAP de Velicer",
             subtitle = paste("Le nombre optimal de facteurs est:", opt_factors),
             x = "Nombre de facteurs",
             y = "Critère MAP") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"),
              axis.title = element_text(face = "bold"))
    } else {
      # Afficher un message d'erreur
      plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
      text(1, 1, "Impossible de calculer le test MAP.\nEssayez avec une matrice définie positive.")
    }
  })
  
  output$map_test_result_cormat <- renderPrint({
    req(map_result_reactive_cormat())
    result <- map_result_reactive_cormat()
    
    # Vérifier si le test a réussi
    if (!is.null(result$map)) {
      # Ajouter 1 comme première valeur pour le critère MAP complet
      map_values <- c(1, result$map)
      
      # Déterminer le nombre optimal de facteurs (le point le plus bas)
      opt_factors <- which.min(map_values) - 1
      
      cat("Test MAP de Velicer\n")
      cat("------------------\n")
      cat("Le nombre optimal de facteurs selon le critère MAP est:", opt_factors, "\n")
      cat("\n")
      cat("Valeurs du critère MAP:\n")
      print(round(map_values, 4))
    } else {
      cat("Impossible de calculer le test MAP.\n")
      cat("Essayez avec une matrice de corrélation définie positive ou un nombre d'observations plus élevé.")
    }
  })
  
  # Fonction pour le test parallèle avec la matrice de corrélation
  parallel_results_auto_cormat <- eventReactive(
    c(input$n_datasets_auto_cormat, input$percentile_auto_cormat), 
    {
      req(input$data_type == "cormat")
      req(cormat_data())
      
      withProgress(message = 'Exécution du test parallèle...', {
        tryCatch({
          # Simuler des données à partir de la matrice de corrélation
          # pour l'analyse parallèle
          result <- fa.parallel(cormat_data(), 
                                n.obs = input$n_obs,  # Utiliser le nombre d'observations spécifié
                                fm = input$fm, 
                                n.iter = input$n_datasets_auto_cormat,
                                use = "pairwise.complete.obs",
                                plot = FALSE)  # Ne pas afficher le graphique ici
          
          return(result)
        }, error = function(e) {
          showNotification(paste("Erreur dans le test parallèle:", e$message), type = "error")
          return(NULL)
        })
      })
    }
  )
  
  # Graphique du test parallèle pour la matrice de corrélation
  output$parallel_plot_auto_cormat <- renderPlot({
    req(parallel_results_auto_cormat())
    
    result <- parallel_results_auto_cormat()
    
    # Extraire les valeurs propres observées
    observed <- result$fa.values
    
    # Calculer le percentile spécifié des valeurs simulées
    percentile_value <- input$percentile_auto_cormat / 100
    simulated_percentiles <- apply(result$values, 2, function(x) quantile(x, probs = percentile_value))
    
    # S'assurer que nous avons assez de valeurs simulées
    n_factors <- length(observed)
    if (length(simulated_percentiles) < n_factors) {
      # Compléter avec la dernière valeur si nécessaire
      simulated_percentiles <- c(simulated_percentiles, 
                                 rep(simulated_percentiles[length(simulated_percentiles)], 
                                     n_factors - length(simulated_percentiles)))
    } else {
      # Limiter au nombre de facteurs observés
      simulated_percentiles <- simulated_percentiles[1:n_factors]
    }
    
    # Créer un dataframe pour le graphique
    df <- data.frame(
      Factor = c(1:n_factors, 1:n_factors),
      Eigenvalue = c(observed, simulated_percentiles),
      Type = c(rep("Données observées", n_factors), 
               rep(paste0("Données simulées (P", input$percentile_auto_cormat, ")"), n_factors))
    )
    
    # Déterminer le nombre de facteurs à retenir (où valeurs observées > percentile simulé)
    greater_than_random <- observed > simulated_percentiles
    factors_to_retain <- sum(greater_than_random)
    
    # Nous devons créer le mapping des couleurs avant de l'utiliser
    color_values <- c("blue", "red")
    names(color_values) <- c("Données observées", paste0("Données simulées (P", input$percentile_auto_cormat, ")"))
    
    # Créer le graphique avec ggplot2
    p <- ggplot(df, aes(x = Factor, y = Eigenvalue, group = Type, color = Type)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      scale_color_manual(values = color_values) +
      labs(title = "Analyse parallèle",
           subtitle = paste("Comparaison des valeurs propres observées avec le", 
                            input$percentile_auto_cormat, "e percentile de", 
                            input$n_datasets_auto_cormat, "jeux de données simulés"),
           x = "Facteur", y = "Valeur propre") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_blank())
    
    # Ajouter la ligne verticale et l'annotation pour le nombre de facteurs à retenir
    if (factors_to_retain > 0) {
      p <- p + 
        geom_vline(xintercept = factors_to_retain + 0.5, linetype = "dashed", color = "darkgreen") +
        annotate("text", x = factors_to_retain + 1, y = max(df$Eigenvalue) * 0.9, 
                 label = paste("Retenir", factors_to_retain, "facteurs"), 
                 color = "darkgreen", fontface = "bold")
    }
    
    return(p)
  })
  
  # Résultats textuels du test parallèle pour la matrice de corrélation
  output$parallel_summary_auto_cormat <- renderPrint({
    req(parallel_results_auto_cormat())
    
    result <- parallel_results_auto_cormat()
    
    # Extraire les valeurs propres observées
    observed <- result$fa.values
    
    # Calculer le percentile spécifié des valeurs simulées
    percentile_value <- input$percentile_auto_cormat / 100
    simulated_percentiles <- apply(result$values, 2, function(x) quantile(x, probs = percentile_value))
    
    # S'assurer que nous avons assez de valeurs simulées
    n_factors <- length(observed)
    if (length(simulated_percentiles) < n_factors) {
      # Compléter avec la dernière valeur si nécessaire
      simulated_percentiles <- c(simulated_percentiles, 
                                 rep(simulated_percentiles[length(simulated_percentiles)], 
                                     n_factors - length(simulated_percentiles)))
    } else {
      # Limiter au nombre de facteurs observés
      simulated_percentiles <- simulated_percentiles[1:n_factors]
    }
    
    # Déterminer le nombre de facteurs à retenir (où valeurs observées > percentile simulé)
    greater_than_random <- observed > simulated_percentiles
    factors_to_retain <- sum(greater_than_random)
    
    cat("Résultats de l'analyse parallèle\n")
    cat("------------------------------\n")
    cat("Nombre de facteurs à retenir selon l'analyse parallèle:", factors_to_retain, "\n")
    cat("(facteurs dont la valeur propre observée > percentile", input$percentile_auto_cormat, "des valeurs simulées)\n\n")
    
    # Notons aussi la valeur renvoyée par fa.parallel() qui peut être différente
    if (result$nfact != factors_to_retain) {
      cat("Note: La fonction fa.parallel() suggère", result$nfact, "facteurs, mais l'examen des valeurs\n")
      cat("montre que", factors_to_retain, "facteurs ont des valeurs propres supérieures au percentile", 
          input$percentile_auto_cormat, "des données simulées.\n\n")
    }
    
    cat("Comparaison des valeurs propres:\n")
    
    # S'assurer que les dimensions correspondent
    n_to_compare <- min(length(observed), length(simulated_percentiles), 10)
    
    # Créer d'abord les noms de colonnes
    col_names <- c("Facteur", "Valeur propre observée", 
                   paste0("Valeur simulée (P", input$percentile_auto_cormat, ")"), 
                   "Différence", "À retenir")
    
    # Créer le dataframe avec des noms de colonnes standard
    comparison <- data.frame(
      "Facteur" = 1:n_to_compare,
      "Valeur_propre_observee" = observed[1:n_to_compare],
      "Valeur_simulee" = simulated_percentiles[1:n_to_compare],
      "Difference" = observed[1:n_to_compare] - simulated_percentiles[1:n_to_compare],
      "A_retenir" = ifelse(observed[1:n_to_compare] > simulated_percentiles[1:n_to_compare], "Oui", "Non")
    )
    
    # Renommer les colonnes après la création
    names(comparison) <- col_names
    
    print(comparison)
  })
  

  
  # Graphique EGA avec options contrôlables pour matrice de corrélation
  output$ega_plot_cormat <- renderPlot({
    req(ega_results())
    
    # Afficher le graphique EGA avec les options choisies par l'utilisateur
    plot(ega_results(), 
         title = "Analyse Graphique Exploratoire (EGA)",
         vsize = input$ega_node_size_cormat,   # taille des nœuds contrôlable
         label.size = input$ega_label_size_cormat,  # taille des étiquettes contrôlable
         edge.alpha = input$ega_edge_alpha_cormat,  # transparence des arêtes contrôlable 
         layout = input$ega_layout_cormat  # layout contrôlable
    )
  })
  
  # Résumé EGA pour matrice de corrélation
  output$ega_summary_cormat <- renderUI({
    req(ega_results())
    
    result <- ega_results()
    
    # Calculer les statistiques pour chaque dimension
    dim_table <- table(result$wc)
    
    # Créer un dataframe pour les variables par dimension
    var_dims <- data.frame(
      Variable = names(result$wc),
      Dimension = result$wc
    )
    
    # Trier par dimension, puis par nom de variable
    var_dims <- var_dims[order(var_dims$Dimension, var_dims$Variable), ]
    
    # Créer une liste de variables par dimension pour l'affichage en accordion
    dim_vars_list <- split(var_dims$Variable, var_dims$Dimension)
    
    # Couleurs pour les dimensions (palette RColorBrewer Set1)
    dim_colors <- RColorBrewer::brewer.pal(min(9, result$n.dim), "Set1")
    if(result$n.dim > 9) {
      # Étendre la palette si nécessaire
      dim_colors <- colorRampPalette(dim_colors)(result$n.dim)
    }
    
    # Créer l'interface utilisateur
    tagList(
      div(
        class = "mb-4",
        h3("Résultats de l'Analyse Graphique Exploratoire (EGA)"),
        
        # Informations principales dans une card
        card(
          card_header("Informations générales"),
          div(
            class = "p-3",
            layout_columns(
              col_widths = c(6, 6),
              div(
                h4("Nombre de dimensions détectées"),
                div(
                  class = "display-4 text-center mb-3",
                  result$n.dim
                ),
                
                hr(),
                
                # Ajouter l'indice TEFI
                h5("Indice d'ajustement TEFI"),
                div(
                  class = "text-center",
                  div(
                    style = paste0(
                      "display: inline-block; padding: 8px 16px; border-radius: 5px; ",
                      "background-color: ", 
                      ifelse(result$TEFI < 0.3, "#d4edda", 
                             ifelse(result$TEFI < 0.5, "#fff3cd", "#f8d7da")),
                      "; font-weight: bold; font-size: 18px;"
                    ),
                    round(result$TEFI, 4)
                  ),
                  div(
                    class = "small text-muted mt-2",
                    "Indice de qualité d'ajustement (valeurs plus basses indiquent un meilleur ajustement)"
                  )
                )
              ),
              
              # Graphique de répartition des variables par dimension
              div(
                h4("Répartition des variables par dimension"),
                plotOutput("ega_dimension_distribution_cormat", height = "200px")
              )
            )
          )
        ),
        
        # Tableau détaillé des variables par dimension
        card(
          card_header("Distribution des variables par dimension"),
          div(
            class = "p-3",
            tabsetPanel(
              # Onglet 1: Tableau des variables
              tabPanel("Tableau des variables",
                       div(
                         class = "mt-3",
                         DTOutput("ega_dimension_table_cormat")
                       )
              ),
              
              # Onglet 2: Visualisation par dimension (accordion)
              tabPanel("Visualisation par dimension",
                       div(
                         class = "mt-3",
                         accordion(
                           id = "dim_accordion_cormat",
                           multiple = TRUE,
                           # Créer un panneau d'accordéon pour chaque dimension
                           lapply(1:result$n.dim, function(dim) {
                             accordion_panel(
                               title = paste0("Dimension ", dim, " (", length(dim_vars_list[[as.character(dim)]]), " variables)"),
                               div(
                                 style = paste0("border-left: 5px solid ", dim_colors[dim], "; padding-left: 15px;"),
                                 # Afficher les variables de cette dimension avec des puces
                                 lapply(dim_vars_list[[as.character(dim)]], function(var) {
                                   div(
                                     class = "mb-1",
                                     icon("circle"), " ", var
                                   )
                                 })
                               )
                             )
                           })
                         )
                       )
              ),
              
              # Onglet 3: Statistiques des dimensions
              tabPanel("Statistiques des dimensions",
                       div(
                         class = "mt-3",
                         h5("Nombre de variables par dimension"),
                         layout_columns(
                           col_widths = c(5, 7),
                           datatable(
                             data.frame(
                               Dimension = 1:length(dim_table),
                               `Nombre de variables` = as.numeric(dim_table),
                               `Pourcentage` = paste0(round(100 * as.numeric(dim_table) / sum(dim_table), 1), "%")
                             ),
                             options = list(dom = 't', ordering = FALSE),
                             rownames = FALSE
                           ),
                           plotOutput("ega_dim_barplot_cormat", height = "200px")
                         )
                       )
              )
            )
          )
        ),
        
        # Section supplémentaire pour l'interprétation
        card(
          card_header("Interprétation des résultats"),
          div(
            class = "p-3",
            p("L'analyse EGA a identifié ", tags$strong(result$n.dim), " dimensions distinctes dans les données:"),
            
            div(
              class = "mb-3 p-3 bg-light rounded",
              h5("Signification des dimensions"),
              p("Les dimensions identifiées par EGA représentent des clusters cohérents de variables, reflétant potentiellement des construits latents dans les données. Elles peuvent être comparées aux facteurs identifiés par l'analyse factorielle traditionnelle."),
              
              h5("Comment interpréter l'indice TEFI"),
              p("L'indice TEFI (Total Entropy Fit Index) mesure la qualité d'ajustement du modèle EGA:"),
              tags$ul(
                tags$li("Valeurs plus basses indiquent un meilleur ajustement"),
                tags$li("Généralement interprété en comparant différents modèles"),
                tags$li("Pas de seuil universel, mais les modèles avec TEFI < 0.3 sont souvent considérés comme bien ajustés")
              )
            )
          )
        )
      )
    )
  })
  
  # Ajouter ces nouveaux outputs pour les visualisations
  output$ega_dimension_distribution_cormat <- renderPlot({
    req(ega_results())
    
    # Créer un dataframe pour le graphique
    dim_counts <- as.data.frame(table(ega_results()$wc))
    names(dim_counts) <- c("Dimension", "Count")
    dim_counts$Dimension <- as.numeric(as.character(dim_counts$Dimension))
    
    # Couleurs pour les dimensions (palette RColorBrewer Set1)
    dim_colors <- RColorBrewer::brewer.pal(min(9, ega_results()$n.dim), "Set1")
    if (ega_results()$n.dim > 9) {
      dim_colors <- colorRampPalette(dim_colors)(ega_results()$n.dim)
    }
    
    # Créer un graphique en secteurs
    ggplot(dim_counts, aes(x = "", y = Count, fill = factor(Dimension))) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = dim_colors) +
      labs(fill = "Dimension") +
      theme_void() +
      theme(legend.position = "right")
  })
  
  output$ega_dimension_table_cormat <- renderDT({
    req(ega_results())
    
    # Créer un dataframe pour le tableau
    var_dims <- data.frame(
      Variable = names(ega_results()$wc),
      Dimension = ega_results()$wc
    )
    
    # Trier par dimension, puis par nom de variable
    var_dims <- var_dims[order(var_dims$Dimension, var_dims$Variable), ]
    
    # Créer un tableau interactif
    datatable(
      var_dims,
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50),
        order = list(list(1, 'asc')),
        dom = 'lfrtip'
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'Dimension',
        backgroundColor = styleEqual(
          1:ega_results()$n.dim,
          RColorBrewer::brewer.pal(min(9, ega_results()$n.dim), "Pastel1") %>%
            {if(ega_results()$n.dim > 9) colorRampPalette(.)(ega_results()$n.dim) else .}
        )
      )
  })
  
  output$ega_dim_barplot_cormat <- renderPlot({
    req(ega_results())
    
    # Créer un dataframe pour le graphique
    dim_counts <- as.data.frame(table(ega_results()$wc))
    names(dim_counts) <- c("Dimension", "Count")
    dim_counts$Dimension <- as.numeric(as.character(dim_counts$Dimension))
    
    # Couleurs pour les dimensions (palette RColorBrewer Set1)
    dim_colors <- RColorBrewer::brewer.pal(min(9, ega_results()$n.dim), "Set1")
    if (ega_results()$n.dim > 9) {
      dim_colors <- colorRampPalette(dim_colors)(ega_results()$n.dim)
    }
    
    # Créer un graphique en barres
    ggplot(dim_counts, aes(x = factor(Dimension), y = Count, fill = factor(Dimension))) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_manual(values = dim_colors) +
      labs(x = "Dimension", y = "Nombre de variables") +
      theme_minimal() +
      theme(legend.position = "none",
            panel.grid.minor = element_blank(),
            axis.title = element_text(face = "bold"))
  })
  
  # Matrice de réseau EGA pour matrice de corrélation
  output$ega_network_matrix_cormat <- renderDT({
    req(ega_results())
    
    # Obtenir la matrice du réseau
    network_matrix <- ega_results()$network
    
    # Convertir en format triangulaire inférieur (sans diagonale)
    network_lower <- network_matrix
    network_lower[upper.tri(network_lower, diag = TRUE)] <- NA
    
    # Convertir en dataframe pour l'affichage
    network_df <- as.data.frame(network_lower)
    
    # Créer un datatable avec la matrice triangulaire inférieure
    datatable(network_df, 
              options = list(
                scrollX = TRUE, 
                pageLength = 20,
                searching = FALSE,  # Désactiver la recherche
                lengthChange = FALSE  # Désactiver le changement de nombre d'éléments par page
              )) %>%
      formatRound(columns = 1:ncol(network_df), digits = 4) %>%
      formatStyle(
        columns = 1:ncol(network_df),
        backgroundColor = styleInterval(
          c(-0.5, -0.2, 0.2, 0.5), 
          c("#d73027", "#f46d43", "#ffffff", "#74add1", "#4575b4")
        )
      )
  })
  
  # UI pour afficher les network loads pour matrice de corrélation
  output$netloads_ui_cormat <- renderUI({
    req(ega_results())
    
    # Si les network loads sont calculés, les afficher
    if (!is.null(netloads_result())) {
      loads <- netloads_result()
      
      tagList(
        hr(),
        h4("Network Loads"),
        tabsetPanel(
          tabPanel("Loads standardisés", 
                   DTOutput("netloads_std_table_cormat")
          ),
          tabPanel("Loads non-standardisés", 
                   DTOutput("netloads_unstd_table_cormat")
          ),
          tabPanel("Loads pivotés", 
                   if(!is.null(loads$rotated)) {
                     tagList(
                       h5("Charges pivotées"),
                       DTOutput("netloads_rotated_table_cormat"),
                       h5("Corrélations entre dimensions"),
                       DTOutput("netloads_phi_table_cormat")
                     )
                   } else {
                     p("Aucune rotation n'a été appliquée. Choisissez une méthode de rotation dans les paramètres pour afficher les charges pivotées.")
                   }
          )
        )
      )
    } else {
      # Message d'attente au cas où le calcul prend du temps
      tagList(
        div(
          style = "text-align: center; padding: 30px;",
          h4("Calcul des Network Loads en cours..."),
          p("Les network loads sont calculés automatiquement. Veuillez patienter quelques instants.")
        )
      )
    }
  })

  
  # Rendu des tableaux netloads pour matrice de corrélation
  output$netloads_std_table_cormat <- renderDT({
    req(netloads_result())
    datatable(netloads_result()$std, 
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatRound(columns = 1:ncol(netloads_result()$std), digits = 3)
  })
  
  output$netloads_unstd_table_cormat <- renderDT({
    req(netloads_result())
    datatable(netloads_result()$unstd, 
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatRound(columns = 1:ncol(netloads_result()$unstd), digits = 3)
  })
  
  output$netloads_rotated_table_cormat <- renderDT({
    req(netloads_result())
    req(netloads_result()$rotated)
    
    datatable(netloads_result()$rotated$loadings, 
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatRound(columns = 1:ncol(netloads_result()$rotated$loadings), digits = 3)
  })
  
  output$netloads_phi_table_cormat <- renderDT({
    req(netloads_result())
    req(netloads_result()$rotated)
    
    datatable(netloads_result()$rotated$Phi, 
              options = list(scrollX = TRUE, pageLength = 10)) %>%
      formatRound(columns = 1:ncol(netloads_result()$rotated$Phi), digits = 3)
  })
  
  # Handler pour télécharger les network loads pour matrice de corrélation
  output$download_netloads_cormat <- downloadHandler(
    filename = function() {
      paste("network_loads_cormat_", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      req(netloads_result())
      
      # Créer un workbook Excel
      wb <- openxlsx::createWorkbook()
      
      # Ajouter les feuilles
      openxlsx::addWorksheet(wb, "Loads standardisés")
      openxlsx::writeData(wb, "Loads standardisés", netloads_result()$std)
      
      openxlsx::addWorksheet(wb, "Loads non-standardisés")
      openxlsx::writeData(wb, "Loads non-standardisés", netloads_result()$unstd)
      
      if(!is.null(netloads_result()$rotated)) {
        openxlsx::addWorksheet(wb, "Loads pivotés")
        openxlsx::writeData(wb, "Loads pivotés", netloads_result()$rotated$loadings)
        
        openxlsx::addWorksheet(wb, "Corrélations entre dimensions")
        openxlsx::writeData(wb, "Corrélations entre dimensions", netloads_result()$rotated$Phi)
      }
      
      # Sauvegarder le workbook
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

}

# Lancer l'application
shinyApp(ui = ui, server = server)