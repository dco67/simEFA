library(shiny)
library(bslib)
library(EFAtools)
library(dplyr)
library(ggplot2)
library(DT)

ui <- page_sidebar(
  title = "Analyse Factorielle Exploratoire avec EFAtools",
  sidebar = sidebar(
    fileInput("data_file", "Importer les données (CSV ou Excel)",
              accept = c(".csv", ".xls", ".xlsx")),
    selectInput("method", "Méthode d'extraction",
                choices = c("PAF", "ML", "ULS"),
                selected = "PAF"),
    numericInput("n_factors", "Nombre de facteurs", 
                 value = 2, min = 1, max = 20),
    selectInput("rotation", "Méthode de rotation",
                choices = c("none", "varimax", "promax", "oblimin", "quartimax", 
                            "equamax", "geominT", "bentlerT", "bifactorT", 
                            "quartimin", "simplimax", "bentlerQ", "geominQ", "bifactorQ"),
                selected = "none"),
    selectInput("type", "Type d'implémentation",
                choices = c("EFAtools", "psych", "SPSS", "none"),
                selected = "EFAtools"),
    selectInput("init_comm", "Méthode pour estimer les communalités initiales (pour PAF)",
                choices = c("NA" = NA, "smc", "mac", "unity"),
                selected = NA),
    numericInput("max_iter", "Nombre maximum d'itérations",
                 value = NA, min = 1, max = 1000),
    numericInput("criterion", "Critère de convergence",
                 value = NA, min = 0.000001, max = 0.1, step = 0.00001),
    selectInput("criterion_type", "Type de critère de convergence",
                choices = c("NA" = NA, "max_individual", "sum"),
                selected = NA),
    checkboxInput("normalize", "Normalisation de Kaiser", value = TRUE),
    selectInput("order_type", "Méthode de tri des facteurs",
                choices = c("NA" = NA, "eigen", "ss_factors"),
                selected = NA),
    actionButton("run_analysis", "Exécuter l'analyse", class = "btn-primary")
  ),
  
  layout_columns(
    card(
      card_header("Aperçu des données"),
      card_body(DTOutput("data_preview"))
    ),
    card(
      card_header("Résultats de l'analyse factorielle"),
      card_body(
        tabsetPanel(
          tabPanel("Variance expliquée", verbatimTextOutput("variance_explained")),
          tabPanel("Loadings", verbatimTextOutput("loadings")),
          tabPanel("Communalités", verbatimTextOutput("communalities")),
          tabPanel("Structure", verbatimTextOutput("structure")),
          tabPanel("Correlation des facteurs", verbatimTextOutput("factor_cors")),
          tabPanel("Scree plot", plotOutput("scree_plot")),
          tabPanel("Biplots", plotOutput("biplot")),
          tabPanel("Résumé complet", verbatimTextOutput("summary"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Réactif pour les données
  data <- reactive({
    req(input$data_file)
    
    ext <- tools::file_ext(input$data_file$name)
    
    if (ext == "csv") {
      df <- read.csv(input$data_file$datapath)
    } else if (ext %in% c("xls", "xlsx")) {
      df <- readxl::read_excel(input$data_file$datapath)
    } else {
      stop("Format de fichier non pris en charge")
    }
    
    return(df)
  })
  
  # Aperçu des données
  output$data_preview <- renderDT({
    req(data())
    datatable(data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  # Création d'un réactif pour les résultats EFA
  efa_results <- eventReactive(input$run_analysis, {
    req(data())
    
    # Gérer les valeurs NA
    max_iter_val <- if(is.na(input$max_iter)) NA else input$max_iter
    criterion_val <- if(is.na(input$criterion)) NA else input$criterion
    criterion_type_val <- if(input$criterion_type == "NA") NA else input$criterion_type
    init_comm_val <- if(input$init_comm == "NA") NA else input$init_comm
    order_type_val <- if(input$order_type == "NA") NA else input$order_type
    
    # Exécuter l'analyse EFA
    tryCatch({
      EFA(x = data(), 
          n_factors = input$n_factors,
          method = input$method,
          rotation = input$rotation,
          type = input$type,
          max_iter = max_iter_val,
          init_comm = init_comm_val,
          criterion = criterion_val,
          criterion_type = criterion_type_val,
          normalize = input$normalize,
          order_type = order_type_val)
    }, error = function(e) {
      showNotification(paste("Erreur dans l'analyse:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Afficher la variance expliquée
  output$variance_explained <- renderPrint({
    req(efa_results())
    print(efa_results()$vars_accounted)
  })
  
  # Afficher les loadings
  output$loadings <- renderPrint({
    req(efa_results())
    print(efa_results()$rot_loadings)
  })
  
  # Afficher les communalités
  output$communalities <- renderPrint({
    req(efa_results())
    comm_df <- data.frame(
      Variable = rownames(efa_results()$h2),
      Communalité = efa_results()$h2
    )
    print(comm_df)
  })
  
  # Afficher la structure factorielle (si disponible)
  output$structure <- renderPrint({
    req(efa_results())
    if (!is.null(efa_results()$Phi) && input$rotation != "none") {
      structure <- efa_results()$rot_loadings %*% efa_results()$Phi
      colnames(structure) <- colnames(efa_results()$rot_loadings)
      rownames(structure) <- rownames(efa_results()$rot_loadings)
      print(structure)
    } else {
      print("Structure factorielle non disponible pour cette rotation")
    }
  })
  
  # Afficher les corrélations entre facteurs (si disponible)
  output$factor_cors <- renderPrint({
    req(efa_results())
    if (!is.null(efa_results()$Phi) && input$rotation != "none") {
      print(efa_results()$Phi)
    } else {
      print("Corrélations entre facteurs non disponibles pour cette rotation")
    }
  })
  
  # Scree plot
  output$scree_plot <- renderPlot({
    req(efa_results())
    
    # Extraire les valeurs propres
    eigen_values <- efa_results()$eigenvalues
    
    # Créer un data frame pour ggplot
    eigen_df <- data.frame(
      Factor = 1:length(eigen_values),
      Eigenvalue = eigen_values
    )
    
    # Créer le scree plot
    ggplot(eigen_df, aes(x = Factor, y = Eigenvalue)) +
      geom_line() +
      geom_point(size = 3) +
      theme_minimal() +
      labs(title = "Scree Plot", x = "Facteur", y = "Valeur propre") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = 1:length(eigen_values))
  })
  
  # Biplot des facteurs
  output$biplot <- renderPlot({
    req(efa_results())
    
    # Récupérer les loadings
    loadings <- efa_results()$rot_loadings
    
    # Sélectionner les deux premiers facteurs pour le plot
    if (ncol(loadings) >= 2) {
      loadings_df <- data.frame(
        Variable = rownames(loadings),
        Factor1 = loadings[, 1],
        Factor2 = loadings[, 2]
      )
      
      # Créer le biplot
      ggplot(loadings_df, aes(x = Factor1, y = Factor2, label = Variable)) +
        geom_point(color = "blue") +
        geom_text(vjust = -0.5, hjust = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_fixed() +  # Échelle égale pour les deux axes
        theme_minimal() +
        labs(title = "Biplot des deux premiers facteurs",
             x = paste("Facteur 1 (", round(efa_results()$vars_accounted$Proportion[1] * 100, 1), "%)", sep = ""),
             y = paste("Facteur 2 (", round(efa_results()$vars_accounted$Proportion[2] * 100, 1), "%)", sep = ""))
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Au moins 2 facteurs sont nécessaires pour un biplot") +
        theme_void()
    }
  })
  
  # Résumé complet
  output$summary <- renderPrint({
    req(efa_results())
    print(efa_results())
  })
}

shinyApp(ui, server)
