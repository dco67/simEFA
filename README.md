# simEFA — L'analyse factorielle exploratoire

**Simuler une structure factorielle connue, puis voir si l'analyse factorielle exploratoire la retrouve : extraction, nombre de facteurs, rotation et interprétation.**

simEFA est une application pédagogique **R/Shiny**, en français, consacrée à l'analyse factorielle exploratoire (AFE). Son laboratoire de simulation crée des données dont on contrôle la structure : nombre de facteurs, variables par facteur, force des saturations principales et secondaires, items ambigus. On peut aussi charger un jeu de données classique ou importer son propre fichier.

- [Utiliser l'application en ligne](https://shiny.dcstats.org/simEFA/)
- [Consulter le manuel pédagogique](https://dcstats.org/simEFA/)

## Objectifs pédagogiques

- vérifier si les données se prêtent à une AFE (KMO, test de Bartlett);
- choisir le nombre de facteurs avec plusieurs critères, sans s'en remettre à un seul;
- comprendre l'extraction des facteurs, les saturations et les communautés;
- comparer les rotations orthogonales (Varimax) et obliques (Oblimin…);
- nommer les facteurs et évaluer leur fiabilité;
- comparer la solution obtenue à la structure qui a réellement généré les données.

## Fonctionnalités

| Onglet | Contenu |
| --- | --- |
| Données | Laboratoire de simulation (situations prêtes à explorer), jeux classiques, import (CSV, Excel, SPSS…), inversion d'items. |
| Population simulée | Ce que la simulation a réellement généré. |
| Statistiques descriptives | Résumé, distributions, matrice de corrélation. |
| Analyse factorielle | Tests préliminaires, résultats principaux, comparaison des modèles, rotation dynamique, fiabilité des facteurs. |
| Interprétation | Lecture guidée des saturations et des facteurs. |

L'AFE est réalisée sur les corrélations de Pearson; l'application signale quand des items Likert très asymétriques appelleraient plutôt une analyse sur corrélations polychoriques.

## Exécuter l'application localement

```r
install.packages(c("shiny", "bslib", "shinyjs", "DT", "ggplot2", "plotly", "corrplot", "psych",
                   "lavaan", "haven", "readxl", "writexl"))
shiny::runApp("simEFA")
```

## Licence

GPL-3.0. Voir le fichier `LICENSE`.

## Auteur

Daniel Coulombe, Ph.D. — <https://dcstats.org>
