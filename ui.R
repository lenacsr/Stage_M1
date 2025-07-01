#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyWidgets)

library(readr)
library(deconica)
library(tibble)
library(FactoMineR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(fgsea)
library(knitr)
library(tidyr)
library(tidyverse)

library(factoextra)
library(cluster)
library(NbClust)

library(shinyjs)
library(shinycssloaders)
library(shiny)
library(plotly)
library(stringr)

library(kableExtra)
library(bslib)
library(DT)
library(shinycustomloader)

all_basenames_38 <- unique(score_grouped_XBP1_RIDD$BaseName)
all_basenames_pur <- unique(score_grouped_XBP1_RIDD_pur$BaseName)
all_basenames_SEA <- unique(score_grouped_CancerSEA$BaseName)
all_basenames_Hallmarks <- unique(score_grouped_Hallmark$BaseName)


calculator_fill <- HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-calculator-fill" viewBox="0 0 16 16">
                           <path d="M2 2a2 2 0 0 1 2-2h8a2 2 0 0 1 2 2v12a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2zm2 .5v2a.5.5 0 0 0 .5.5h7a.5.5 0 0 0 .5-.5v-2a.5.5 0 0 0-.5-.5h-7a.5.5 0 0 0-.5.5m0 4v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5h-1a.5.5 0 0 0-.5.5M4.5 9a.5.5 0 0 0-.5.5v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5zM4 12.5v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5h-1a.5.5 0 0 0-.5.5M7.5 6a.5.5 0 0 0-.5.5v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5zM7 9.5v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5h-1a.5.5 0 0 0-.5.5m.5 2.5a.5.5 0 0 0-.5.5v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5zM10 6.5v1a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-1a.5.5 0 0 0-.5-.5h-1a.5.5 0 0 0-.5.5m.5 2.5a.5.5 0 0 0-.5.5v4a.5.5 0 0 0 .5.5h1a.5.5 0 0 0 .5-.5v-4a.5.5 0 0 0-.5-.5z"/>
                           </svg>')
indice_coude <- HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-graph-down-arrow" viewBox="0 0 16 16">
  <path fill-rule="evenodd" d="M0 0h1v15h15v1H0zm10 11.5a.5.5 0 0 0 .5.5h4a.5.5 0 0 0 .5-.5v-4a.5.5 0 0 0-1 0v2.6l-3.613-4.417a.5.5 0 0 0-.74-.037L7.06 8.233 3.404 3.206a.5.5 0 0 0-.808.588l4 5.5a.5.5 0 0 0 .758.06l2.609-2.61L13.445 11H10.5a.5.5 0 0 0-.5.5"/>
</svg>')

indice_silhouette <- HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-person-standing" viewBox="0 0 16 16">
  <path d="M8 3a1.5 1.5 0 1 0 0-3 1.5 1.5 0 0 0 0 3M6 6.75v8.5a.75.75 0 0 0 1.5 0V10.5a.5.5 0 0 1 1 0v4.75a.75.75 0 0 0 1.5 0v-8.5a.25.25 0 1 1 .5 0v2.5a.75.75 0 0 0 1.5 0V6.5a3 3 0 0 0-3-3H7a3 3 0 0 0-3 3v2.75a.75.75 0 0 0 1.5 0v-2.5a.25.25 0 0 1 .5 0"/>
</svg>')

ui <-
  fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-image: url('bg.jpg');
        background-size: cover;
        background-position: center;
        background-repeat: no-repeat;
        background-attachment: fixed;
        color: white;
      }
      
      .plot-container {
        position: relative;
        min-height: 400px;
      }
      
      /* C'est pour avoir le summary de TCGA en blanc */
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_filter input,
      .dataTables_wrapper .dataTables_length select,
      .dataTables_wrapper .dataTables_paginate .paginate_button {
        color: white !important;
      }
      /* Style pour le label Search */
      .dataTables_wrapper .dataTables_filter label {
        color: white !important;
      }
      .dataTables_wrapper table.dataTable tr th,
      .dataTables_wrapper table.dataTable tr td {
        color: white !important;
      }
      /* Style pour les boutons Previous et Next */
      .dataTables_wrapper .dataTables_paginate .paginate_button.previous{
        color: white !important;
      }
      .dataTables_wrapper .dataTables_paginate .paginate_button.next {
        color: white !important;
      }
       /* Style pour le texte Show et entries */
      .dataTables_wrapper .dataTables_length label {
        color: white !important;
      }
      /* Style pour les autres éléments */
      .dataTables_wrapper .dataTables_info {
        color: white !important;
      }

      /* Style pour les panels */
      .well {
        background-color: rgba(248, 248, 248, 0.9);
      }
      
      .kable caption {
        color: white !important; /* Changer la couleur du texte de la légende */
      }
      
      /* Change tab background color */
      .nav-tabs > li > a {
        
        color: white;
      }

      /* Change active tab background color */
      .nav-tabs > li.active > a {
        background-color: orchid !important;
        color: white !important;
      }
      
      /*avoir une value_box blanche*/
      .custom-value-box {
      border: 2px solid white !important;
      margin: 5px !important;
      }
      
      
      .custom-loader {
    width: 80px;
    max-width: 100px;
    height: auto;
    display: block;
    margin-left: auto;
    margin-right: auto;
      }
  
    

    
    "),
              
    )),
  div(class = "main-content",
  titlePanel("Etude autour de la déconvolution de glioblastomes"),
  
  
  tabsetPanel(
    
    # Welcome Tab
    tabPanel("Accueil",
             
             h2("Introduction et contexte"),
             fluidPage(
               
             p("Le ",strong("glioblastome"), "est la tumeur du cerveau la plus fréquente et agressive. Une des raisons de sa létalité est l'hétérogénéité des types cellulaires qu'elle contient. La tumeur contient donc des cellules différentes que l'on peut mettre dans des 'catégories' différentes selon leur expression transcriptomique, c'est-à-dire l'ensemble des ARN que la cellule produit. Les différents",strong("états cellulaires"),"  du glioblastome se nomment:"),
            tags$ul(
              tags$li("AC-like"),
              tags$li("MES-like"),
              tags$li("OPC-like"),
              tags$li("NPC-like")),
            p("Certains chercheurs pensent aussi avoir trouvé des nouveaux états cellulaires : ODC-like et stem-like."),
            p("A partir du", strong("transcriptome"), "(l'ensemble des ARN produit) d'une cellule, on peut aussi retrouver des ",strong("voies de signalisations"),", qui désignent une suite de réactions chimiques où un groupe de molécule d'une cellule travaille ensemble pour contrôler des fonctions d'une cellule comme l'apoptose (mort programmé de la cellule) ou la division de la cellule par exemple. Deux voies de signalisations qui vont particulièrement nous intéresser sont", strong("XBP1s")," et ",strong("RIDD"),". Ces voies de signalisation sont activées en cas de stress du",strong ("réticulum endoplasmique")," (un organite de la cellule). Le réticulum endoplasmique peut être stressé en cas de cancer, ou d'une carence en nutriment, ou bien même un choc thermique dans ce cas, il active ",strong("IRE1"), "qui est la protéine transmembranaire qui lance ces voies de signalisation. Ainsi, lorsque le réticulum endoplasmique est stressé elle active la voie de signalisation XBP1s qui permet d'augmenter la production de protéine utile dans cette période de stress. Elle peut aussi activer la voie de signalisation RIDD qui pour 'soulager' la cellule stressée va détruire les ARNs afin de réduire la production de protéine et de lui donner moins de travail en période de stress."),
            p("Ainsi pour mieux comprendre et traiter les glioblastomes, il est important de connaître dans quel 'état' est chaque cellule et quelle voie de signalisation elle active. Cependant, lorsque l'on prélève des échantillons de tumeurs et qu'on récupère leurs transcriptomes, on n'a pas le transcriptome de chaque cellule. Mais on a l'ensemble des transcriptomes des cellules de l'échantillon mélangés, ne formant qu'un seul gros transcriptome."),
            p("Nous allons donc tenter de", strong("déconvoluer")," dans nos échantillons TCGA (avec une expression génétique pour chaque échantillon), c'est-à-dire retrouver la quantité des différents états cellulaires et voie de signalisation au sein d'un échantillon."),
            p("Un objectif complémentaire de ce projet est de déterminer la qualité des signatures de XBP1s et de RIDD, c'est-à-dire déterminer si la signature représente bien la voie de signalisation. Ici, quand on parle de signature, on parle de l'ensemble des gènes qui sont produit spécifiquement lorsque la voie de signalisation est activé ou lorsque la cellule est dans un de ses états cellulaires. Déterminer si la signature représente bien la voie de signalisation nous permet de mieux interpréter et ou valider/invalider les résultats de la déconvolution."),
            p("Notre objectif principal reste de déterminer la proportion des états cellulaires ainsi que des voies de signalisation et d'étudier ces résultats."),
            p("Nous allons étudier des voies de signalisations qui sont issues des études ",strong("CancerSEA"), "et" ,strong("Hallmarks"),"."),
            
            titlePanel("CancerSEA Signatures"),
              
            htmlOutput("cancerSEA_sig_table")%>%
              withSpinner(type =6, color = "#7F00FF"),
            br(),
              
            titlePanel("Hallmarks Signatures"),
            htmlOutput("Hallmarks_sig_table")%>%
              withSpinner(type =6, color = "#7F00FF"),
            br(),
            
            p("Nous allons étudier les signatures de XBP1s et RIDD issues de la signature globale IRE138 ainsi que les signatures XBP1s pur et RIDD pur. Les chercheurs ont trouvé deux signatures différentes pour XBP1s et RIDD. La première, celle issue de IRE1 38, vient du fait que des chercheurs aient remarqué que 38 ARN était produit lorsque IRE1 était activé. Ils ont fait une liste du plus up regulated (qui augmente le plus en quantité d'ARN) au plus down regulated (qui baisse le plus en quantité d'ARN). Ils ont pris les 19 premiers (les plus up regulated) pour en faire la signature de XBP1s (qu'on va appeler ici XBP1 19) et les 19 derniers (les plus down regulated) pour en faire la signature de RIDD (qu'on va appeler ici RIDD 19)."),
            p("Pour les signatures de XBP1s pur et RIDD pur, on a pris un échantillon qu'on a privé d'environnement et on a 'activé' IRE1 de là, on a vu que 40 gènes étaient up regulated c'est ceux qu'on a désigné comme signature XBP1s et 37 gènes était down regulated, c'est ceux qu'on a désigné comme signature de RIDD."),
            
            titlePanel("XBP1S 19 et RIDD 19 Signatures"),
            htmlOutput("XBP1s_RIDD_19_sig_table")%>%
                withSpinner(type =6, color = "#7F00FF"),
            br(),
            
            titlePanel("XBP1S pur et RIDD pur Signatures"),
            htmlOutput("XBP1s_RIDD_pur_sig_table")%>%
                withSpinner(type =6, color = "#7F00FF"),
            br(),
            
            p("Afin de voir l'importance de nos différentes voies de signalisations dans nos échantillons nous allons utiliser une méthode appelée l'ICA qui est une méthode mathématique qui nous permet d'isoler des signaux latents dans un ensemble de signaux mixtes. Nous avons utilisé le package deconica d'Urszula Czerwinska pour faire l'ICA puis nous avons fait une gene set enrichment analysis qui nous a permis ensuite de déterminer un score d'abondance des états cellulaires et des voies de signalisations dans chaque échantillon."),
            p("Puis on va étudier ce que ses scores révèlent chez chaque patient et on va tenter de valider la qualité de nos signatures."),
            ),
    ),
          
  
  # Tab 2
  tabPanel("Déconvolution ICA",
           h2("Déconvolution de nos échantillons TCGA avec l'ICA"),
           p("Dans cette partie, nous allons donc chercher à estimer les proportions des différentes voies de signalisations que l'on a mentionnées précédemment au sein d'échantillons TCGA que l'on va importer."),
           
           DTOutput("resume_TCGA")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           
           p("Nous allons ici 'préparer' notre jeu de donnée pour pouvoir lancer avec notre algorithme d'ICA. Le nom des gènes de TCGA était du type ensembl et nos signatures ne le sont pas donc on passe tous les noms des gènes au même format que nos signatures."),
           
           p("Une ICA (Independent Component Analysis) est une technique non-déterministe pour séparer les signaux mixtes en leurs composantes indépendantes et non-gaussienne. Elle vise à trouver la transformation linéaire des données qui maximise l’indépendance mutuelle entre les composantes. Ainsi dans notre contexte de déconvolution cellulaire on part de nos données avec l'expression des gènes de chaque échantillon (appelée la matrice X en ICA) et on veut trouver deux matrices l'une avec les profils cellulaire et l'autre avec les proportion de ses cellules. Le terme profile cellulaire comme il est utilisé dans la thèse d'Urszula Czerwinska (et comme on l'utilise puisque nous avons été très fortement inspirés par sa thèse) désigne une classification des cellules basés sur leurs expressions transcriptomiques."),
           
           br(),
           
           img(src = "ICA.png", height = "500px"),
           
           br(),
           br(),
           
           p("Cette méthode n'est pas déterministe donc les résultats à la sortie sont différents. Pour minimiser ce problème, on va chercher à maximiser la stabilité des résultats en utilisant donc un autre algorithme qui trouve la stabilité pour un certains nombre de composantes (qui contient des profils cellulaire) que l'on choisit à la sortie de notre algorithme d'ICA."),
           
           p("Néanmoins, ici, je n'ai pas réussi à refaire l'algorithme. Le package deconica nous permet de le faire, mais il attend une version de matlab qui ne correspondait pas a celle que j'avais sur mon ordinateur."),
           
           p("Ainsi, on a estimé 'à la louche' qu'un bon nombre de composantes pour notre ICA serait 25 pour avoir plus de composantes que de voies de signalisations, mais" ,strong("ce n'est pas la procédure idéal à suivre et c'est un point à améliorer"),"."),
           
           p("On utilise donc l'algorithme de deconica run_fastica sur nos TCGA 'nettoyés' avec 25 composantes."),
           
           verbatimTextOutput("CodeCleanChunk")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("Lorsqu'on utilise le package deconica il est pratique de mettre nos signatures dans des dataframe avec des listes nommé. En effet, plusieurs fonctions du package ont besoin que les signatures soient de cette forme. Afin de continuer l'analyse chaque df est de la forme :df_plusieurs_sig <- (sig1=(gene = list_gene_1),sig2=(gene = list_gene-2),...,sign=(gene = list_gene_n))")
           
  ),
  
  # Tab 3
  tabPanel("XBP1 et RIDD",
           
          
           h2("Déconvolution des voies de signalisations XBP1s et RIDD au sein de nos échantillons TCGA"),
           
           p("Ici on s'intéresse à deux types de signatures pour XBP1s et RIDD:"),
           
           tags$ul(
             tags$li("XBP1s et RIDD issues de la signature globale d'IRE1 38"),
             tags$li("Les signatures XBP1s pur et RIDD pur")),
           br(),
           
           p("On va interpréter les résultats de notre déconvolution avec une analyse d'enrichissement. Elle va venir nous permettre d'associer une composante issues de l'ICA précédemment réalisée à une voie de signalisation/état cellulaire. L'analyse d'enrichissement ou 'gene set enrichment analysis' (GSEA) est un outil pour résumer l'information de tableaux d'expression de gène. Le GSEA nous permet d'avoir un tableau du type."),
           
           uiOutput("exemple_GSEA_XR")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           
           p(HTML("Ici la p-value est issue d'un test où H<sub>0</sub>  (l'hypothèse nulle) est, 'la distribution de cette liste n'est due qu'au hasard'. Elle suivrait donc une loi normale. Lorsque la p-value est en dessous d'un seuil que l'on se fixe (souvent 5 % ou 1 %) alors on rejette H<sub>0</sub> et on peut en conclure que ce gène à cette place dans la liste grâce à son expression (c'est-à-dire s'il est up ou down).")),
           
           p("Dans les graphiques de cette section, chaque point représente le résultat du test d'enrichissement, plus le point est gros, plus le résultat est significatif. La couleur indique la direction de l'enrichissement (rouge = enrichi, bleu = appauvri)."),
           
           p("",strong("NES")," = Normalized Enrichment Score, qui indique la force de l'enrichissement. Un NES positif indique que la signature est enrichie dans la composante ICA, tandis qu'un NES négatif indique que la signature est appauvrie. On a normalisé parce que la pvaleur est ajustée et pour tenir compte de la taille de la signature et du nombre total de gènes."),
           
           p("p-value ajustée (",strong("padj"),") = pvaleur ajustée pour tenir compte des tests multiples."),
           
           titlePanel("Table de résultats du GSEA pour XBP1s 19 et RIDD 19"),
           
           withLoader(uiOutput("result_GSEA_IRE1_38"),type = "image",loader =  "cat.gif"),
             
           br(),
           
           titlePanel("Table de résultat du GSEA pour XBP1s pur et RIDD pur"),
           
           uiOutput("result_GSEA_IRE1_pur")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           
           titlePanel("Représentation graphique des résultats de GSEA"),
           
           plotlyOutput("result_GSEA_IRE1_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           p("Afin de retrouver dans nos composantes issues de l'ICA les voies de signalisations qui nous intéressent ici. On va tout d'abord résumer l'information de chaque composante en ne prenant que ses 10 gènes les plus exprimés. De plus, on va calculer un score d'abondance pour chaque composante par échantillon (combien de fois on retrouve la composante dans chaque échantillon). Ainsi on va associer grâce a un test GSEA les voies de signalisations aux composantes qui les correspondent le plus et on dit que le score d'abondance de la composante dans l'échantillon correspond à la place de la voie de signalisation correspondante dans l'échantillon."),
           
           br(),
           
           titlePanel("Score d'abondance de XBP1s 19 et RIDD 19 au sein des échantillons TCGA"),
           checkboxGroupInput("selected_basenames", "Choisir les voies de signalisations :",
                              choices = all_basenames_38,
                              selected = all_basenames_38),
           
           plotlyOutput("result_ab_IRE138_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           
           titlePanel("Scores d'abondance de XBP1s pur et RIDD pur au sein des échantillons TCGA"),
           checkboxGroupInput("selected_basenames", "Choisir les voies de signalisations :",
                              choices = all_basenames_pur,
                              selected = all_basenames_pur),
           
           plotlyOutput("result_ab_IRE1pur_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           
           h2("Etude exploratoire de l'importance de nos différentes voies de signalisation au sein des échantillons"),
           
           p("Dans cette section, nous allons utiliser les résultats de notre ICA pour voir s'il y a des groupes d'échantillons qui ressortent particulièrement. On peut également voir s'il y a une tendance générale, et même peut-être tirer des conclusions sur la qualité des signatures XBP1s et RIDD."),
           
           h3("Etude de tendance générale au sein de nos groupes de voies de signalisation"),
           
           p("Ici, nous allons chercher à déterminer des groupes/des tendances avec tout d'abord une ACP pour voir le lien entre les différentes variables. Puis on va faire du clustering avec la CAH et les kmeans afin de voir si l'on ne retrouve pas des individus qui ont tendance à se regrouper."),
           
           h4("ACP"),
           
           titlePanel("Cercle des variables de l'ACP"),
           plotOutput("acp_ire1")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           p("Ici, on ne montre que le cercle des variables, car la représentation graphique des individus de l'ACP ne montrait aucune tendance particulière ou intéressante."),
           p("Sur le cercle des variables, on voit que les down sont liés entre eux et les up sont liés entre eux."),
           
           h4("kmeans"),
           p("Avant d'appliquer les kmeans on va chercher le nombre de groupes idéals qu'on doit utiliser (le k), pour cela on utilise la moyenne de divers indices avec la fonction Nbclust."),
           
           br(),
           layout_column_wrap(
             width = 1,
             value_box(
               title = "Nombre de groupe optimal selon plusieurs indices (silhouette, coude,etc...)",
               showcase = indice_coude,
               value = withSpinner(uiOutput("k_XR"), type = 6, color = "#7F00FF"),
               theme = "black",
               class = "custom-value-box"
             ),
             
           ),
           br(),
           titlePanel("Représentation graphique des groupes trouvés avec kmeans"),
           
           plotOutput("kmeans_ire1")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("Due à la nature non-déterministe de l'ICA les résultats peuvent être complètement différent d'une initialisation à une autre"),
           
           h4("CAH"),
           p("Avant d'appliquer la CAH on va chercher le nombre de groupes idéals qu'on doit utiliser (le k), pour cela on utilise la moyenne de divers indices avec la fonction Nbclust."),
           
           br(),
           layout_column_wrap(
             width = 1,
             value_box(
               title = "Nombre de groupe optimal selon plusieurs indices (silhouette, coude,etc...)",
               showcase = indice_coude,
               value = withSpinner(uiOutput("k_XR_CAH"), type = 6, color = "#7F00FF"),
               theme = "black",
               class = "custom-value-box"
             ),
             
           ),
           br(),
           titlePanel("Représentation graphique des groupes trouvés avec la CAH"),
           
           plotOutput("cah_ire1")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("Due à la nature non-déterministe de l'ICA les résultats peuvent être complètement différent d'une initialisation à une autre"),
           
           p("",strong("En conclusion : "),"On ne voit pas de groupe se former à l'issue de notre déconvolution, mais l'on voit que les variables down régulés sont liés entre elle, de même pour les up régulés.")
           
           
           
    
  ),
  
  
  # Tab 4
  tabPanel("CancerSEA",
           h2("Déconvolution des voies de signalisations issues de CancerSEA au sein de nos échantillons TCGA"),
           
           p("Ici, nous nous intéressons aux signatures des voies de signalisation qu'on retrouve dans Cancer SEA sont : Angiogenèse, Apoptose, Cycle cellulaire, Differentiation, Dégradation de l'ADN, Réparation de l'ADN, l'EMT, l'Hypoxie, l'Inflamation, l'Invasion, la Metastase, la Proliferation, Quiescence, Stemness (La souchitude d'une cellule)"),
           
           p("On va interpréter les résultats de notre déconvolution avec une analyse d'enrichissement. Elle va venir nous permettre d'associer une composante issue de l'ICA précédemment réalisée à une voie de signalisation/état cellulaire. L'analyse d'enrichissement ou 'gene set enrichment analysis' (GSEA) est un outil pour résumer l'information de tableaux d'expression de gènes. Le GSEA nous permet d'avoir un tableau du type."),
           
           uiOutput("exemple_GSEA_CS")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p(HTML("Ici la p-value est issue d'un test où H<sub>0</sub>  (l'hypothèse nulle) est, 'la distribution de cette liste n'est due qu'au hasard'. Elle suivrait donc une loi normale. Lorsque la p-value est en dessous d'un seuil que l'on se fixe (souvent 5 % ou 1 %) alors on rejette H<sub>0</sub> et on peut en conclure que ce gène à cette place dans la liste grâce à son expression (c'est-à-dire s'il est up ou down).")),
           
           p("Dans les graphiques de cette section, chaque point représente le résultat du test d'enrichissement, plus le point est gros, plus le résultat est significatif. La couleur indique la direction de l'enrichissement (rouge = enrichi, bleu = appauvri)."),
           
           p("",strong("NES")," = Normalized Enrichment Score, qui indique la force de l'enrichissement. Un NES positif indique que la signature est enrichie dans la composante ICA, tandis qu'un NES négatif indique que la signature est appauvrie. On a normalisé parce que la pvaleur est ajustée et pour tenir compte de la taille de la signature et du nombre total de gènes."),
           
           p("p-value ajustée (",strong("padj"),") = pvaleur ajustée pour tenir compte des tests multiples."),
           
           titlePanel("Table de résultats du GSEA pour les signatures de CancerSEA"),
           
           uiOutput("result_GSEA_SEA")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           titlePanel("Représentation graphique des résultats de GSEA"),
           
           plotlyOutput("result_GSEA_SEA_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
           br(),
           
           p("Afin de retrouver dans nos composantes issues de l'ICA les voies de signalisations qui nous intéressent ici. On va tout d'abord résumer l'information de chaque composante en ne prenant que ses 10 gènes les plus exprimés. De plus, on va calculer un score d'abondance pour chaque composante par échantillon (combien de fois on retrouve la composante dans chaque échantillon). Ainsi on va associer grâce a un test GSEA les voies de signalisations aux composantes qui les correspondent le plus et on dit que le score d'abondance de la composante dans l'échantillon correspond à la place de la voie de signalisation correspondante dans l'échantillon."),
           titlePanel("Scores d'abondance des signatures de CancerSEA au sein des échantillons TCGA"),
           pickerInput("selected_basenames", "Veuillez sélectionner une pathway :",
                       choices = all_basenames_SEA,
                       selected = all_basenames_SEA,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE,
                                      `live-search` = TRUE,
                                      `selected-text-format` = "count > 3")),
           
           
           plotlyOutput("result_ab_SEA_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           h2("Etude de l'importance de nos différentes voies de signalisation au sein des échantillons"),
           
           p("Dans cette section, nous allons utiliser les résultats de notre ICA pour voir s'il y a des groupes d'échantillons qui ressortent particulièrement. On peut également voir s'il y a une tendance générale."),
           
           h3("Etude de tendance générale au sein de nos groupes de voies de signalisation"),
           
           p("Ici, nous allons chercher à déterminer si on trouve des groupes/des tendances au sein de nos différents ensembles de voies de signalisations.

On va tout d'abord procéder à une ACP pour voir le lien entre les différentes variables. Puis on va faire du clustering avec la CAH et les kmeans afin de voir si l'on ne retrouve pas des individus qui ont tendance à se regrouper."),
           
           h4("ACP"),
           
           titlePanel("Représentation graphiques des individus "),
           
           plotOutput("ind_SEA")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("",strong("Interprétation "),": On ne voit pas de groupe particulier."),
           
           titlePanel("Représentation graphiques des variables "),
           
           plotOutput("var_SEA")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("",strong("Interprétation :"),"Là, on voit que toutes les variables sont au moins corrélées à une ou plusieurs variables. Mais on ne remarque pas vraiment de groupe non plus."),
           
           h4("kmeans"),
           
           p("Avant d'appliquer les kmeans on va chercher le nombre de groupes idéals qu'on doit utiliser (le k), pour cela, on utilise la moyenne de divers indices avec la fonction Nbclust."),
           
           br(),
           layout_column_wrap(
             width = 1,
             value_box(
               title = "Nombre de groupe optimal selon plusieurs indices (silhouette, coude,etc...)",
               showcase = indice_coude,
               value = withSpinner(uiOutput("k_SEA"), type = 6, color = "#7F00FF"),
               theme = "black",
               class = "custom-value-box"
             ),
             
           ),
           
           
           
           titlePanel("Représentation graphique des groupes trouvés avec kmeans"),
           
           plotOutput("kmeans_SEA")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("Due à la nature non-déterministe de l'ICA les résultats peuvent être complètement différent d'une initialisation à une autre"),
           
           h4("CAH"),
           p("Avant d'appliquer la CAH on va chercher le nombre de groupes idéals qu'on doit utiliser (le k), pour cela, la moyenne de divers indices avec la fonction Nbclust."),
           
           br(),
           layout_column_wrap(
             width = 1,
             value_box(
               title = "Nombre de groupe optimal selon plusieurs indices (silhouette, coude,etc...)",
               showcase = indice_coude,
               value = withSpinner(uiOutput("k_SEA_CAH"), type = 6, color = "#7F00FF"),
               theme = "black",
               class = "custom-value-box"
             ),
             
           ),
           br(),
           titlePanel("Représentation graphique des groupes trouvés avec la CAH"),
           
           plotOutput("cah_SEA")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("Due à la nature non-déterministe de l'ICA les résultats peuvent être complètement différent d'une initialisation à une autre"),
           
           
  
           
           
  ),
  tabPanel("Hallmarks",
           
           h2("Déconvolution des voies de signalisations issues de Hallmarks au sein de nos échantillons TCGA"),
           
           p("Les signatures des voies de signalisation qu'on retrouve dans Hallmarks sont : angiogenesis, cell replicative immortality, change of cellular energetics, escaping immune response, escaping programmed cell death, genome instability and mutations, invasion and metastasis, proliferative signalling, suppression of growth, tumour promoting inflamation"),
           
           p("On va interpréter les résultats de notre déconvolution avec une analyse d'enrichissement. Elle va venir nous permettre d'associer une composante issue de l'ICA précédemment réalisée à une voie de signalisation/état cellulaire. L'analyse d'enrichissement ou 'gene set enrichment analysis' (GSEA) est un outil pour résumer l'information de tableaux d'expression de gène. Le GSEA nous permet d'avoir un tableau du type."),
           
           uiOutput("exemple_GSEA_HM")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p(HTML("Ici la p-value est issue d'un test où H<sub>0</sub>  (l'hypothèse nulle) est, 'la distribution de cette liste n'est due qu'au hasard'. Elle suivrait donc une loi normale. Lorsque la p-value est en dessous d'un seuil que l'on se fixe (souvent 5 % ou 1 %) alors on rejette H<sub>0</sub> et on peut en conclure que ce gène à cette place dans la liste grâce à son expression (c'est-à-dire s'il est up ou down).")),
           
           p("Dans les graphiques de cette section, chaque point représente le résultat du test d'enrichissement, plus le point est gros, plus le résultat est significatif. La couleur indique la direction de l'enrichissement (rouge = enrichi, bleu = appauvri)."),
           
           p("",strong("NES")," = Normalized Enrichment Score, qui indique la force de l'enrichissement. Un NES positif indique que la signature est enrichie dans la composante ICA, tandis qu'un NES négatif indique que la signature est appauvrie. On a normalisé parce que la pvaleur est ajustée et pour tenir compte de la taille de la signature et du nombre total de gènes."),
           
           p("p-value ajustée (",strong("padj"),") = pvaleur ajustée pour tenir compte des tests multiples."),
           
           titlePanel("Table de résultat du GSEA pour les signatures de Hallmarks"),
           
           uiOutput("result_GSEA_Hallmark")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           titlePanel("Représentation graphique des résultats de GSEA"),
           
           plotlyOutput("result_GSEA_Hallmark_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
            br(),
           
           p("Afin de retrouver dans nos composantes issues de l'ICA les voies de signalisations qui nous intéresse ici. On va tout d'abord résumer l'information de chaque composante en ne prenant que ses 10 gènes les plus exprimés. De plus on va calculer un score d'abondance pour chaque composante par échantillon (combien de fois on retrouve la composante dans chaque échantillon). Ainsi on va associer grace a un test GSEA les voies de signalisations aux composantes qui les correspondent le plus et on dit que le score d'abondance de la composante dans l'échantillon correspond à la place de la voie de signalisation correspondante dans l'échantillon."),
           
          titlePanel("Score d'abondance des signature Hallmarks au sein des échantillons TCGA"),
           pickerInput("selected_basenames", "Veuillez sélectionner une pathway :",
                       choices = all_basenames_Hallmarks,
                       selected = all_basenames_Hallmarks,
                       multiple = TRUE,
                       options = list(`actions-box` = TRUE,
                                      `live-search` = TRUE,
                                      `selected-text-format` = "count > 3")),
           
           plotlyOutput("result_ab_Hallmarks_graph")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           h3("Etude de tendance générale au sein de nos groupes de voies de signalisation"),
           
           p("Ici, nous allons chercher à déterminer si on trouve des groupes/des tendances au sein de nos différents ensembles de voies de signalisations.

On va tout d'abord procéder à une ACP pour voir le lien entre les différentes variables. Puis on va faire du clustering avec la CAH et les kmeans afin de voir si l'on ne retrouve pas des individus qui ont tendance à se regrouper."),
           
           h4("ACP"),
           
           titlePanel("Représentation graphiques des individus "),
           
           plotOutput("ind_Hallmarks")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("",strong("Inteprétation"),"On ne voit pas de groupe particulier"),
           
           titlePanel("Représentation graphiques des variables "),
           
           plotOutput("var_Hallmarks")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("",strong("Inteprétation"),"Là, on voit que toutes les variables sont au moins corrélées à une ou plusieurs variables. Mais on ne remarque pas vraiment de groupe non plus."),
           
           h4("kmeans"),
           
          p("Avant d'appliquer les kmeans on va chercher le nombre de groupes idéals qu'on doit utiliser (le k), pour cela, la moyenne de divers indices avec la fonction Nbclust."),
          
          br(),
          layout_column_wrap(
            width = 1,
            value_box(
              title = "Nombre de groupe optimal selon plusieurs indices (silhouette, coude,etc...)",
              showcase = indice_coude,
              value = withSpinner(uiOutput("k_H"), type = 6, color = "#7F00FF"),
              theme = "black",
              class = "custom-value-box"
            ),
            
          ),
          br(),
           
          
           
           titlePanel("Représentation graphique des groupes trouvés avec kmeans"),
           
           plotOutput("kmeans_Hallmarks")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
          p("Due à la nature non-déterministe de l'ICA les résultats peuvent être complètement différent d'une initialisation à une autre"),
           
           h4("CAH"),
          p("Avant d'appliquer la CAH on va chercher le nombre de groupes idéals qu'on doit utiliser (le k), pour cela, la moyenne de divers indices avec la fonction Nbclust."),
          br(),
          layout_column_wrap(
            width = 1,
            value_box(
              title = "Nombre de groupe optimal selon plusieurs indices (silhouette, coude,etc...)",
              showcase = indice_coude,
              value = withSpinner(uiOutput("k_H_CAH"), type = 6, color = "#7F00FF"),
              theme = "black",
              class = "custom-value-box"
            ),
            
          ),
          br(),
           titlePanel("Représentation graphique des groupes trouvés avec la CAH"),
           
           plotOutput("cah_Hallmarks")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("",strong("Inteprétation"),"On ne voit pas de groupes se former particulièrement"),
          
  
           
           ),
  tabPanel("Qualité signature XBP1s et RIDD",
           h2("Etude et comparaison de la qualité de la signature de XBP1s et RIDD"),
           p("Nous allons ici entamer un travail pour voir si les signatures XBP1s pur et RIDD pur et XBP1s 19 et RIDD 19 s'activent vraiment lorsque les voies de signalisations XBP1s et RIDD s'activent. On va également regarder si l'une des signatures est plus performante que l'autre à nous montrer que XBP1s et RIDD sont activés. "),
           
           h3("Tableau de contingence"),
           p("Tout d'abord, nous allons voir les effectifs de chaque signature au sein de nos échantillons TCGA."),
           
           titlePanel("Tableau de contingence des occurrences de chaque signature"),
           uiOutput("contingence_occurence")%>%
             withSpinner(type =6, color = "#7F00FF"),
           p("",strong("Note : "),"Il y a 173 patients dans notre jeu de donnés"),
           br(),
           titlePanel("Tableau de contingence des scores cumulés de chaque signature"),
           uiOutput("contingency_scores")%>%
             withSpinner(type =6, color = "#7F00FF"),
           h3("Distribution des scores lorsqu'IRE1 est actif"),
           p("Ici, nous allons simplement utiliser les scores des signatures XBP1s up et RIDD down, car ce sont des signes qu'IRE1 est activés. Or, nos échantillons sont des patients atteints du Cancer donc on s'attendrait à ce qu'IRE1 soit activé.

Nous allons donc simplement faire la somme des scores de XBP1 UP et RIDD DOWN pour chaque échantillon."),
           titlePanel("Distribution des scores de XBP1s 19 et RIDD 19 lorsqu'IRE1 est actif"),
           plotlyOutput("activite_38")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           titlePanel("Distribution des scores de XBP1s pur et RIDD pur lorsqu'IRE1 est actif"),
           
           plotlyOutput("activite_pur")%>%
             withSpinner(type =6, color = "#7F00FF"),
           
           p("On voit qu'IRE1 parait beaucoup plus actif lorsqu'on utilise les signatures issues de la signature globale d'IRE1 38 que lorsqu'on utilise les signatures de XBP1s pur et RIDD pur. En effet, sa distribution de score est bien plus élevée que celle de XBP1s pur et RIDD pur."),
           
           p("On peut le visualiser ici :"),
           titlePanel("Boxplots de la distribution des scores des deux signatures différentes de XBP1s et RIDD lorsqu'IRE1 est actif"),
           plotOutput("comparaison_activite")%>%
             withSpinner(type =6, color = "#7F00FF"),
           p("",strong("Légende :"),"en orange, c'est les signatures XBP1s 19 et RIDD 19 et en turquoise, c'est XBP1s pur et RIDD pur."),
           
           h3("Etude de la cohérence biologique au sein des scores"),
           p("Ici, on va chercher à savoir si les scores qu'on a attribué aux échantillons respectent une cohérence biologique. En effet, on s'attend à ce qu'un patient qui a beaucoup de cellules dans sa tumeur avec XBP1s UP régulé ai peu de cellule DOWN régulé et vice versa. Ainsi quand un échantillon a un score XBP1s UP élevé il est censé avoir un score XBP1s DOWN faible et de même pour l'inverse. Ainsi nous avons élaboré un indice qui va regarder si les rangs de XBP1s UP et DOWN sont inversés et vice versa. "),
           titlePanel("indice d'inversion de rang (basé sur la corrélation)"),
           
           
           layout_column_wrap(
             width = 1/2,
             value_box(
               title = "Indice inversion de rang entre XBP1S 19 DOWN et UP",
               showcase = calculator_fill,
               value = withSpinner(uiOutput("my_box1"), type = 6, color = "#7F00FF"),
               theme = "black",
               class = "custom-value-box"
             ),
             value_box(
               title = "Indice inversion de rang entre XBP1S pur DOWN et UP",
               showcase = calculator_fill,
               value = withSpinner(uiOutput("my_box2"), type = 6, color = "#7F00FF"),
               theme = "black",
               class = "custom-value-box"
             )
            ),
           br(),
           
           p("Cet indice nous permet de voir à quel point deux listes de même taille ont un ordre inversé de rangs. On fait cela parce qu'on attendrait normalement que pour un échantillon lorsque XBP1S DOWN est fortement exprimé alors XBP1s UP est faiblement exprimé.  Autrement dit on s'attendrait qu'un échantillon qui est un rang haut dans la liste des score XBP1s UP est un rang bas dans la liste des score XBP1s DOWN. Ainsi lorsque notre indice vaut -1 les deux listes sont parfaitement dans le même ordre, lorsqu'il vaut 0 il n'y a aucun lien entre les ordres et lorsqu'il vaut 1 alors les listes sont parfaitement inversés. "),
           
           p("Nos indices ici sont proches de zéro, ce qui montre qu'il n'y pas de lien entre les rangs des scores de nos signatures lorsqu'elles sont UP ou DOWN. Ainsi cela prouve que les résultats n'ont pas de cohérence biologique. Cela peut être due à plusieurs facteurs : la méthode utilisé pour la déconvolution, la méthode utilisé pour élaborer les scores et/ou la qualité des signatures. "),
           
          
          
          #p("Ici pour représenter le bruit je vais vous montrer l'inversion de rangs pour les signatures issues d'IRE1 38 et celle de XBP1s pur et RIDD pur."),
          
          #layout_column_wrap(
          #plotOutput("inversion_rang_38")%>%
          #  withSpinner(type =6, color = "#7F00FF"),
          #plotOutput("inversion_rang_pur")%>%
          #  withSpinner(type =6, color = "#7F00FF"),
          #),
          
          h3("Etude de nos échantillons lorsqu'IRE1 est actif"),
          p("Lorsqu'IRE1 est actif XBP1s est UP et RIDD est DOWN on va étudier les échantillons où l'on retrouve cela. Et on va tenter de refaire la pipeline précédement étudié en ne gardant que les composantes qui pour le même type de signature ont de manière significative du XBP1 UP lorsque que le RIDD est DOWN ou l'inverse."),
          p("Attention due a la nature non déterministe de l'ICA il se peut que cette section soit parfaitement vide car aucun des composante ne correspond à ce critère tout simplement."),
          titlePanel("Représentation graphique de la gsea pour ire1 actif"),
          plotlyOutput("gsea_ire1_actif")%>%
            withSpinner(type =6, color = "#7F00FF"),
          titlePanel("Distribution des score ire1 actif"),
          plotlyOutput("score_ire1actif")%>%
            withSpinner(type =6, color = "#7F00FF"),
          
        titlePanel("Tableau de contingence des occurences"),
          uiOutput("contingence_occurence_actif")%>%
            withSpinner(type =6, color = "#7F00FF"),
        p("",strong("Note : "),"Il y a 173 patients dans notre jeu de donnés"),
          titlePanel("Tableau des contingence des scores"),
          uiOutput("contingence_score_actif")%>%
            withSpinner(type =6, color = "#7F00FF"),
          
          
  )
  
  ),
    
  ),

  )
  
  
    
    
  

