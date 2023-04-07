require(shiny)

ui <- fluidPage(
  theme = shinythemes::shinytheme("united"),
  titlePanel("PBMC data in epileptic patients and controls"),

  tabsetPanel(
    tabPanel("Data Overview", h3("Select cluster"),
             sidebarLayout(
               sidebarPanel(
                 radioButtons("Variable", h4("Select variable to visualize"),
                              choices = list("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","Has TCR 10x data"),
                              selected = "Cell Type (L1)"),
                 checkboxGroupInput("inCheckboxGroup", "Input checkbox",
                                    unlist(vals), selected = unlist(vals)[1]),
                 downloadButton('downloadPlot_UMAP', 'Download UMAP Plot')

               ),
               mainPanel(plotOutput("UMAP")),
             )

    ),
    tabPanel("Marker genes",
             h3("Cell Type (Level 1)"),
             fluidRow(
               column(2,
                      selectInput("DEG_L1", h4("Cell type"),
                                  choices = lapply(names(roc_l1_groups),function(x) x), selected = 1),
                      hr(),
                      selectInput("DEG_L1_marker", h4("Cell type marker"),
                                  choices = lapply(names(roc_l1_genes),function(x) x)),
                      hr(),
                      downloadButton('downloadPlot_UMAP_Vln_feature_L1', 'Download Plots'),
                      hr(),
                      downloadButton('downloadTable_ROC_results_L1', 'Download ROC Table')

               ),
               column(3,plotOutput("UMAP_L1")
               ),
               column(3,plotOutput("UMAP_feature_L1")
               ),
               column(4,plotOutput("Vln_feature_L1")
               )),hr(),
             dataTableOutput("DEG_L1_Table"),
             hr(),
             h3("Cell Type (Level 2)"),
             fluidRow(
               column(2,
                      selectInput("DEG_L2", h4("Cell sub-type"),
                                  choices = lapply(names(roc_l2_groups),function(x) x), selected = 1),
                      hr(),
                      selectInput("DEG_L2_marker", h4("Cell sub-type marker"),
                                  choices = lapply(names(roc_l2_genes),function(x) x)),
                      hr(),
                      downloadButton('downloadPlot_UMAP_Vln_feature_L2', 'Download Marker Plots'),
                      hr(),
                      downloadButton('downloadTable_ROC_results_L2', 'Download ROC Table')

               ),
               column(3,plotOutput("UMAP_L2")
               ),
               column(3,plotOutput("UMAP_feature_L2")
               ),
               column(4,plotOutput("Vln_feature_L2")
               )),hr(),
             dataTableOutput("DEG_L2_Table")
    ),
    tabPanel("MASC Analysis",
             h3("Scroll through table for specific statistics"),
             fluidRow(column(12,plotOutput("MASC"))),
            hr(),
                       fluidRow(
                         column(3,
                                radioButtons("Cell_Type_Level",h4("Cell Type Level"),
                                             choices = list("L1","L2"),
                                             selected = "L1"),hr(),
                                downloadButton('downloadPlot_MASC', 'Download Plot'),hr(),
                                downloadButton('downloadTable_MASC', 'Download Table')
                         ),
                         column(9,dataTableOutput("MASC_results")))),
    tabPanel("Healthy controls vs epileptic patients",
             h3("Select major cell type"),
             fluidRow(
               column(2,
                      radioButtons("HC_vs_EP_Level",label = h4("Cell Type Level"),
                                   choices = c("L1","L2"), selected = "L1"), hr(),
                      radioButtons("HC_Epil_Cluster",label="Cell type",
                                  choices=lapply(names(roc_l1_groups),function(x) x) , selected = 1),
                      hr(),
                      downloadButton('downloadPlots_HC_vs_EP', 'Download Plots'),
               ),
               column(10,plotOutput("Volcano_HC_EP")
               )),
             fluidRow(
               column(2,
                      radioButtons("Vln_Cols", h4("Number of columns"),
                                   choices= lapply(1:5, function(x) x), selected=4)),
               column(10,plotOutput("Vln_Plots_HC_EP"))

             )),
    tabPanel("Group pairwise comparisons",
             h3("Broad cell type"),
             fluidRow(
               column(2,
                      radioButtons("L1_Pairwise", label = h4("Comparison"),
                                   choices = c("HC vs MC","HC vs MR","MR vs MC"),
                                   selected = "HC vs MR"),
                      hr(),
                      radioButtons("Direction_L1",h4("Direction"),
                                   choices = c("Up-regulated","Down-regulated"),
                                   selected = "Down-regulated"),
                      hr(),
                      radioButtons("DEG_L1_choice", h4("Cell Type"),
                                  choices= l1_names,
                                  selected = "B"),
                      hr(),
                      downloadButton("downloadPlot_Pairwise_L1",'Download Plot')),
               column(10,
                      plotOutput("DEG_Violins_L1"),
                      dataTableOutput("DEG_L1_Table_Pair"))
             ),
             hr(),
             h3("Sub cell type"),
             fluidRow(
               column(2,
                      radioButtons("L2_Pairwise", label = h4("Comparison"),
                                  choices = c("HC vs MC","HC vs MR","MR vs MC"),
                                  selected = "HC vs MR"),
                      hr(),
                      radioButtons("Direction_L2",h4("Direction"),
                                  choices = c("Up-regulated","Down-regulated"),
                                  selected = "Down-regulated"),
                      hr(),
                      radioButtons("DEG_L2_choice", h4("Sub Cell Type"),
                                  choices= l2_names,
                                  selected = "Bintermediate"),
                      hr(),
                      downloadButton("downloadPlot_Pairwise_L2",'Download Plot')
               ),
               column(10,
                      plotOutput("DEG_Violins_L2"),
                      dataTableOutput("DEG_L2_Table_Pair"))
             )),
    tabPanel("Tau associations",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("Variable_tau", h4("Select variable to visualize"),
                              choices = list("Cell Type (L1)",
                                             "Cell Type (L2)"),
                              selected = "Cell Type (L1)"),
                 checkboxGroupInput("inCheckboxGroup_tau", "Input checkbox",
                                    unlist(vals[2:3]), selected = unlist(vals[2:3])[1]),
                 downloadButton('downloadPlot_Tau', 'Download Plot')

               ),
               mainPanel(plotOutput("Tau"))
             )

    ),
    tabPanel("Explore markers",
             fluidRow(
               column(6,textInput("marker_A",h3("Gene A"),value = "CD8A"),
                      plotOutput("DimPlot_A"),
                      hr(),
                      selectInput("marker_A_variable",h4("Variable of interest"),
                                  choices= c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","Has TCR 10x data"),
                                  selected = "Condition"),
                      plotOutput("Vln_Plot_A"),
                      hr(),
                      downloadButton("dowloadPlot_Grid","Download Figures")),
               column(6,textInput("marker_B",h3("Gene B"),value = "CD79A"),
                      plotOutput("DimPlot_B"),
                      hr(),
                      selectInput("marker_B_variable",h4("Variable of interest"),
                                  choices =c("Cell Type (L1)","Cell Type (L2)",
                                             "Sex", "Condition",
                                             "Seurat Clusters","Has TCR 10x data"),
                                  selected = "Condition" ),
                      plotOutput("Vln_Plot_B"))
             ))
  )
)
