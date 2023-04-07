require(shiny)

server <- function(input, output, session) {

  session$onSessionEnded(function() {
    stopApp()
  })

  # vals

  observe({
    var <- switch(input$Variable,
                  "Cell Type (L1)" = 1,
                  "Cell Type (L2)" = 2,
                  "Sex" = 3,
                  "Condition" = 4,
                  "Seurat Clusters" = 5,
                  "Has TCR 10x data" = 6)
    x <- vals[[var]]
    # Can also set the label and select items
    updateCheckboxGroupInput(session, "inCheckboxGroup",
                             label = "Variable of interest options",
                             choices = x,
                             selected = x[1]
    )
  })

  observe({
    var_tau <- switch(input$Variable_tau,
                      "Cell Type (L1)" = 1,
                      "Cell Type (L2)" = 2)
    x_var <- vals[[var_tau]]
    # Can also set the label and select items
    updateCheckboxGroupInput(session, "inCheckboxGroup_tau",
                             label = "Cell type level",
                             choices = x_var,
                             selected = x_var[1]
    )
  })

  observe({
    var_marker_L1 <- match(input$DEG_L1, names(roc_l1_groups))
    x_marker_L1 <- roc_l1_genes_list[[var_marker_L1]]
    # Can also set the label and select items
    updateSelectInput(session, "DEG_L1_marker",
                             label = "Cell type markers",
                             choices = x_marker_L1,
                             selected = x_marker_L1[1]
    )
  })

  observe({
    var_marker_L2 <- match(input$DEG_L2, names(roc_l2_groups))
    x_marker_L2 <- roc_l2_genes_list[[var_marker_L2]]
    # Can also set the label and select items
    updateSelectInput(session, "DEG_L2_marker",
                      label = "Cell sub-type markers",
                      choices = x_marker_L2,
                      selected = x_marker_L2[1]
    )
  })

  observe({

    if(input$HC_vs_EP_Level =="L1"){
      x_var <- unique(hc_vs_ep_df_l1 %>% filter(p_adj <= 0.05) %>% pull(cell_type) %>% unique())
    } else{
      x_var <- unique(hc_vs_ep_df_l2 %>% filter(p_adj <= 0.05) %>% pull(cell_type) %>% unique())
    }
    # Can also set the label and select items
  updateRadioButtons(session, "HC_Epil_Cluster",
                             label = "Cell types",
                             choices = x_var,
                             selected = x_var[1]
    )
  })

  sub_types <- reactive({
    input$inCheckboxGroup
  })
  sub_types_tau <- reactive({
    input$inCheckboxGroup_tau
  })

  sub_types_epil_ctrl <- reactive({
    input$HC_Epil_Cluster
  })

  sub_types_DEG_L1 <- reactive({
    input$DEG_L1_marker
  })
  sub_types_DEG_L2 <- reactive({
    input$DEG_L2_marker
  })


  ###### UMAP ###
  UMAP_plot <- function(){
    var <- switch(input$Variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "gender",
                  "Condition" = "condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "Has TCR 10x data" = "is_t_cell")
    pbmc$selected_column <- as.character(pbmc@meta.data[[var]])
    if(!("ALL" %in% sub_types())){
      pbmc$selected_column[!(pbmc$selected_column %in% sub_types())] <- "Unselected"
    }
    gg <- DimPlot(pbmc, group.by = "selected_column",label=T)+ ggtitle(input$Variable)
    return(gg)
  }

  output$UMAP <- renderPlot({
    UMAP_plot()
  }, height= 500,width=600)

  output$downloadPlot_UMAP <- downloadHandler(
    filename = "UMAP.png",
    content = function(file) {
      png(file, height= 8, width = 9,res=1200, units = "in")
      print(UMAP_plot())
      dev.off()
    }
  )

####################  Marker genes

  umap_feature_L1_plot <- function(){
    plot <- FeaturePlot(pbmc, features = sub_types_DEG_L1())
    return(plot)
  }

  vln_feature_L1_plot <- function(){
    plot <- VlnPlot(pbmc, features = sub_types_DEG_L1(), group.by="predicted.celltype.l1")+NoLegend()
    return(plot)
  }

  DEG_L1_table <- function(){
    return(roc_l1 %>% filter(cluster == input$DEG_L1))
  }

  output$UMAP_L1 <- renderPlot({
    pbmc$selected_column <- paste0("Not ", input$DEG_L1)
    pbmc$selected_column[pbmc$predicted.celltype.l1==input$DEG_L1] <- input$DEG_L1
    gg <- DimPlot(pbmc, group.by = "selected_column",label=T) +ggtitle("")+NoLegend()
    return(gg)
  })

  output$UMAP_feature_L1 <- renderPlot({
    umap_feature_L1_plot()
  })
  output$Vln_feature_L1 <- renderPlot({
    vln_feature_L1_plot()
  })
  output$DEG_L1_Table <- renderDataTable({
    DEG_L1_table()
  })

  output$downloadPlot_UMAP_Vln_feature_L1 <- downloadHandler(
    filename = paste0("UMAP_Vln_Feature",sub_types_DEG_L1(),".png"),
    content = function(file) {
      png(file, height= 8, width = 10,res=1200, units = "in")
      print(umap_feature_L1_plot() + vln_feature_L1_plot() )
      dev.off()
    }
  )

  output$downloadTable_ROC_results_L1 <- downloadHandler(
    filename = paste0(input$DEG_L1,"_Marker_Genes_table.csv"),
    content = function(file) {
      write.csv(DEG_L1_table(), file)
    }
  )
#####Marker L2###

  umap_feature_L2_plot <- function(){
    plot <- FeaturePlot(pbmc, features = sub_types_DEG_L2())
    return(plot)
  }

  vln_feature_L2_plot <- function(){
    plot <- VlnPlot(pbmc, features = sub_types_DEG_L2(), group.by="predicted.celltype.l2")+NoLegend()
    return(plot)
  }

  DEG_L2_table <- function(){
    return(roc_l2 %>% filter(cluster == input$DEG_L2))
  }

  output$UMAP_L2 <- renderPlot({
    pbmc$selected_column <- paste0("Not ", input$DEG_L2)
    pbmc$selected_column[pbmc$predicted.celltype.l2==input$DEG_L2] <- input$DEG_L2
    gg <- DimPlot(pbmc, group.by = "selected_column",label=T)+ggtitle("")+NoLegend()
    return(gg)
  })

  output$UMAP_feature_L2 <- renderPlot({
    umap_feature_L2_plot()
  })
  output$Vln_feature_L2 <- renderPlot({
    vln_feature_L2_plot()
  })
  output$DEG_L2_Table <- renderDataTable({
    DEG_L2_table()
  })

  output$downloadPlot_UMAP_Vln_feature_L2 <- downloadHandler(
    filename = paste0("UMAP_Vln_Feature",sub_types_DEG_L2(),".png"),
    content = function(file) {
      png(file, height= 8, width = 10,res=1200, units = "in")
      print(umap_feature_L2_plot() + vln_feature_L2_plot() )
      dev.off()
    }
  )

  output$downloadTable_ROC_results_L2 <- downloadHandler(
    filename = paste0(input$DEG_L2,"_Marker_Genes_table.csv"),
    content = function(file) {
      write.csv(DEG_L2_table(), file)
    }
  )

  ######## MASC ###

  MASC_plot <- function(){
    plot <- masc_df %>%
      filter(level == input$Cell_Type_Level) %>%
      dplyr::mutate(cluster = paste0(str_remove(cluster,"cluster"),"\n(",size,")")) %>%
      dplyr::mutate(ref = factor(ref, levels = c(0,1,2), labels = c("Control","Medically Controlled","Medically Refractory"))) %>%
      dplyr::filter(size  > 100) %>%
      ggplot(aes(x=cluster, y = diagnosis1.OR, color = ref, group= ref))+
      geom_point(aes(size = Freq/100), position = position_dodge(width=.7))+
      geom_errorbar(aes(ymin=diagnosis1.OR.95pct.ci.lower,
                        ymax=diagnosis1.OR.95pct.ci.upper), color = "black",width=.2, position=position_dodge(width=.7))+
      # scale_color_brewer(palette="Set1")+
      theme_bw()+
      ylab("OR")+
      coord_cartesian(ylim=c(0,3))+
      geom_hline(yintercept = 1, linetype = "dashed")+
      facet_wrap(~level, scales = "free", nrow=2)+
      theme(axis.text.x = element_text(angle = 35, vjust=.5))+
      scale_color_brewer(palette="Set1", name = "Condition")

    return(plot)
  }

  MASC_table <- function(){
    df <- masc_df %>%
      filter(level==input$Cell_Type_Level)
    colnames(df) <- c("Cell Type","size","model.pval","diagnosis.ref.OR","OR.ci.lower","OR.ci.upper","level","Ref","n")
    df <- df %>%
      dplyr::select(-level, -size) %>%
      mutate(Ref = factor(Ref, levels= 0:2,labels =c("Healthy Control","Med. Controlled","Med. Refractory"))) %>%
      dplyr::rename('signif' = 'model.pval') %>%
      mutate(signif = "n.s.") %>%
      relocate(Ref, `Cell Type`,n)
    df[,5:7] <-round(as.matrix(df[,5:7]),4)
    df$signif[!sapply(1:nrow(df), function(i) between(1, df$OR.ci.lower[i],df$OR.ci.upper[i]))] <- "*"
    return(df)
  }

  output$MASC <- renderPlot({
    MASC_plot()
  }, height=400,width=900)

  output$MASC_results <- renderDataTable({
    MASC_table()
  },
  filter = list(position = 'top', clear = FALSE)
  )

  output$downloadPlot_MASC <- downloadHandler(
    filename = "Shiny_MASC_plot.png",
    content = function(file) {
      png(file, height= 8, width = 9,res=1200, units = "in")
      print(MASC_plot())
      dev.off()
    }
  )

  output$downloadTable_MASC <- downloadHandler(
    filename = "Shiny_MASC_table.csv",
    content = function(file) {
      write.csv(MASC_table(), file)
    }
  )
  ################

  volc_hc_ep_plot <- function(){

    if(input$HC_vs_EP_Level=="L1"){
      hc_vs_ep_df <- hc_vs_ep_df_l1
    } else{
      hc_vs_ep_df <- hc_vs_ep_df_l2
    }

    df <- hc_vs_ep_df %>%
      filter(cell_type %in% sub_types_epil_ctrl()) %>%
      mutate(differential = "no change") %>%
      mutate(differential = ifelse(p_adj <= 0.05 & logFC >0, "Up-regulated",differential)) %>%
      mutate(differential = ifelse(p_adj <= 0.05 & logFC < 0,"Down-regulated",differential)) %>%
      filter(differential != "no change")
    plt <- hc_vs_ep_df %>%
      filter(cell_type %in% sub_types_epil_ctrl()) %>%
      filter(!(gene %in% df$gene)) %>%
      mutate(differential = "no change") %>%
      ggplot(aes(x=logFC, y= -log10(p_adj), color = differential))+
      geom_point()+
      geom_point(data =df, aes(x= logFC, y= -log10(p_adj), color = differential)) +
      ggrepel::geom_label_repel(data= df %>% top_n(15,abs(logCPM)), aes(label = gene),
                                show.legend = F,max.overlaps=100)+
      geom_hline(yintercept=-log10(0.05), linetype="dashed",color ="gray")+
      theme_bw()+
      scale_color_manual(values=c("Up-regulated"="darkred", "Down-regulated"="blue","no change"="black"))+
      ggtitle(sub_types_epil_ctrl())+
      theme(legend.position = "down")
    return(plt)
  }

  vln_hc_ep_plot <- function(){

    if(input$HC_vs_EP_Level=="L1"){
      genes <- hc_vs_ep_df_l1 %>%
        filter(cell_type %in% sub_types_epil_ctrl(), p_adj <= 0.05) %>%
        top_n(12, abs(logFC)) %>% pull(gene)
      vln_plt <- VlnPlot(subset(pbmc,subset = predicted.celltype.l1_NO_SPACES %in% sub_types_epil_ctrl()),
                         features = genes, ncol=as.numeric(input$Vln_Cols), group.by = "Condition")
    } else{
      genes <- hc_vs_ep_df_l2 %>%
        filter(cell_type %in% sub_types_epil_ctrl(), p_adj <= 0.05) %>%
        top_n(15, abs(logCPM)) %>% arrange(desc(logCPM)) %>% pull(gene)

      vln_plt <- VlnPlot(subset(pbmc,subset = predicted.celltype.l2_NO_SPACES %in% sub_types_epil_ctrl()),
                         features = genes, ncol=as.numeric(input$Vln_Cols), group.by = "Condition")
    }
    return(vln_plt)
  }
  output$Volcano_HC_EP <- renderPlot({
    volc_hc_ep_plot()
  })
  output$Vln_Plots_HC_EP <- renderPlot({
    vln_hc_ep_plot()
  })

  output$downloadPlots_HC_vs_EP <- downloadHandler(
    filename = "Shiny_HC_vs_EP_plot.csv",
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(volc_hc_ep_plot()+ vln_hc_ep_plot())
      dev.off()
    }
  )
  ###### TAU######
  output$tau_text <- renderText({"Works only for Cell Type categories."})

  Tau_plot <- function(){

    if(input$Variable_tau=="Cell Type (L1)"){
      data_input <- l1_l2 %>% filter(is.na(L1))
    } else {
      data_input <- l1_l2 %>% filter(!is.na(L1))
    }

    if(!("ALL" %in% str_remove_all(sub_types_tau()," "))){
      if(is.null(sub_types_tau())){
        data_input <- data_input
      } else{
        data_input <- data_input %>%
          filter(cell_type %in% str_remove_all(sub_types_tau()," "))
      }
    }
    # plt_1 <- data_input %>%
    #   ggplot(aes(x=Condition, y= value, color = cell_type, group = cell_type))+
    #   geom_line(aes(linetype = delta_p_sig,group=groups))+
    #   #  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
    #   ggrepel::geom_label_repel(data = data_input %>% filter(Condition=="MC") %>%
    #                               mutate(cell_type =str_remove(cell_type,"_prop")),
    #                             aes(label = cell_type, color=cell_type),
    #                             size=3,min.segment.length = 1, show.legend = F)+
    #   geom_hline(yintercept = 0)+
    #   geom_point(aes(shape = value_p_sig))+
    #   scale_linetype(name = "Significant Change")+
    #   theme_bw()+
    #   ylab("Association with tau level changes")+
    #   xlab("Group")+
    #   guides(color = "none")+
    #   theme(legend.position = "bottom")
    #
    plt_1 <- data_input %>%
      ggplot(aes(x=Condition, y=value,group=Condition))+
      geom_bar(aes(color=Condition, fill=Condition),position = position_dodge(width = 1),stat="identity")+
      geom_point(position = position_dodge(width=1), color = "black") +
      geom_errorbar(aes(ymax=upper, ymin=lower),width=.2, color = "black", position = position_dodge(width=1))+
      geom_hline(yintercept =0, color = "black", linetype="longdash")+
      facet_grid(~tau)+
      theme_bw()+
      ylab("Tau associations with cell type proportion")+
      scale_color_brewer(palette="Set1")+
      scale_fill_brewer(palette="Set1")+
      facet_grid(tau~cell_type,scales ="free_x")+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position="bottom")

    if(nrow(data_input %>% filter(change_sig==0))>=0){
      plt_1 <- plt_1 + add_pvalue(data=data_input %>% filter(change_sig==0) %>% distinct(cell_type,tau, delta_p) %>%
                                    mutate(ConditionA="MC",
                                           ConditionB="MR",
                                           y.position=1,
                                           label = round(delta_p,4)),
                                  xmin="ConditionA",xmax="ConditionB",label="label",y.position="y.position")
    }

    return(plt_1)
  }

  output$Tau <-renderPlot({
    Tau_plot()
  }, height=500,width= 600)

  output$downloadPlot_Tau <- downloadHandler(
    filename = "Shiny_Tau_plot.png",
    content = function(file) {
      png(file, height= 8, width = 9,res=1200, units = "in")
      print(Tau_plot())
      dev.off()
    }
  )






  ###############
  observe({
    ft_df <- summarized_data %>% filter(exclude == input$L1_Pairwise,
                                        direction ==input$Direction_L1,
                                        level =="L1") %>% pull(cell_type)
    print(ft_df)
    updateRadioButtons(session, "DEG_L1_choice",
                       label = "Cell type",
                      choices = sort(unique(ft_df)),
                      selected =  sort(unique(ft_df))[1])
  })

  pair_degs_L1 <- reactive({
    input$DEG_L1_choice
  })

  observe({
    ft_df_l2 <- summarized_data %>% filter(exclude==input$L2_Pairwise,
                                           direction==input$Direction_L2,
                                           level =="L2") %>% pull(cell_type)
    print(ft_df_l2)
    updateRadioButtons(session, "DEG_L2_choice",
                       label = "Cell sub-type",
                      choices = sort(unique(ft_df_l2)),
                      selected =  sort(unique(ft_df_l2))[1])
  })

  pair_degs_L2 <- reactive({
    input$DEG_L2_choice
  })


  deg_violins_l1 <- function(){
    features_input <- features_df(input$L1_Pairwise, input$Direction_L1,"L1") %>%
      filter(cell_type == str_remove_all(pair_degs_L1()," ")) %>% arrange(desc(logCPM)) %>% pull(gene)
    vln <- VlnPlot(subset(pbmc, predicted.celltype.l1_NO_SPACES==pair_degs_L1()),
                   features = features_input, ncol=5, group.by ="Condition")
    return(vln)
  }

  deg_violins_l2 <- function(){
    features_input <- features_df(input$L2_Pairwise, input$Direction_L2,"L2") %>%
      filter(cell_type == str_remove_all(pair_degs_L2()," ")) %>% arrange(desc(logCPM)) %>% pull(gene)
    vln <- VlnPlot(subset(pbmc, predicted.celltype.l2_NO_SPACES==pair_degs_L2()),
                   features = features_input, ncol=5, group.by="Condition")
    return(vln)
  }
  deg_l1_table_pair <- function(){
    tab <- features_df(input$L1_Pairwise, input$Direction_L1,"L1")%>%
      filter(cell_type == str_remove_all(pair_degs_L1()," ")) %>%
      arrange(p_adj, desc(logCPM)) %>% dplyr::select(-any_of(c("cell_type","level","exclude","direction")))
    return(tab)
  }
  deg_l2_table_pair <- function(){
    tab <- features_df(input$L2_Pairwise, input$Direction_L2,"L2")%>%
      filter(cell_type == str_remove_all(pair_degs_L2()," ")) %>%
      arrange(p_adj, desc(logCPM)) %>% dplyr::select(-any_of(c("cell_type","level","exclude","direction")))
    return(tab)
  }


  output$DEG_Violins_L1 <- renderPlot({
    deg_violins_l1()
  })
  output$DEG_Violins_L2 <- renderPlot({
    deg_violins_l2()
  })

  output$DEG_L1_Table_Pair <- renderDataTable({
    deg_l1_table_pair()
  }, filter = list(position = 'top', clear = FALSE))

  output$DEG_L2_Table_Pair <- renderDataTable({
    deg_l2_table_pair()
  }, filter = list(position = 'top', clear = FALSE))

  output$downloadPlot_Pairwise_L1 <- downloadHandler(
    filename = paste0("Shiny_Vln_plots",pair_degs_L1(),"_",input$L1_Pairwise,"_",input$Direction_L1,".png"),
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(deg_violins_l1())
      dev.off()
    }
  )
  output$downloadPlot_Pairwise_L2 <- downloadHandler(
    filename = paste0("Shiny_Vln_plots",pair_degs_L2(),"_",input$L2_Pairwise,"_",input$Direction_L2,".png"),
    content = function(file) {
      png(file, height= 8, width = 11,res=1200, units = "in")
      print(deg_violins_l2())
      dev.off()
    }
  )

  ####################

  dim_plot_A <- function(){
    fp <-FeaturePlot(pbmc, features = input$marker_A)
    return(fp)
  }
  dim_plot_B <- function(){
    fp <-FeaturePlot(pbmc, features = input$marker_B)
    return(fp)
  }

  vln_plot_A <- function(){
    var <- switch(input$marker_A_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "gender",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "Has TCR 10x data" = "is_t_cell")
    pairwise.wilcox.test(GetAssay(pbmc, "SCT")[input$marker_A,]  %>% as.vector(),
                         pbmc@meta.data[[var]], "BH") -> res
    res$p.value %>% as.data.frame() %>% mutate(compA = rownames(res$p.value)) %>%
      pivot_longer(-compA,names_to = "comp") %>% filter(value <= .05)  %>% top_n(-6, value) -> lis_df
    lis <- lapply(1:nrow(lis_df), function(i) c(pull(lis_df[i,"comp"]),pull(lis_df[i,"compA"])))

    vln <- VlnPlot(pbmc, features = input$marker_A, group.by = var, y.max=5)+NoLegend()+
      stat_compare_means(comparisons = list(unique(pbmc@meta.data[[var]])), label = "p.signif",na.rm=T, hide.ns = T)
    return(vln)
  }
  vln_plot_B <- function(){
    var <- switch(input$marker_B_variable,
                  "Cell Type (L1)" = "predicted.celltype.l1",
                  "Cell Type (L2)" = "predicted.celltype.l2",
                  "Sex" = "gender",
                  "Condition" = "Condition",
                  "Seurat Clusters" = "seurat_clusters",
                  "Has TCR 10x data" = "is_t_cell")
    pairwise.wilcox.test(GetAssay(pbmc, "SCT")[input$marker_B,]  %>% as.vector(),
                         pbmc@meta.data[[var]], "BH") -> res
    res$p.value %>% as.data.frame() %>% mutate(compA = rownames(res$p.value)) %>%
      pivot_longer(-compA,names_to = "comp") %>% filter(value <= .05)  %>% top_n(-6, value) -> lis_df
    lis <- lapply(1:nrow(lis_df), function(i) c(pull(lis_df[i,"comp"]),pull(lis_df[i,"compA"])))

    vln <- VlnPlot(pbmc, features = input$marker_B, group.by = var, y.max=5)+ NoLegend()+
      stat_compare_means(comparisons = list(unique(pbmc@meta.data[[var]])), label = "p.signif",na.rm=T)
    return(vln)
  }

  output$DimPlot_A <- renderPlot({
    dim_plot_A()
  })
  output$DimPlot_B <- renderPlot({
    dim_plot_B()
  })

  output$Vln_Plot_A <- renderPlot({
    vln_plot_A()
  })
  output$Vln_Plot_B <- renderPlot({
    vln_plot_B()
  })

  output$dowloadPlot_Grid <- downloadHandler(
    filename = paste0(input$marker_A,"_",input$marker_B,"_grid.png"),
    content = function(file) {
      ggsave(file, (dim_plot_B()+dim_plot_B())/(vln_plot_A()+vln_plot_B()),
             height=12, width=12)
    }
  )


  }
