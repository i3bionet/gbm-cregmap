#' FAQ UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_FAQ_ui <- function(id){
    # `NS(id)` returns a namespace function, which was save as `ns` and will
    # invoke later.
    ns <- NS(id)
    tagList(
        uiOutput(ns("tutorialUI"))
    )
}

#' FAQ Server Functions
#'
#' @noRd
mod_FAQ_server <- function(id, r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns
        output$tutorialUI <- renderUI({
            div(
                fluidRow(style = "padding-left: 3%; padding-right: 3%; margin:1%",
                         p(
                             HTML( paste0("Welcome to <b> ", tolower(r$cregmapVersionUpper) , ".cregmap </b>")),
                             HTML("<i> (The co-regulatory influence map of ",r$cregmapVersionUpper, " heterogeneity and evolution— From cells
                             in vitro to tumors and back again)</i>")
                             ,
                             style = "font-size:150%;")
                         # style = "font-size:150%;"
                ),
                HTML(
                    paste0("The ",r$cregmapVersionUpper, "-cRegMap
                <a href='https://", tolower(r$cregmapVersionUpper) , ".cregmap.com/'>(https://",
                           tolower(r$cregmapVersionUpper) , ".cregmap.com/)</a>
                is an interactive online tool that provides researchers with swift access to a
                unified coregulatory influence network view of bladder cancer heterogeneity and plasticity.
                The tool includes three main tabs: “CoRegNet”, “CoRegMap”, “CoRegQuery” and
                two supplementary tabs: “Data sets” and this one, the “Tutorial”.
                <br/>
                <br/>
                &emsp;With the <b>“CoRegNet”</b> tab,the coregulatory network inferred from the transcriptome of cancer-derived
                cell lines is visualized. CoRegNet network representations are intuitively understandable and facilitate
               the integration of several layers of information over nodes and edges. Nodes represent transcription
               factors and cofactors (TFs/coTFs). Node color (red= high; blue = low) indicates the influence of the
               corresponding TF/coTF. Node size is proportional to the number of targets of the TF/coTF and the
               intensity of the color indicates the strength of its influence. The coregulatory interactions between
               nodes are indicated as follows: predicted co-operativity interactions are shown in gray, and interactions
               for which there is published evidence are shown in blue (for protein-protein interactions, ppi) or red
               (for transcriptional regulation, tfbs). This tab is designed to allow users to explore the coregulatory
               network and to identify core active coregulators in single sample (e.g., cell line) or particular phenotype
               (e.g., tumor subtype) for further analysis. The ",r$cregmapVersionUpper, "-CoRegNet network is visualized via a
                <a href='https://shiny.rstudio.com/'>Shiny</a> and
                <a href='https://cytoscape.org/'>Cytoscape</a> javascript applet
                in the R package <a href='https://datastorm-open.github.io/visNetwork/'>visNetwork</a>. <br>
                   <br/>
                   <br/>
                   <br/>
                   <h3>Interface Outline</h3>
                   <br/>
                   <h4>Data sets</h4>
                   <br/>
                   <center><img src='/img/Tutorial_imgs/Fig1_dataset.png' width='75%'/></center>
                   <br/>
                   <p>
                   In this tab the user can find tables that include the available data which
               have been used as input in the CoRegNet / CoRegMap tabs.
               The user can search the data tables to identify information of
               interest regarding them.
               </p>
                   <br/>
                   <h4>CoRegNet</h4>
                   <br/>
                   <center><img src='/img/Tutorial_imgs/Fig2_coregnet.png' width='75%'/></center>
                   <br/>
                   <p>
                   In this tab the user can find the visualized ",r$cregmapVersionUpper,
                           " Co-Regulatory Network.
               Several layers of information are being displayed over nodes and edges.
               Nodes represent transcription factors and cofactors (TFs/coTFs).
               <br/>
                   The color of the nodes
               (<b style='color:red;'>red</b> = high; <b style='color:blue;'>blue</b> = low) indicates the
               influence of the corresponding TF/coTF.
               Node size is proportional to the number of targets of
               the TF/coTF and the intensity of the color indicates the strength of its influence.
               The co-regulatory interactions between nodes are indicated as follows:
                   interactions defined solely by the CoRegNet algorithm are shown in gray,
               and interactions for which there is published evidence are shown
               in blue (for protein-protein interactions, ppi) or
               red (for transcriptional regulation, transcription factor binding sites; TFBs).
               <br/>
                   This tab is designed to allow users to explore
               the coregulatory network and to identify core active coregulators
               in a single sample (e.g., cell line) or particular phenotype
               (e.g., tumor subtype) for further analysis.
               </p>
                   <br/>
                   <p>
                   Intuitively the user can start from selecting the cohort from the options of:
                   </p>
                   <ul>
                   <li>Meta-cohort of ",r$cregmapVersionUpper, " tumors (initially selected)</li>
                   <li>CCLE ",r$cregmapVersionUpper, " cell lines, or</li>
                   </ul>
                   <p>
                   and then proceed choosing the annotation with either:
                   </p>
                   <ul>
                   <li>Verhaak (initially selected)</li>
                   <li>",r$cregmapVersionUpper, "-cRegMap</li>
                   <li>G-CIMP status</li>
                   </ul>
                   <p>
                   Next the user can choose to visualize the subtypes on the network
               with respect to the previous selection regarding the annotation.
               </p>
                   <p>
                   Inside the diagram the user can select to visualize only
               the gene id they are interested in at the specific moment.
               They can store the current network and then download it as HTML
               </p>
                <h4>Contact Information</h4>
                <br/>
                <p>
                For any questions, suggestions or additional information please contact:
                <b>mohamed.elati <at> univ-lille.fr</b>
                </p>
                <br/>
                ")),
                # style ="padding:10px 50px; text-align:justify; font-size:large"
                style ="padding-left:3%; padding-right: 3%; text-align:justify; font-size:large"
            )
        })
        browserComp <- data.frame(OS = c("Linux","MacOS","Windows"),
                                  Version = c("to fill","to fill","to fill"),
                                  Chrome = c("to fill", "to fill", "to fill"),
                                  Firefox = c("to fill", "to fill", "to fill"),
                                  'Microsoft Edge' = c("to fill", "to fill", "to fill"),
                                  Safari = c("to fill", "to fill", "to fill"),
                                  Opera = c("to fill", "to fill", "to fill"))
        output$browserTable <- formattable::renderFormattable({formattable::formattable(browserComp, list())})
    })
}

## To be copied in the UI
# mod_FAQ_ui("FAQ_1")

## To be copied in the server
# mod_FAQ_server("FAQ_1", r=r)
