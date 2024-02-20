#' welcomePage UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_welcomePage_ui <- function(id){
    ns <- NS(id)
    tagList(
        fluidRow(
            div(
                style ="padding-left:3%; padding-right:3%; padding-bottom:3%; font-size:large",
                fluidRow(
                    style = "margin:2%; padding-left: 2%; padding-right: 2%;",
                    p(
                        HTML( "Welcome to <b> gbm.cregmap </b>"),
                        HTML("<i> (The co-regulatory influence map of GBM heterogeneity and evolution— From cells
                             in vitro to tumors and back again)</i>"),
                        HTML(", a powerful web-based tool to provide researchers with rapid access
                             to a unified coregulatory influence network view of GBM cancer heterogeneity and evolution. ")
                        ,
                        style = "font-size:150%;")
                ),
                fluidRow(
                    style = "margin:2%;",
                    ## I used tags$table to have equal height columns
                    tags$table(style = "width: 100%",
                               tags$tr(
                                   tags$td(
                                       style = "width: 70%",
                                       HTML("<center><img src='/img/figureWelcomePage.png' width='70%'/></center>")
                                   ),
                                   tags$td(
                                       style = "margin:2%;vertical-align: center;",
                                       HTML("
                                       <p>Developed in the DISCO team by
                                       <li>Geoffrey P. </li>
                                       <li>Aurelien D. </li>
                                       <li>Mohamed E. </li></p>
                                       <br/>
                                       <p>Contact information: <br/>
                                       Mohamed Elati <br/>
                                       ONCOLille Cancer Institute <br/>
                                       CNRS UMR 9020 - Inserm UMR 1277 CANTHER <br/>
                                       Team DISCO (DIgital, Systems and COmputational cancer) <br/>
                                       Email : <a href=mailto:“mohamed.elati@univ-lille.fr”>mohamed.elati@univ-lille.fr</a> <br/>
                                       Web : team DISCO website </p>"
                                       )
                                   )
                               )
                    )
                ),
                fluidRow(
                    style = "margin-top:4%;font-size:120%;",
                    HTML("<center>
                    <p>
                    Please acknowledge GBM-cRegMap in your publications by citing the following references:
                    </p>
                    <br/>
                    <p>
                    GBM-cRegMap:  <a href=''>XXX et al. XXXX</a> <br/>
                    CoRegNet: <a href='https://academic.oup.com/bioinformatics/article/31/18/3066/240755?login=true'>Nicolle et al. 2015</a> <br/>

                    Regulator Influence: <a href='https://ieeexplore.ieee.org/document/6406597'>Nicolle et al. 2012</a>
                    &
                    <a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2481-y'>Dhifli et al. 2019</a> <br/>

                    LICORN: <a href='https://ieeexplore.ieee.org/document/6803994'>Chebil et al. 2014</a>
                    &
                    <a href='https://academic.oup.com/nar/article/41/3/1406/2902980'>Elati et al. 2007</a>
                    </p>
                         </center>")
                )
            )
        )
        # )
    )
}

#' welcomePage Server Functions
#'
#' @noRd
mod_welcomePage_server <- function(id, r=r){
    moduleServer( id, function(input, output, session){
        ns <- session$ns



    })
}

## To be copied in the UI
# mod_welcomePage_ui("welcomePage_1")

## To be copied in the server
# mod_welcomePage_server("welcomePage_1", r=r)
