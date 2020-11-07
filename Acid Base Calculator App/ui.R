#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# App calculates titration diagrams as well as logC pH and C pH diagrams
# for acids and bases with up to three dissociation stages.
# App programmed by Dr. Kai Exner 2020 / 11 - drkaiexner@gmail.com

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Acid-Base Calculator"),
   
   # Sidebar to collect user input 
   sidebarLayout(
     
      sidebarPanel(
        
        radioButtons("acidOrBase",
                     "Acid or base to be studied",
                     choices = c("acid", "base"),
                     selected = "acid",
                     inline = TRUE
          ),
        
        radioButtons("stages",
                     "no. of dissociation stages acid",
                     choices = c(1, 2, 3),
                     selected = 1,
                     inline = TRUE
          ),
        
        sliderInput("pK1",
                    "pK1",
                    min = -10,
                    max = 14,
                    value = 0
          ),

        conditionalPanel(condition = "input.stages > 1",
                         sliderInput(
                           "pK2",
                           "pK2",
                           min = -10,
                           max = 14,
                           value = 0
                         )
          ),

        conditionalPanel(condition = "input.stages == 3",
                         sliderInput(
                           "pK3",
                           "pK3",
                           min = -10,
                           max = 14,
                           value = 0
                         )
          ),

        sliderInput("C00",
                    "Conc. analyte solution [mol / L]",
                    min = 0.01,
                    max = 1,
                    value = 0.1
          ),
        
        sliderInput("CT00",
                    "Conc. titrant [mol / L]",
                    min = 0.01,
                    max = 1,
                    value = 0.1
          )
        
      ),
      
      # Tabs for individual plots and ReadMe
      mainPanel(
         tabsetPanel(type = "tabs",
           
           tabPanel("Titration Diagram", plotOutput("titrationDiagram")),
           
           tabPanel("logC pH Diagram", plotOutput("logCpHDiagram")),
           
           tabPanel("C pH Diagram", plotOutput("cpHDiagram")),
           
           tabPanel("ReadMe",
             
             HTML("<h3>App Description</h3><br>
               <p>This app calculates titration diagrams as well as logC pH and 
               C pH diagrams for acids and bases with up to three dissociation 
               stages. In the latter two diagrams the concentration of all 
               species present in solution is plotted against the pH of the 
               solution. All calculations done are exact in the sense that no 
               approximations are utilized to solve the underlying equations. 
               Yet, concentrations are used as such, no corrections for 
               activities are introduced.<br>
               </p>
               <p>User input:<br>
               - decision whether an acid or base should be studied  <br>
               - number of dissociation stages  <br>
               - pK values for individual dissociation stages  <br>
               - concentration of the acid resp. base solution  <br>
               - concentration of the base resp. acid to be used for the 
               titration  <br><br>
               </p>
               <p>The titrant is assumed to be a strong acid resp. base.<br></p>
               The code for the app as well as a tutorial 
               that leads through the background of acid-base calculations 
               incl. R code for all calculations and plots in the 
               tutorial are also available on GitHub.<br><br>
               For questions about this app and suggestions please contact 
               drkaiexner@gmail.com
               <p>
               
               </p>")
             )
           )
        )
     )
  )
