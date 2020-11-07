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

# Define server logic od App
server <- function(input, output) {
   
   output$titrationDiagram <- renderPlot({
     
     C00 = input$C00    # Initial concentration of analyte solution. Effect 
                        # of dilution during addition of titrant will be taken
                        # into account in calculations
     
     V0 = 1    # Volume of analyte solution set to 1. Volume of titrant given 
               # relative to volume of analyte solution.
     
     CT00 = input$CT00    # Concentration of titrant
     
     pK1 = input$pK1    # pK values for titration curve can be taken as such, no
                        # matter whether acid or base is titrated
     
     # If user choses acid / base with less than three dissociation stages,
     # pK values of stages "not in use" are set to infinity
     
     pK2 = ifelse(input$stages > 1, input$pK2, Inf) 
     
     pK3 = ifelse(input$stages == 3, input$pK3, Inf) 
     
     K1 = 10^-pK1
     
     K2 = 10^-pK2
     
     K3 = 10^-pK3
     
     KW = 1E-14
     
     pH = seq(from = 0 , to = 14, by = 0.01)
     
     H = 10^-pH
     
     OH = 10^-14 / H
     
     
     if (input$acidOrBase == "acid") {
       
       H_OH = H
       
       fact = 1    # Equations for titration of acids and bases only differ
                   # by sign of [H+] and Kw / [H+] = [OH-]. Variable fact
                   # therefore switches between cases.
       
     } else {
       
       H_OH = OH
       
       fact = -1
       
     }
     
     K1H = K1 / H_OH
     
     K12H2 = K1 * K2 / H_OH^2
     
     K123H3 = K1 * K2 * K3 / H_OH^3
     
     equivIons = (K1H + 2 * K12H2 + 3 * K123H3) * C00 / 
       (1 + K1H + K12H2 + K123H3)
     
     # Calculate volume of titrant that needs to be added to achieve given pH
     # Negative values will occur for pH values that are less (acids) or higher
     # (bases) than the initial pH of the analyte solution.
     # pH values beyond the pH of the titrant will also lead to unreasonable
     # volumes. Both cases will be handled by capping the VT values.
     
     VT = V0 / (CT00 + fact * H - fact * KW / H) * 
       (equivIons + fact * KW / H - fact * H)
     
     curveData = data.frame(pH, VT)
     
     nProt = 1
     
     if (K2 != 0) nProt = 2
     
     if (K2 * K3 != 0) nProt = 3
     
     # Define max volume of titrant based on dissociation stages of analyte.
     # Excess to correspond to one additional dissociation stage.
     
     VTMax = (nProt + 1) * (V0 * C00) / CT00
     
     # Cap unreasonalbe VT values
     
     curveData = curveData[(curveData$VT > 0 & curveData$VT <= VTMax), ]
     
     plot(curveData$VT, curveData$pH, "l", 
       xlab = "vol. titrant / vol. analyte solution", ylab = "pH")
     
     grid()
     
   })
   
   output$logCpHDiagram <- renderPlot({
     
     # All calculations done based on acids as species. pK values entered for
     # bases are converted into pK values for corresponding acids. Code set
     # up for up to three dissociation stages. Handling of fewer than three
     # dissociation stages by setting "unused stages" to pK = infinity
     
     par(xpd = NA)    # Prevent clipping of plot at axis margins
     
     Kw = 1E-14
     
     pH = seq(from = 0, to = 14, length.out = 450)
     
     H = 10^-pH
     
     OH = Kw / H
     
     C0 = input$C00
     
     if (input$acidOrBase == "acid") {
       
        pK1 = input$pK1
        
        pK2 = ifelse(input$stages > 1, input$pK2, Inf)
        
        pK3 = ifelse(input$stages == 3, input$pK3, Inf)
     
       } else {
         
         # Transform pK_B into pK_A values of corresponding acids.
         # Order of must be changed: strongest base has weakest corresponding
         # acid.
         
        if (input$stages == 1) {
          
          pK1 = 14 - input$pK1
          
          pK2 = Inf
          
          pK3 = Inf
          
        }
      
        if (input$stages == 2) {
          
          pK1 = 14 - input$pK2
          
          pK2 = 14 - input$pK1
          
          pK3 = Inf
          
        }
         
        if (input$stages == 3) {
          
          pK1 = 14 - input$pK3
          
          pK2 = 14 - input$pK2
          
          pK3 = 14 - input$pK1
          
         }
         
       }
     
     K1 = 10 ^ -pK1
     
     K2 = 10 ^ -pK2
     
     K3 = 10 ^ -pK3
     
     
     H3A = C0 / (1 + K1 / H + K1 * K2 / H^2 + K1 * K2 * K3 / H^3)
     
     H2A = K1 / H * H3A
     
     HA = K1 * K2 / H^2 * H3A
     
     A = K1 * K2 * K3 / H^3 * H3A
     
     ymin = ifelse(log10(C0) > -16, -16, -20)
     
     plot(pH, log10(H), "l", xlab = "", ylab = "", 
       xlim = c(0, 14), ylim = c(ymin, 1), axes  = FALSE)
     
     pHMarker = c(0:14)
     
     cMarker = c(ymin:1)
     
     # Create markers "by hand" rather than by using grid().
     
     for (i in seq_along(cMarker)) {
       
       lines(c(0, 14), c(cMarker[i], cMarker[i]), col = "grey")
       
     }
     
     for (i in seq_along(pHMarker)) {
       
       lines(c(pHMarker[i], pHMarker[i]), c(ymin, 1), col = "grey")
       
     }
     
     # Mark pH = pK (buffer region) and pH = 1/2 (pK1 + pK2) (pH at equivalence)
     # point for species with more than one dissociation stages. 
     
     if (pK1 >= 0) lines(c(pK1, pK1), c(ymin, 1), col = "blue", lwd = 1)
     
     lines(c(pK2, pK2), c(ymin, 1), col = "blue", lwd = 1)
     
     lines(c(pK3, pK3), c(ymin, 1), col = "blue", lwd = 1)
     
     lines(c(0.5 * (pK1 + pK2), 0.5 * (pK1 + pK2)), c(ymin, 1), 
       col = "blue", lwd = 1, lty = "dotted")
     
     lines(c(0.5 * (pK2 + pK3), 0.5 * (pK2 + pK3)), c(ymin, 1), 
       col = "blue", lwd = 1, lty = "dotted")
     
     lines(pH, log10(H), "l")
     
     lines(pH, log10(OH), "l")
     
     
     selectH3A = which(H3A > 10^ymin)
     
     selectH2A = which(H2A > 10^ymin)
     
     selectHA = which(HA > 10^ymin)
     
     selectA = which(A > 10^ymin)
     
     lines(pH[selectH3A], log10(H3A[selectH3A]), col = "red", lwd = 2)
     
     lines(pH[selectH2A], log10(H2A[selectH2A]), 
       col = "red", lwd = 2, lty = "dashed")
     
     lines(pH[selectHA], log10(HA[selectHA]), 
       col = "red", lwd = 2, lty = "dotdash")
     
     lines(pH[selectA], log10(A[selectA]), 
       col = "red", lwd = 2, lty = "dotted")
     
     axis(2)
     
     mtext("log(C)", side = 2, line = 3)
     
     axis(3)
     
     mtext("pH", side = 3, line = 3)
     
     box()
     
     legendTextAcid = list(
       
       c(expression("["*"HA"*"]"), expression("["*A^"-"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*H[2]*A*"]"), expression("["*H*A^"-"*"]"), 
         expression("["*A^"2-"*"]"), 
         
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*H[3]*A*"]"), expression("["*H[2]*A^"-"*"]"), 
         expression("["*HA^"2-"*"]"), expression("["*A^"3-"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]"))
       
       )
     
     legendTextBase = list(
       
       c(expression("["*BH^"+"*"]"), expression("["*"B"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*BH[2]^{2*"+"}*"]"), expression("["*BH^"+"*"]"), 
         expression("["*"B"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*BH[3]^{3*"+"}*"]"), expression("["*BH[2]^{2*"+"}*"]"), 
         expression("["*BH^"+"*"]"), expression("["*"B"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]"))
       
     )
     
     legendLwd = c(2, 2, 2, 2, 2)
     
     legendCol = list(
       
       c("red", "red", "black"),
       
       c("red", "red", "red", "black"),
       
       c("red", "red", "red", "red", "black")
       
     )
     
     legendLty = list(
       
       c("solid", "dashed", "solid"),
       
       c("solid", "dashed", "dotdash", "solid"),
       
       c("solid", "dashed", "dotdash", "dotted", "solid")
       
     )
     
     if (pK2 == Inf & pK3 == Inf) legendLength = 2
     
     if (pK2 != Inf & pK3 == Inf) legendLength = 3
     
     if (pK2 != Inf & pK3 != Inf) legendLength = 4
     
     if (input$acidOrBase == "acid") legendText = legendTextAcid 
     else legendText = legendTextBase
     
     legend(0, ymin - 2, 
       legend = legendText[[legendLength - 1]], 
       lwd = legendLwd[1:legendLength], 
       col = legendCol[[legendLength - 1]],
       lty = legendLty[[legendLength - 1]], 
       bty = "n",
       horiz = TRUE)
   })
   
   output$cpHDiagram <- renderPlot({
     
     # All calculations done based on acids as species. pK values entered for
     # bases are converted into pK values for corresponding acids. Code set
     # up for up to three dissociation stages. Handling of fewer than three
     # dissociation stages by setting "unused stages" to pK = infinity
     
     par(xpd = NA)    # Prevent clipping of plot region at axis margins
     
     Kw = 1E-14
     
     pH = seq(from = 0, to = 14, length.out = 900)
     
     H = 10^-pH
     
     OH = Kw / H
     
     C0 = input$C00
     
     if (input$acidOrBase == "acid") {
       
       pK1 = input$pK1
       
       pK2 = ifelse(input$stages > 1, input$pK2, Inf)
       
       pK3 = ifelse(input$stages == 3, input$pK3, Inf)
       
     } else {
       
       # Transform pK_B into pK_A values of corresponding acids.
       # Order of must be changed: strongest base has weakest corresponding
       # acid.
       
       if (input$stages == 1) {
         
         pK1 = 14 - input$pK1
         
         pK2 = Inf
         
         pK3 = Inf
         
       }
       
       if (input$stages == 2) {
         
         pK1 = 14 - input$pK2
         
         pK2 = 14 - input$pK1
         
         pK3 = Inf
         
       }
       
       if (input$stages == 3) {
         
         pK1 = 14 - input$pK3
         
         pK2 = 14 - input$pK2
         
         pK3 = 14 - input$pK1
       }
       
     }
    
     
     K1 = 10 ^ -pK1
     
     K2 = 10 ^ -pK2
     
     K3 = 10 ^ -pK3
     
     
     H3A = C0 / (1 + K1 / H + K1 * K2 / H^2 + K1 * K2 * K3 / H^3)
     
     H2A = K1 / H * H3A
     
     HA = K1 * K2 / H^2 * H3A
     
     A = K1 * K2 * K3 / H^3 * H3A
     
     
     selectH = which(H <= C0)
     
     selectOH = which(OH <= C0)
     
     plot(pH[selectH], H[selectH], "l", xlab = "", ylab = "", 
       xlim = c(0, 14), ylim = c(0, C0), axes  = FALSE)
     
     pHMarker = c(0:14)
     
     ymin = 0
     
     cMarker = seq(from = ymin, to = C0, length.out = 10)
     
     # Create markers "by hand" rather than by using grid().
     
     for (i in seq_along(cMarker)) {
       
       lines(c(0, 14), c(cMarker[i], cMarker[i]), col = "grey")
       
     }
     
     for (i in seq_along(pHMarker)) {
       
       lines(c(pHMarker[i], pHMarker[i]), c(ymin, C0), col = "grey")
       
     }
     
     # Mark pH = pK (buffer region) and pH = 1/2 (pK1 + pK2) (pH at equivalence)
     # point for species with more than one dissociation stages.      
     
     if (pK1 >= 0) lines(c(pK1, pK1), c(ymin, C0), col = "blue", lwd = 1)
     
     lines(c(pK2, pK2), c(ymin, C0), col = "blue", lwd = 1)
     
     lines(c(pK3, pK3), c(ymin, C0), col = "blue", lwd = 1)
     
     lines(c(0.5 * (pK1 + pK2), 0.5 * (pK1 + pK2)), c(ymin, C0), 
       col = "blue", lwd = 1, lty = "dotted")
     
     lines(c(0.5 * (pK2 + pK3), 0.5 * (pK2 + pK3)), c(ymin, C0), 
       col = "blue", lwd = 1, lty = "dotted")
     
     lines(pH[selectH], H[selectH], "l")
     
     lines(pH[selectOH], OH[selectOH], "l")
     
     lines(pH, H3A, col = "red", lwd = 2)
     
     lines(pH, H2A, col = "red", lwd = 2, lty = "dashed")
     
     lines(pH, HA, col = "red", lwd = 2, lty = "dotdash")
     
     lines(pH, A, col = "red", lwd = 2, lty = "dotted")
     
     axis(2)
     
     mtext("C", side = 2, line = 3)
     
     axis(3)
     
     mtext("pH", side = 3, line = 3)
     
     box()
     
     legendTextAcid = list(
       c(expression("["*"HA"*"]"), expression("["*A^"-"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*H[2]*A*"]"), expression("["*H*A^"-"*"]"), 
         expression("["*A^"2-"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*H[3]*A*"]"), expression("["*H[2]*A^"-"*"]"), 
         expression("["*HA^"2-"*"]"), expression("["*A^"3-"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]"))
     )
     
     legendTextBase = list(
       
       c(expression("["*BH^"+"*"]"), expression("["*"B"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*BH[2]^{2*"+"}*"]"), expression("["*BH^"+"*"]"), 
         expression("["*"B"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]")),
       
       c(expression("["*BH[3]^{3*"+"}*"]"), expression("["*BH[2]^{2*"+"}*"]"), 
         expression("["*BH^"+"*"]"), expression("["*"B"*"]"), 
         expression("["*H^"+"*"]"*" & "*"["*OH^"-"*"]"))
       
     )
     
     legendLwd = c(2, 2, 2, 2, 2)
     
     legendCol = list(
       
       c("red", "red", "black"),
       
       c("red", "red", "red", "black"),
       
       c("red", "red", "red", "red", "black")
       
     )
     
     legendLty = list(
       
       c("solid", "dashed", "solid"),
       
       c("solid", "dashed", "dotdash", "solid"),
       
       c("solid", "dashed", "dotdash", "dotted", "solid")
       
     )
     
     if (pK2 == Inf & pK3 == Inf) legendLength = 2
     
     if (pK2 != Inf & pK3 == Inf) legendLength = 3
     
     if (pK2 != Inf & pK3 != Inf) legendLength = 4
     
     if (input$acidOrBase == "acid") legendText = legendTextAcid 
       else legendText = legendTextBase
     
     legend(0, 0 - C0 / 10, 
       legend = legendText[[legendLength - 1]], 
       lwd = legendLwd[1:legendLength], 
       col = legendCol[[legendLength - 1]],
       lty = legendLty[[legendLength - 1]], 
       bty = "n",
       horiz = TRUE)
     
   })
}
