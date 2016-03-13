library(shiny)
library(shinyRGL)

shinyUI(fluidPage(
  titlePanel("Cardinality data fitted to multiplicative model"),
  sidebarPanel( "Select Time Point",width = 2, radioButtons("radio", label="",c("T1","T2"))),
  mainPanel(
    "Click and drag to spin plot. Zoom in and out with two-finger swipe (or equivalent mouse action).",
    "Green points indicate CP-Knowers and red points indicate non-CP-Knowers. Note how the prediction plane flexes less at T2 compared to T1. This shows the diminished influence of OTS at the later time point.",
    "See ",a("mult_appendix",href="")," for model details.",
    webGLOutput("threedplot",width=800,height=800)
  )
))