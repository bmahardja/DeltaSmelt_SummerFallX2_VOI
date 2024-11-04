#SampleFlowChart.R----

# Brian Mahardja
# Start Date: Nov 4, 2024

#About----
#Project:Delta Smelt Summer-fall X2 VOI study flow chart visualization 

#Libraries ----
library(tidyverse)
#library(treemap)
#library(data.tree)
library(DiagrammeR)
#library(igraph)
#library(networkD3)
#library(collapsibleTree)
#library(readr)
library(DiagrammeRsvg)
library(rsvg)
library(tiff)

# Set working directory ----
setwd("~/GitHub/DeltaSmelt_SummerFallX2_VOI/Figures/Workflow")

# Create a flow chart
flow_chart <- grViz("
digraph flowchart {
  node [fontname = Helvetica, shape = rectangle]
  
  Start [label = 'CalSim 3 (California Water Operation Model)']
  Decision1 [label = 'X2 location']
  EndPoint1 [label = 'Water Cost Calculation (million cubic meter)' style=filled fillcolor=white, color=red, penwidth=2]
  Process1a [label = 'X2-Salinity Submodel']
  Process1b [label = 'Salinity-Zooplankton Submodel']
  Process2 [label = 'Delta Smelt Distribution Submodel']
  Process3a [label = 'Delta Smelt IBMR v1']
  Process3b [label = 'Delta Smelt IBMR v2']
  EndPoint2 [label = 'Delta Smelt Population Growth (lambda)' style=filled fillcolor=white, color=red, penwidth=2]

  // Subgraph for IBMR
  subgraph cluster_ibmr {
    style=dotted;
    label = 'IBMR';
    color=black; // Border color of the dotted box
    Process1a;
    Process1b;
    Process2;
    Process3a;
    Process3b;
  }
  
  Start -> EndPoint1
  Start -> Decision1
  Decision1 -> Process1a
  Process1a -> Process1b
  Process1b -> Process3a [label = 'On/Off']
  Process1b -> Process3b [label = 'On/Off']
  Decision1 -> Process2
  Process2 -> Process3a [label = 'On/Off']
  Process2 -> Process3b [label = 'On/Off']
  Process3a -> EndPoint2
  Process3b -> EndPoint2
}
")

# Render the flow chart
flow_chart

# Convert the flow chart to SVG
svg_output <- DiagrammeRsvg::export_svg(flow_chart)

# Save the SVG output to a temporary file
temp_svg_file <- tempfile(fileext = ".svg")
writeLines(svg_output, temp_svg_file)

# Read the SVG and convert it to a PNG format
png_output <- rsvg::rsvg(temp_svg_file)

# Save the PNG output as a TIFF file with LZW compression
tiff("flowchart.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
grid::grid.raster(png_output)  # Use grid to render the raster image
dev.off()  # Close the graphics device

# Optionally, remove the temporary SVG file
unlink(temp_svg_file)