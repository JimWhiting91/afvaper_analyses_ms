library(tidyverse)  # for %>% pipes
library(DiagrammeR)
library(DiagrammeRsvg)  # for conversion to svg
library(rsvg)  # for saving svg

afvaper_workflow <- grViz("
  digraph graph2 {

  graph [layout = dot, rankdir = TB]

  # node definitions with substituted label text
  node [shape = oval,color=grey,fontname=arial,style=filled]
  input_vcf [label = 'Single Chromosome VCF (vcfR object)']
  input_popmap [label = 'Popmap']
  input_vectors [label = 'Vector List']
  a [label = 'calc_AF_vectors()']
  b [label = 'calc_AF_vectors(null_perms)']
  c [label = 'eigen_analyse_vectors()']
  d [label = 'find_null_cutoffs()']
  e [label = 'eigen_pvals()']
  f [label = 'eigenval_plot()']
  g [label = 'signif_eigen_windows()']
  h [label = 'summarise_window_parallelism()']
  # i [label = 'merge_eigen_res()']
  # j [label = 'vcfR2AF()']
  # k [label = 'sum_eigenvals()']

  # Define edge defintions
  # a -> b -> c -> d
  input_vcf -> a;
  input_popmap -> a;
  input_vectors -> a;
  input_vcf -> b;
  input_popmap -> b;
  input_vectors -> b;
  a -> c;
  b -> d;
  c -> e;
  b -> e;
  c -> f;
  d -> f;
  d -> g;
  c -> g -> h;
  c -> h;
  }")

# Save to SVG + PNG
afvaper_workflow_svg <- export_svg(afvaper_workflow)
afvaper_workflow_svg <- charToRaw(afvaper_workflow_svg)
rsvg_pdf(afvaper_workflow_svg,"~/Exeter/angle_grinder/figs/FigureX_afvaper_workflow.pdf")
rsvg_png(afvaper_workflow_svg,"~/Exeter/angle_grinder/figs/FigureX_afvaper_workflow.png")

