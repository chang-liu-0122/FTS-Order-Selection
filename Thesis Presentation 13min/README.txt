Order Selection and Forecasting in Long-Memory Functional Time Series
13-minute thesis presentation

Main files
  main.tex                         LaTeX Beamer source
  main_handout.tex                 One-page-per-logical-slide handout wrapper
  CNU13.sty                        CNU theme adapted from the supplied decks
  Order_Selection_Thesis_13min.pdf Final presentation PDF
  Order_Selection_Thesis_13min_handout.pdf
                                   Final handout PDF
  TALK_TIMING.txt                  13-minute pacing guide

Build
  latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build_latest main.tex
  latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build_latest_handout main_handout.tex

The deck contains 23 main slides plus 16 backup slides. The animated presentation
has 82 PDF pages; the handout has 39 logical slides. The appendixnumberbeamer
package keeps the main-slide counter at 23. Suggested timings and concise talk cues
are stored as comments and hidden Beamer notes in main.tex.

All figures in assets/ were copied from the thesis project and the two supplied
reference presentation folders, so the project is self-contained.
